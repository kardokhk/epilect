#!/usr/bin/env python3
"""Epilect command-line interface.

Subcommands
-----------
  classify    extract per-read methylation calls from a BAM (via modkit) and
              run the 2-strain likelihood-ratio classifier.
  phase       methylation-only or joint (methylation + SNP) per-read/allele
              phasing on a BAM.
  resolvable  tier conspecific strain pairs by methylation resolvability,
              either from a REBASE-style motif table or from a single pair of
              motif sets.

External tools
--------------
modkit is located via ``$MODKIT`` if set, otherwise the bundled symlink
``envs/bin/modkit`` relative to the repo, otherwise ``modkit`` on ``$PATH``.
samtools/bcftools/whatshap are taken from ``$PATH`` (the conda env).
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, Optional, Tuple

from . import __version__
from .core import classify as core_classify
from .core import density, pair_tier

# --------------------------------------------------------------------------- #
# Tool discovery
# --------------------------------------------------------------------------- #
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _find_modkit() -> str:
    cand = os.environ.get("MODKIT")
    if cand and (os.path.isabs(cand) or shutil.which(cand)):
        return cand
    bundled = os.path.join(_REPO_ROOT, "envs", "bin", "modkit")
    if os.path.exists(bundled):
        return bundled
    found = shutil.which("modkit")
    if found:
        return found
    raise SystemExit(
        "ERROR: modkit not found. Set $MODKIT to the binary path, install it "
        "(see environment.yml / README), or place it at envs/bin/modkit."
    )


def _need(tool: str) -> str:
    found = shutil.which(tool)
    if not found:
        raise SystemExit(f"ERROR: required tool '{tool}' not found on $PATH.")
    return found


# --------------------------------------------------------------------------- #
# Shared: parse modkit `extract calls` TSV into per-read channel counts
# --------------------------------------------------------------------------- #
# 0-based columns in modkit extract calls TSV:
#   0 read_id, 11 read_length, 13 call_code, 17 canonical_base, 19 fail, 23 motifs
_COL_READ = 0
_COL_RLEN = 11
_COL_CODE = 13
_COL_CANON = 17
_COL_FAIL = 19


def _parse_motif_spec(spec: str):
    """Parse a ``--motifs`` value like 'GATC:1:A:a,CCWGG:0:C:m'.

    Each entry is MOTIF:OFFSET:CANON_BASE:METH_CODE. Returns
      (modkit_args, channel_defs)
    where channel_defs maps motif -> (canon_base, meth_code).
    """
    modkit_args = []
    channels: Dict[str, Tuple[str, str]] = {}
    for entry in spec.split(","):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(":")
        if len(parts) != 4:
            raise SystemExit(
                f"ERROR: bad --motifs entry '{entry}'. "
                "Expected MOTIF:OFFSET:CANON_BASE:METH_CODE "
                "(e.g. GATC:1:A:a,CCWGG:0:C:m)."
            )
        motif, offset, canon, code = parts
        modkit_args += ["--motif", motif, offset]
        channels[motif] = (canon.upper(), code)
    if not channels:
        raise SystemExit("ERROR: --motifs is empty.")
    return modkit_args, channels


def _run_extract(modkit: str, bam: str, reference: str, motif_args, out_tsv: str,
                 threads: int) -> None:
    cmd = [modkit, "extract", "calls", "--reference", reference,
           "--preload-references", *motif_args, "-t", str(threads), bam, out_tsv]
    sys.stderr.write("[epilect] " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True)


def _aggregate_calls(tsv: str, channels: Dict[str, Tuple[str, str]]):
    """read_id -> {motif: (n_sites, n_methylated)} from a modkit calls TSV.

    A call is attributed to a channel by its canonical base; methylation is
    judged by whether the call_code matches the channel's methylation code.
    Failed (low-confidence) calls are skipped, matching the research scripts.
    """
    # canonical_base -> (motif, meth_code)
    by_canon = {canon: (motif, code) for motif, (canon, code) in channels.items()}
    agg: Dict[str, Dict[str, list]] = {}
    rlen: Dict[str, int] = {}
    with open(tsv) as fh:
        next(fh, None)  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= _COL_CANON:
                continue
            if f[_COL_FAIL] == "true":
                continue
            chan = by_canon.get(f[_COL_CANON])
            if chan is None:
                continue
            motif, meth_code = chan
            rid = f[_COL_READ]
            rec = agg.setdefault(rid, {})
            n, m = rec.get(motif, (0, 0))
            rec[motif] = (n + 1, m + (1 if f[_COL_CODE] == meth_code else 0))
            if rid not in rlen:
                try:
                    rlen[rid] = int(f[_COL_RLEN])
                except (ValueError, IndexError):
                    rlen[rid] = 0
    return agg, rlen


# --------------------------------------------------------------------------- #
# `classify`
# --------------------------------------------------------------------------- #
# Default SDSE GATC/CCWGG 2-strain emission model (per_read_classify.py:11-14).
DEFAULT_EMISSION = {
    "UT9728": {"GATC": 0.903, "CCWGG": 0.005},
    "UT10237": {"GATC": 0.041, "CCWGG": 0.873},
}
DEFAULT_MOTIFS = "GATC:1:A:a,CCWGG:0:C:m"


def _load_emission(path: str) -> Dict[str, Dict[str, float]]:
    with open(path) as fh:
        em = json.load(fh)
    if not isinstance(em, dict) or len(em) != 2:
        raise SystemExit(
            "ERROR: emission JSON must map exactly two strain names to "
            '{motif: P(methylated)} dicts, e.g. '
            '{"A": {"GATC": 0.9}, "B": {"GATC": 0.04}}.'
        )
    return em


def _estimate_emission(agg, channels) -> Dict[str, Dict[str, float]]:
    """2-cluster emission estimate when no model is supplied.

    Computes a per-read methylation fraction over all channels, splits reads at
    the median into a low- and high-methylation cluster, and uses each
    cluster's per-channel methylated fraction as the emission for a synthetic
    strain. This is a deterministic, unsupervised fallback so `classify` runs
    end-to-end without a known model; it is NOT the validated SDSE model.
    """
    motifs = list(channels.keys())
    fracs = []
    for rid, rec in agg.items():
        n = sum(rec.get(m, (0, 0))[0] for m in motifs)
        meth = sum(rec.get(m, (0, 0))[1] for m in motifs)
        if n:
            fracs.append((meth / n, rid))
    if not fracs:
        raise SystemExit("ERROR: no informative reads to estimate emission from.")
    fracs.sort()
    mid = len(fracs) // 2
    lo_ids = {rid for _, rid in fracs[:mid]}
    hi_ids = {rid for _, rid in fracs[mid:]}

    def cluster_emission(ids):
        out = {}
        for m in motifs:
            n = sum(agg[r].get(m, (0, 0))[0] for r in ids)
            meth = sum(agg[r].get(m, (0, 0))[1] for r in ids)
            out[m] = (meth / n) if n else 0.0
        return out

    return {"cluster_lo": cluster_emission(lo_ids),
            "cluster_hi": cluster_emission(hi_ids)}


def cmd_classify(args) -> int:
    modkit = _find_modkit()
    motif_args, channels = _parse_motif_spec(args.motifs)

    workdir = tempfile.mkdtemp(prefix="epilect_classify_")
    calls_tsv = os.path.join(workdir, "calls.tsv")
    try:
        _run_extract(modkit, args.bam, args.reference, motif_args, calls_tsv,
                     args.threads)
        agg, _rlen = _aggregate_calls(calls_tsv, channels)
    finally:
        if not args.keep_calls:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            sys.stderr.write(f"[epilect] kept calls TSV: {calls_tsv}\n")

    if args.estimate:
        emission = _estimate_emission(agg, channels)
        sys.stderr.write("[epilect] estimated 2-cluster emission model:\n")
        sys.stderr.write(json.dumps(emission, indent=2) + "\n")
    else:
        emission = _load_emission(args.emission)

    labels = list(emission.keys())
    la, lb = labels[0], labels[1]
    ea, eb = emission[la], emission[lb]

    rows = []
    counts = collections.Counter()
    for rid, rec in agg.items():
        n_info = sum(rec.get(m, (0, 0))[0] for m in channels)
        if n_info == 0:
            label, llr = "unclassifiable", 0.0
        else:
            label, llr = core_classify(rec, ea, eb, la, lb)
        rows.append((rid, label, llr, n_info))
        counts[label] += 1

    out = args.out
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["read_id", "strain", "llr", "n_informative_sites"])
        for rid, label, llr, n in rows:
            w.writerow([rid, label, f"{llr:.4f}", n])

    tot = len(rows)
    n_cls = sum(v for k, v in counts.items() if k != "unclassifiable")
    print(f"# wrote {out}")
    print(f"# strains: {la} vs {lb}  (model: "
          f"{'estimated' if args.estimate else args.emission})")
    print(f"# reads={tot:,}  classifiable={n_cls:,} "
          f"({100*n_cls/tot if tot else 0:.1f}%)  "
          f"unclassifiable={counts['unclassifiable']:,}")
    for k in labels:
        print(f"#   {k}: {counts[k]:,}")
    return 0


# --------------------------------------------------------------------------- #
# `phase`
# --------------------------------------------------------------------------- #
def _read_labels(path: str) -> Dict[str, str]:
    """Load read_id -> strain from a TSV (output of `classify` or a labels file).

    Accepts either the `classify` output (read_id, strain, ...) or a 2-column
    read_id/strain TSV. Skips 'unclassifiable' rows.
    """
    out = {}
    with open(path) as fh:
        first = fh.readline()
        has_header = "read_id" in first.lower()
        if not has_header:
            fh.seek(0)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            rid, strain = f[0], f[1]
            if strain in ("unclassifiable", "none", ""):
                continue
            out[rid] = strain
    return out


def _meth_only_phase(bam, vcf, labels, out):
    """Methylation-only allele phasing (meth_only_phasing.py core).

    For each het site, group covering reads by their methylation label, take
    each strain's majority base, and infer which strain carries ALT. Reports
    how many het sites were phasable. No truth needed.
    """
    import gzip
    import pysam

    # het sites: 0-based pos -> (REF, ALT)
    het = {}
    opener = gzip.open if vcf.endswith(".gz") else open
    with opener(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            het[int(f[1]) - 1] = (f[3], f[4])
    if not het:
        raise SystemExit("ERROR: no variants parsed from VCF.")

    strains = sorted(set(labels.values()))
    bamf = pysam.AlignmentFile(bam, "rb")
    contig = bamf.references[0]
    lo, hi = min(het), max(het)
    phased = []
    n_sites = 0
    for col in bamf.pileup(contig, lo, hi + 1, truncate=True,
                           min_base_quality=0, stepper="samtools"):
        pos = col.reference_pos
        hr = het.get(pos)
        if hr is None:
            continue
        n_sites += 1
        refb, altb = hr
        tally = {s: collections.Counter() for s in strains}
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            s = labels.get(pr.alignment.query_name)
            if s is None or s not in tally:
                continue
            base = pr.alignment.query_sequence[pr.query_position]
            tally[s][base] += 1
        # which strain's majority base is ALT?
        alt_strain = None
        majorities = {}
        for s in strains:
            if sum(tally[s].values()) == 0:
                majorities[s] = None
            else:
                majorities[s] = tally[s].most_common(1)[0][0]
        present = [s for s in strains if majorities[s] is not None]
        if len(present) < 2:
            continue
        carriers = [s for s in present if majorities[s] == altb]
        if len(carriers) != 1:
            continue  # ambiguous (none or all carry ALT)
        alt_strain = carriers[0]
        phased.append((contig, pos + 1, refb, altb, alt_strain))
    bamf.close()

    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "alt_strain"])
        for row in phased:
            w.writerow(row)
    print(f"# wrote {out}")
    print(f"# mode=meth-only  het_sites={n_sites:,}  "
          f"phased={len(phased):,} "
          f"({100*len(phased)/n_sites if n_sites else 0:.1f}%)")
    print(f"# strains: {', '.join(strains)}")
    return 0


def _joint_phase(bam, reference, vcf, labels, out, threads):
    """Joint methylation + SNP phasing (join_phasing.py core).

    Run whatshap to produce locally-phased blocks + per-read haplotags, then use
    methylation labels to ORIENT each whatshap block to a global strain, giving
    a genome-wide ALT-strain call per phased het site.
    """
    whatshap = _need("whatshap")
    bcftools = _need("bcftools")
    samtools = _need("samtools")

    workdir = tempfile.mkdtemp(prefix="epilect_phase_")
    try:
        # whatshap needs a bgzipped+indexed VCF
        vcf_in = vcf
        if not vcf_in.endswith(".gz"):
            vz = os.path.join(workdir, "in.vcf.gz")
            with open(vz, "wb") as o:
                subprocess.run([bcftools, "view", "-Oz", vcf_in], check=True, stdout=o)
            vcf_in = vz
        subprocess.run([bcftools, "index", "-f", vcf_in], check=True)

        phased_vcf = os.path.join(workdir, "phased.vcf.gz")
        sys.stderr.write("[epilect] whatshap phase ...\n")
        subprocess.run([whatshap, "phase", "-o", phased_vcf,
                        "--reference", reference, "--ignore-read-groups",
                        vcf_in, bam], check=True)
        subprocess.run([bcftools, "index", "-f", phased_vcf], check=True)

        htag_tsv = os.path.join(workdir, "haplotags.tsv")
        sys.stderr.write("[epilect] whatshap haplotag ...\n")
        subprocess.run([whatshap, "haplotag", "--reference", reference,
                        "--ignore-read-groups",
                        "--output-haplotag-list", htag_tsv,
                        "-o", os.path.join(workdir, "tagged.bam"),
                        phased_vcf, bam], check=True)

        # parse phased het sites: pos -> (PS, alt_hap)
        phased = {}
        res = subprocess.run([bcftools, "view", "-H", phased_vcf],
                             capture_output=True, text=True).stdout
        for line in res.splitlines():
            f = line.split("\t")
            if len(f) < 10:
                continue
            pos = int(f[1])
            d = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = d.get("GT", "")
            if "|" not in gt:
                continue
            alt_hap = "H1" if gt == "1|0" else ("H2" if gt == "0|1" else None)
            if alt_hap:
                phased[pos] = (d.get("PS", "NA"), alt_hap)

        # parse haplotag list: read -> (hap, PS)
        htag = {}
        with open(htag_tsv) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or f[1] in ("none", "None", ""):
                    continue
                hap = "H1" if f[1] in ("H1", "1") else ("H2" if f[1] in ("H2", "2") else None)
                if hap:
                    htag[f[0]] = (hap, f[2])
    finally:
        # keep VCF result before cleanup
        pass

    # orient each PS block by methylation-label majority
    strains = sorted(set(labels.values()))
    if len(strains) < 2:
        shutil.rmtree(workdir, ignore_errors=True)
        raise SystemExit("ERROR: joint phasing needs >=2 strain labels.")
    sA, sB = strains[0], strains[1]
    tally = collections.defaultdict(lambda: collections.Counter())
    for read, (hap, ps) in htag.items():
        s = labels.get(read)
        if s in (sA, sB):
            tally[(ps, hap)][s] += 1
    orient = {}
    for ps in {ps for (ps, _h) in tally}:
        h1 = tally.get((ps, "H1"), collections.Counter())
        h2 = tally.get((ps, "H2"), collections.Counter())
        s1 = h1[sA] - h1[sB]
        s2 = h2[sA] - h2[sB]
        orient[ps] = ({"H1": sA, "H2": sB} if s1 - s2 >= 0
                      else {"H1": sB, "H2": sA})

    contig = pysam_first_contig(bam)
    n_oriented = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "pos", "ps_block", "alt_strain"])
        for pos, (ps, alt_hap) in sorted(phased.items()):
            o = orient.get(ps)
            if not o:
                continue
            n_oriented += 1
            w.writerow([contig, pos, ps, o[alt_hap]])
    shutil.rmtree(workdir, ignore_errors=True)
    print(f"# wrote {out}")
    print(f"# mode=joint  het_phased(local)={len(phased):,}  "
          f"oriented_genomewide={n_oriented:,}  blocks={len(orient):,}")
    print(f"# strains: {sA}, {sB}")
    return 0


def pysam_first_contig(bam: str) -> str:
    import pysam
    bf = pysam.AlignmentFile(bam, "rb")
    c = bf.references[0]
    bf.close()
    return c


def cmd_phase(args) -> int:
    labels = _read_labels(args.labels)
    if not labels:
        raise SystemExit(f"ERROR: no usable read labels in {args.labels}.")
    if args.mode == "meth-only":
        return _meth_only_phase(args.bam, args.vcf, labels, args.out)
    else:
        if not args.reference:
            raise SystemExit("ERROR: joint mode requires --reference.")
        return _joint_phase(args.bam, args.reference, args.vcf, labels,
                            args.out, args.threads)


# --------------------------------------------------------------------------- #
# `resolvable`
# --------------------------------------------------------------------------- #
def _norm_species(s: str) -> str:
    toks = s.replace("[", "").replace("]", "").split()
    return " ".join(toks[:2]) if len(toks) >= 2 else s


def cmd_resolvable(args) -> int:
    if args.pair:
        files = args.pair
        if len(files) != 2:
            raise SystemExit("ERROR: --pair takes exactly two motif-set files.")
        sets = []
        for p in files:
            with open(p) as fh:
                sets.append({ln.strip().upper() for ln in fh if ln.strip()})
        tier = pair_tier(sets[0], sets[1])
        result = {
            "pair": [os.path.basename(files[0]), os.path.basename(files[1])],
            "tier": tier,
            "a_private": sorted(sets[0] - sets[1]),
            "b_private": sorted(sets[1] - sets[0]),
            "a_private_max_density": (max((density(m) for m in (sets[0] - sets[1])), default=0.0)),
            "b_private_max_density": (max((density(m) for m in (sets[1] - sets[0])), default=0.0)),
        }
        _emit_resolvable(result, args.out)
        print(f"# pair tier = {tier}")
        return 0

    # motif-table mode
    gmot = collections.defaultdict(set)
    gspec = {}
    with open(args.motif_table) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if "genome_id" not in (r.fieldnames or []) or "motif" not in (r.fieldnames or []):
            raise SystemExit(
                "ERROR: --motif-table needs tab-separated columns "
                "'genome_id', 'motif', and 'species'."
            )
        for row in r:
            gid = (row.get("genome_id") or "").strip()
            mot = (row.get("motif") or "").strip().upper()
            if not gid or not mot or mot == "MOTIF":
                continue
            gmot[gid].add(mot)
            gspec[gid] = _norm_species(row.get("species", gid))

    import itertools
    sp = collections.defaultdict(list)
    for gid, mots in gmot.items():
        sp[gspec[gid]].append(frozenset(mots))

    tiers = collections.Counter()
    n_species = 0
    n_pairs = 0
    for species, strains in sp.items():
        ss = [m for m in strains if len(m) >= args.min_motifs]
        if len(ss) < 2:
            continue
        n_species += 1
        sig = collections.Counter(ss)
        sigs = list(sig.items())
        for fs, c in sigs:
            if c >= 2:
                tiers["T0_none"] += c * (c - 1) // 2
        for (fa, ca), (fb, cb) in itertools.combinations(sigs, 2):
            tiers[pair_tier(fa, fb)] += ca * cb
        n = sum(sig.values())
        n_pairs += n * (n - 1) // 2

    T0 = tiers["T0_none"]
    TB = tiers["T1_bidirectional"]
    TO = tiers["T1_onesided"]
    T2 = tiers["T2_coverage"]
    N = T0 + TB + TO + T2
    result = {
        "min_motifs": args.min_motifs,
        "n_species": n_species,
        "n_pairs": N,
        "tiers": {"T1_bidirectional": TB, "T1_onesided": TO,
                  "T2_coverage": T2, "T0_none": T0},
    }
    _emit_resolvable(result, args.out)
    print(f"# min motifs/strain >= {args.min_motifs}: "
          f"species={n_species:,}  conspecific pairs={N:,}")
    if N:
        print(f"#   T1 bidirectional: {TB:,} ({100*TB/N:.1f}%)")
        print(f"#   T1 one-sided:     {TO:,} ({100*TO/N:.1f}%)")
        print(f"#   T2 coverage:      {T2:,} ({100*T2/N:.1f}%)")
        print(f"#   T0 unresolvable:  {T0:,} ({100*T0/N:.1f}%)")
        print(f"#   single-read resolvable (T1): {100*(TB+TO)/N:.1f}%   "
              f"any signal: {100*(TB+TO+T2)/N:.1f}%")
    return 0


def _emit_resolvable(result, out: Optional[str]):
    if out:
        with open(out, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"# wrote {out}")


# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="epilect",
        description="Per-read methylation strain deconvolution and phasing.",
    )
    p.add_argument("--version", action="version",
                   version=f"epilect {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # classify
    c = sub.add_parser("classify",
                       help="extract per-read calls (modkit) and run the "
                            "2-strain LLR classifier")
    c.add_argument("--bam", required=True, help="aligned modBAM (with MM/ML tags)")
    c.add_argument("--reference", required=True, help="reference FASTA")
    c.add_argument("--motifs", default=DEFAULT_MOTIFS,
                   help="comma list of MOTIF:OFFSET:CANON_BASE:METH_CODE "
                        f"(default: {DEFAULT_MOTIFS})")
    g = c.add_mutually_exclusive_group()
    g.add_argument("--emission", help="JSON {strainA:{motif:p}, strainB:{motif:p}}")
    g.add_argument("--estimate", action="store_true",
                   help="estimate a 2-cluster emission model from the data "
                        "instead of supplying --emission")
    c.add_argument("--out", required=True, help="output per-read TSV")
    c.add_argument("-t", "--threads", type=int, default=4)
    c.add_argument("--keep-calls", action="store_true",
                   help="keep the intermediate modkit calls TSV")
    c.set_defaults(func=cmd_classify)

    # phase
    ph = sub.add_parser("phase",
                        help="methylation-only or joint methylation+SNP phasing")
    ph.add_argument("--bam", required=True)
    ph.add_argument("--reference", help="reference FASTA (required for joint mode)")
    ph.add_argument("--vcf", required=True, help="het VCF (bgzipped or plain)")
    ph.add_argument("--labels", required=True,
                    help="per-read strain labels TSV (e.g. classify output)")
    ph.add_argument("--mode", choices=["joint", "meth-only"], default="joint")
    ph.add_argument("--out", required=True, help="output phased-site TSV")
    ph.add_argument("-t", "--threads", type=int, default=4)
    ph.set_defaults(func=cmd_phase)

    # resolvable
    r = sub.add_parser("resolvable",
                       help="tier conspecific pairs by methylation resolvability")
    src = r.add_mutually_exclusive_group(required=True)
    src.add_argument("--motif-table",
                     help="TSV with columns genome_id, species, motif")
    src.add_argument("--pair", nargs=2, metavar=("SET_A", "SET_B"),
                     help="two files, each one motif per line")
    r.add_argument("--min-motifs", type=int, default=3,
                   help="min motifs/strain to include (default 3)")
    r.add_argument("--out", help="optional output JSON")
    r.set_defaults(func=cmd_resolvable)

    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "classify" and not args.estimate and not args.emission:
        # default to the bundled SDSE model if neither given
        sys.stderr.write("[epilect] no --emission/--estimate; using built-in "
                         "SDSE GATC/CCWGG model\n")
        tmp = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(DEFAULT_EMISSION, tmp)
        tmp.close()
        args.emission = tmp.name
    try:
        return args.func(args)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"[epilect] external command failed: {e}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
