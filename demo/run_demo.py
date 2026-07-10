#!/usr/bin/env python3
"""Epilect demo: classify the bundled SDSE mock and regenerate a figure panel.

Runs ``epilect classify`` on demo/sdse_demo.bam (a ~1.4k-read downsample of
the C40_k8 SDSE mock: a 1:8 minor:major mix of UT10237:UT9728), then plots the
per-read LLR-score distribution coloured by inferred strain, and reports the
estimated minor-strain abundance vs the truth used to build the mock.

The mock was built at total coverage 40x, skew 1:8, so the expected minor
(UT10237) fraction is 1/(1+8) = 11.11%. The panel shows how many reads land on
the minor strain after classification. This is the headline "noise -> signal"
panel in miniature: methylation alone separates the two strains per read.
"""
from __future__ import annotations

import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BAM = os.path.join(HERE, "sdse_demo.bam")
REF = os.path.join(HERE, "sdse_ref.fna")
OUT_TSV = os.path.join(HERE, "demo_classify.tsv")
PANEL = os.path.join(HERE, "demo_panel.png")

# The demo BAM is a downsample of C40_k8 -> nominal skew 1:8 -> minor frac.
EXPECTED_MINOR_FRAC = 1.0 / (1.0 + 8.0)
MINOR = "UT10237"   # absence-of-GATC strain (minor in this mock)
MAJOR = "UT9728"


def run_classify() -> None:
    cmd = [sys.executable, "-m", "epilect.cli", "classify",
           "--bam", BAM, "--reference", REF,
           "--motifs", "GATC:1:A:a,CCWGG:0:C:m",
           "--out", OUT_TSV, "-t", "4"]
    print("[demo] " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def load_rows():
    rows = []
    with open(OUT_TSV) as fh:
        next(fh)
        for line in fh:
            rid, strain, llr, n = line.rstrip("\n").split("\t")
            rows.append((rid, strain, float(llr), int(n)))
    return rows


def main() -> int:
    run_classify()
    rows = load_rows()
    classifiable = [r for r in rows if r[1] != "unclassifiable"]
    n_total = len(rows)
    n_cls = len(classifiable)
    n_minor = sum(1 for r in classifiable if r[1] == MINOR)
    obs_minor_frac = n_minor / n_cls if n_cls else 0.0

    # ---- panel: LLR score distribution coloured by inferred strain ----
    minor_scores = [r[2] for r in classifiable if r[1] == MINOR]
    major_scores = [r[2] for r in classifiable if r[1] == MAJOR]

    fig, ax = plt.subplots(figsize=(6.2, 4.0), dpi=150)
    bins = 40
    ax.hist(major_scores, bins=bins, color="#3b6ea5", alpha=0.8,
            label=f"{MAJOR} (major)  n={len(major_scores)}")
    ax.hist(minor_scores, bins=bins, color="#c44e52", alpha=0.8,
            label=f"{MINOR} (minor)  n={len(minor_scores)}")
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set_yscale("symlog")
    ax.set_xlabel("per-read log-likelihood ratio  (>0 -> %s)" % MAJOR)
    ax.set_ylabel("reads (symlog)")
    ax.set_title("Epilect demo: per-read strain separation\n(SDSE 1:8 mock)",
                 fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    txt = (f"expected minor fraction = {100*EXPECTED_MINOR_FRAC:.2f}%\n"
           f"observed minor fraction = {100*obs_minor_frac:.2f}%")
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=8, family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    fig.tight_layout()
    fig.savefig(PANEL)
    plt.close(fig)

    print()
    print(f"[demo] reads total={n_total}  classifiable={n_cls} "
          f"({100*n_cls/n_total:.1f}%)")
    print(f"[demo] minor strain = {MINOR}")
    print(f"[demo] EXPECTED minor fraction = {100*EXPECTED_MINOR_FRAC:.2f}%  "
          f"OBSERVED minor fraction = {100*obs_minor_frac:.2f}%")
    print(f"[demo] panel written: {PANEL}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
