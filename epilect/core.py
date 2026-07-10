"""Pure, importable, unit-testable numerics for Epilect.

These functions are refactored verbatim (behaviour-identical) from the original
research scripts so they can be imported and tested in isolation:

  - ``clamp``                : keep an emission probability off {0, 1}.
                               (scripts/per_read_classify.py:16, ef_classify.py:9, ...)
  - ``loglik`` / ``classify``: 2-strain Bernoulli log-likelihood and the
                               likelihood-ratio classifier.
                               (scripts/per_read_classify.py:37-47)
  - ``multi_loglik`` /
    ``multi_assign``         : N-strain log-likelihood, argmax assignment and
                               the confidence margin (gap between best and
                               second-best). (scripts/sa3_deconvolve.py:20-40)
  - ``density`` / ``pair_tier``: REBASE motif density and the resolvability
                               tier of a conspecific strain pair.
                               (scripts/resolvable_fraction.py:16-57)

Nothing here does any I/O; the CLI and the scripts wire these into real data.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, Mapping, Sequence, Tuple

# --------------------------------------------------------------------------- #
# Emission-probability clamping
# --------------------------------------------------------------------------- #
EPS = 1e-3


def clamp(p: float) -> float:
    """Clamp an emission probability into ``[EPS, 1 - EPS]``.

    Keeps probabilities strictly inside the open interval (0, 1) so that
    ``log(p)`` and ``log(1 - p)`` are always finite. Behaviour-identical to
    the ``clamp``/``cl`` helpers in the original scripts.
    """
    return min(1 - EPS, max(EPS, p))


# --------------------------------------------------------------------------- #
# 2-strain likelihood-ratio classifier
# --------------------------------------------------------------------------- #
def loglik(
    counts: Mapping[str, Tuple[int, int]],
    emission: Mapping[str, float],
) -> float:
    """Bernoulli log-likelihood of a read's calls under one strain's model.

    ``counts``   maps each motif/channel name -> ``(n_sites, n_methylated)``
                 observed on the read.
    ``emission`` maps each motif/channel name -> P(methylated | that strain).

    Channels present in ``counts`` but missing from ``emission`` are ignored.
    Equivalent to ``per_read_classify.loglik`` generalised over named channels.
    """
    ll = 0.0
    for motif, (n, m) in counts.items():
        if motif not in emission:
            continue
        p = clamp(emission[motif])
        ll += m * math.log(p) + (n - m) * math.log(1 - p)
    return ll


def classify(
    counts: Mapping[str, Tuple[int, int]],
    emission_a: Mapping[str, float],
    emission_b: Mapping[str, float],
    label_a: str,
    label_b: str,
) -> Tuple[str, float]:
    """Classify a read between two strains by log-likelihood ratio.

    Returns ``(label, llr)`` where ``llr = loglik(A) - loglik(B)``. The read is
    assigned to ``label_a`` when ``llr > 0`` and to ``label_b`` otherwise. This
    matches the tie-handling (``>0``) of the original ``classify`` in
    per_read_classify.py and ``ef_classify.llr``.
    """
    llr = loglik(counts, emission_a) - loglik(counts, emission_b)
    return (label_a if llr > 0 else label_b), llr


# --------------------------------------------------------------------------- #
# N-strain (multi-strain) assignment + confidence margin
# --------------------------------------------------------------------------- #
def multi_loglik(
    counts: Mapping[str, Tuple[int, int]],
    emission: Mapping[str, Mapping[str, float]],
    strains: Sequence[str],
) -> Dict[str, float]:
    """Per-strain Bernoulli log-likelihood over named motif channels.

    ``emission[motif][strain]`` is P(methylated | motif, strain). Motifs in
    ``counts`` absent from ``emission`` are skipped. Mirrors the per-strain
    accumulation in sa3_deconvolve.classify.
    """
    ll = {s: 0.0 for s in strains}
    for motif, (n, m) in counts.items():
        e = emission.get(motif)
        if e is None:
            continue
        for s in strains:
            p = clamp(e[s])
            ll[s] += m * math.log(p) + (n - m) * math.log(1 - p)
    return ll


def multi_assign(
    counts: Mapping[str, Tuple[int, int]],
    emission: Mapping[str, Mapping[str, float]],
    strains: Sequence[str],
) -> Tuple[str, float, Dict[str, float]]:
    """Assign a read to the argmax strain and report the confidence margin.

    Returns ``(best_strain, margin, per_strain_loglik)`` where ``margin`` is the
    log-likelihood gap between the best and second-best strain (the tie/abstain
    knob used in sa3_deconvolve). With a single strain the margin is ``inf``.
    Ties between strains are broken by the order given in ``strains``.
    """
    ll = multi_loglik(counts, emission, strains)
    ordered = sorted(strains, key=lambda s: -ll[s])
    best = ordered[0]
    margin = math.inf if len(ordered) < 2 else ll[best] - ll[ordered[1]]
    return best, margin, ll


# --------------------------------------------------------------------------- #
# REBASE motif density + resolvability tiers
# --------------------------------------------------------------------------- #
# IUPAC degeneracy: how many of the 4 bases each code matches.
DEG: Dict[str, int] = {
    "A": 1, "C": 1, "G": 1, "T": 1,
    "R": 2, "Y": 2, "S": 2, "W": 2, "K": 2, "M": 2,
    "B": 3, "D": 3, "H": 3, "V": 3, "N": 4,
}
READ_BP = 10000  # bp window for "sites per read-length" density
DENSITY_THRESHOLD = 3.0  # sites per 10 kb that count as "single-read dense"


def density(motif: str, read_bp: int = READ_BP) -> float:
    """Expected occurrences of a motif per ``read_bp`` bp, counting both strands.

    Computed from IUPAC degeneracy: P(match at a position) = prod(deg/4).
    Unknown tokens collapse the motif to ~0 density (treated as ultra-rare),
    and an empty motif returns 0.0. Behaviour-identical to
    resolvable_fraction.density.
    """
    m = motif.strip().upper()
    if not m:
        return 0.0
    p = 1.0
    for ch in m:
        d = DEG.get(ch)
        if d is None:
            return 1e-9  # unknown token -> ultra-rare / sparse
        p *= d / 4.0
    return 2 * read_bp * p  # both strands


def pair_tier(
    motifs_a: Iterable[str],
    motifs_b: Iterable[str],
    threshold: float = DENSITY_THRESHOLD,
    read_bp: int = READ_BP,
) -> str:
    """Methylation-resolvability tier for one conspecific strain pair.

    Given the methylated-motif sets of two strains, returns one of:
      - ``"T0_none"``          : identical methylated-motif sets (unresolvable).
      - ``"T1_bidirectional"`` : each strain has a *private* motif dense enough
                                 (>= ``threshold`` sites/``read_bp``) to recover
                                 it as the minor strain (SDSE-like).
      - ``"T1_onesided"``      : only one strain has a private dense motif
                                 (EF-like).
      - ``"T2_coverage"``      : motif sets differ, but only via sparse private
                                 motifs -> resolvable only with coverage, not
                                 single reads.

    Behaviour-identical to resolvable_fraction.pair_tier.
    """
    A = set(motifs_a)
    B = set(motifs_b)
    a_priv = A - B
    b_priv = B - A
    if not a_priv and not b_priv:
        return "T0_none"
    a_dense = bool(a_priv) and max(density(m, read_bp) for m in a_priv) >= threshold
    b_dense = bool(b_priv) and max(density(m, read_bp) for m in b_priv) >= threshold
    if a_dense and b_dense:
        return "T1_bidirectional"
    if a_dense or b_dense:
        return "T1_onesided"
    return "T2_coverage"
