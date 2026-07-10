"""Deterministic unit tests for epilect.core."""
import math

import pytest

from epilect.core import (
    EPS,
    clamp,
    loglik,
    classify,
    multi_loglik,
    multi_assign,
    density,
    pair_tier,
)


# --------------------------------------------------------------------------- #
# clamp
# --------------------------------------------------------------------------- #
def test_clamp_bounds():
    # values pushed to the EPS / 1-EPS edges
    assert clamp(0.0) == EPS
    assert clamp(-5.0) == EPS
    assert clamp(1.0) == 1 - EPS
    assert clamp(2.0) == 1 - EPS
    # interior values pass through untouched
    assert clamp(0.5) == 0.5
    assert clamp(EPS) == EPS
    assert clamp(1 - EPS) == 1 - EPS
    # clamp keeps log() finite at the extremes
    assert math.isfinite(math.log(clamp(0.0)))
    assert math.isfinite(math.log(1 - clamp(1.0)))


# --------------------------------------------------------------------------- #
# loglik
# --------------------------------------------------------------------------- #
def test_loglik_correct_strain_scores_higher():
    # SDSE-style model: strain A methylates GATC, strain B methylates CCWGG.
    emit_a = {"GATC": 0.903, "CCWGG": 0.005}
    emit_b = {"GATC": 0.041, "CCWGG": 0.873}
    # a read that is heavily GATC-methylated and CCWGG-unmethylated -> strain A
    rec = {"GATC": (10, 9), "CCWGG": (8, 0)}
    assert loglik(rec, emit_a) > loglik(rec, emit_b)
    # the mirror-image read -> strain B
    rec_b = {"GATC": (10, 0), "CCWGG": (8, 7)}
    assert loglik(rec_b, emit_b) > loglik(rec_b, emit_a)


def test_loglik_matches_hand_value():
    emit = {"GATC": 0.9}
    rec = {"GATC": (2, 1)}  # 1 meth, 1 unmeth
    expected = 1 * math.log(0.9) + 1 * math.log(0.1)
    assert loglik(rec, emit) == pytest.approx(expected)


def test_loglik_ignores_unknown_channel():
    emit = {"GATC": 0.9}
    assert loglik({"GATC": (1, 1), "OTHER": (5, 5)}, emit) == pytest.approx(math.log(0.9))


# --------------------------------------------------------------------------- #
# classify
# --------------------------------------------------------------------------- #
def test_classify_label_and_sign():
    emit_a = {"GATC": 0.903, "CCWGG": 0.005}
    emit_b = {"GATC": 0.041, "CCWGG": 0.873}
    rec_a = {"GATC": (10, 9), "CCWGG": (8, 0)}
    label, llr = classify(rec_a, emit_a, emit_b, "A", "B")
    assert label == "A"
    assert llr > 0

    rec_b = {"GATC": (10, 0), "CCWGG": (8, 7)}
    label, llr = classify(rec_b, emit_a, emit_b, "A", "B")
    assert label == "B"
    assert llr < 0


def test_classify_tie_goes_to_b():
    # identical emission -> llr exactly 0 -> not > 0 -> label_b
    emit = {"GATC": 0.5}
    label, llr = classify({"GATC": (4, 2)}, emit, emit, "A", "B")
    assert llr == 0.0
    assert label == "B"


# --------------------------------------------------------------------------- #
# multi_assign / multi_loglik
# --------------------------------------------------------------------------- #
def test_multi_assign_argmax_and_margin():
    strains = ["SA62", "SA63", "SA67"]
    # motif m1 is methylated only in SA62; m2 only in SA67.
    emit = {
        "m1": {"SA62": 0.8, "SA63": 0.02, "SA67": 0.02},
        "m2": {"SA62": 0.02, "SA63": 0.02, "SA67": 0.8},
    }
    rec = {"m1": (5, 5), "m2": (5, 0)}  # strongly SA62
    best, margin, ll = multi_assign(rec, emit, strains)
    assert best == "SA62"
    assert margin > 0
    # SA62 should be the unique top, with SA63/SA67 below
    assert ll["SA62"] == max(ll.values())


def test_multi_assign_single_strain_margin_inf():
    best, margin, _ = multi_assign({"m1": (3, 3)}, {"m1": {"X": 0.9}}, ["X"])
    assert best == "X"
    assert margin == math.inf


# --------------------------------------------------------------------------- #
# density
# --------------------------------------------------------------------------- #
def test_density_gatc_vs_nnnn():
    # GATC: 4 fully specified bases -> (1/4)^4 per strand, *2 strands *10kb.
    expected_gatc = 2 * 10000 * (0.25 ** 4)
    assert density("GATC") == pytest.approx(expected_gatc)
    # NNNN: every position matches -> p=1 -> dense (2*10000).
    assert density("NNNN") == pytest.approx(2 * 10000 * 1.0)
    # NNNN must be far denser than GATC
    assert density("NNNN") > density("GATC")
    # empty motif -> 0
    assert density("") == 0.0
    # unknown token collapses to ~0
    assert density("GAT1") == pytest.approx(1e-9)


def test_density_degenerate_codes():
    # GATY: Y is 2-fold degenerate -> (1/4)^3 * (2/4)
    expected = 2 * 10000 * (0.25 ** 3) * 0.5
    assert density("GATY") == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# pair_tier
# --------------------------------------------------------------------------- #
def test_pair_tier_t0_identical():
    s = {"GATC", "CCWGG"}
    assert pair_tier(s, set(s)) == "T0_none"


def test_pair_tier_t1_bidirectional():
    # each strain has a private DENSE motif (GATC is dense at threshold 3/10kb?).
    # GATC density ~ 78/10kb >= 3 -> dense. Use two distinct dense motifs.
    a = {"GATC"}              # private dense
    b = {"CCWGG"}             # CCWGG: 5 positions, one W (2-fold) -> dense too
    # confirm both are dense
    assert density("GATC") >= 3.0
    assert density("CCWGG") >= 3.0
    assert pair_tier(a, b) == "T1_bidirectional"


def test_pair_tier_t1_onesided():
    # one strain has a private dense motif, the other only a sparse private one.
    sparse = "GAATTCGGCC"   # 10 fully-specified bases -> (1/4)^10 -> ~tiny
    assert density(sparse) < 3.0
    a = {"GATC"}            # dense private
    b = {sparse}           # sparse private
    assert pair_tier(a, b) == "T1_onesided"


def test_pair_tier_t2_coverage():
    # sets differ but only via sparse private motifs on both sides.
    sparse1 = "GAATTCGGCC"
    sparse2 = "TTAATTGGCC"
    assert density(sparse1) < 3.0 and density(sparse2) < 3.0
    assert pair_tier({sparse1}, {sparse2}) == "T2_coverage"
