"""Epilect: turning per-read methylation noise into strain-resolving signal.

A small, tested toolkit for nanopore methylation strain deconvolution:
  - per-read 2-strain likelihood-ratio classification,
  - methylation-only and joint (methylation + SNP) phasing,
  - REBASE-based methylation-resolvability tiering of conspecific strain pairs.

The pure, importable, unit-testable numerics live in ``epilect.core``.
The command-line interface lives in ``epilect.cli`` (console script
``epilect``).
"""

__version__ = "0.1.0"

from .core import (
    EPS,
    clamp,
    loglik,
    classify,
    multi_loglik,
    multi_assign,
    density,
    pair_tier,
    DEG,
    DENSITY_THRESHOLD,
)

__all__ = [
    "__version__",
    "EPS",
    "clamp",
    "loglik",
    "classify",
    "multi_loglik",
    "multi_assign",
    "density",
    "pair_tier",
    "DEG",
    "DENSITY_THRESHOLD",
]
