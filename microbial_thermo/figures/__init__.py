"""Figures.

These modules may import from the numerical core; the core must never import
from here. That keeps ``microbial_thermo`` usable headless and keeps plotting
dependencies out of the calculation path.
"""

from .halfreaction import plot_half_reactions

__all__ = ["plot_half_reactions"]
