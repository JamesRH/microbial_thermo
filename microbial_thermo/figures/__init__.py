"""Figures.

These modules may import from the numerical core; the core must never import
from here. That keeps ``microbial_thermo`` usable headless and keeps plotting
dependencies out of the calculation path.
"""

from .explorer import plot_energy_explorer
from .halfreaction import plot_half_reactions
from .interactive_tower import interactive_tower, plot_interactive_tower
from .ladder import plot_affinity_ladder
from .syntrophy import SyntrophyWindow, plot_syntrophy_window, syntrophy_window
from .tower import plot_redox_tower

__all__ = [
    "SyntrophyWindow",
    "interactive_tower",
    "plot_affinity_ladder",
    "plot_energy_explorer",
    "plot_half_reactions",
    "plot_interactive_tower",
    "plot_redox_tower",
    "plot_syntrophy_window",
    "syntrophy_window",
]
