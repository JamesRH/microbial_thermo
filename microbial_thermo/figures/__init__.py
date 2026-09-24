"""Figures.

These modules may import from the numerical core; the core must never import
from here. That keeps ``microbial_thermo`` usable headless and keeps plotting
dependencies out of the calculation path.
"""

from .explorer import plot_energy_explorer
from .halfreaction import plot_half_reactions
from .interactive_tower import interactive_tower, plot_interactive_tower
from .ladder import plot_affinity_ladder
from .light import light_budget, light_budget_table, plot_light_budget
from .syntrophy import SyntrophyWindow, plot_syntrophy_window, syntrophy_window
from .tower import plot_redox_tower

__all__ = [
    "SyntrophyWindow",
    "interactive_tower",
    "light_budget",
    "light_budget_table",
    "plot_affinity_ladder",
    "plot_energy_explorer",
    "plot_half_reactions",
    "plot_interactive_tower",
    "plot_light_budget",
    "plot_redox_tower",
    "plot_syntrophy_window",
    "syntrophy_window",
]
