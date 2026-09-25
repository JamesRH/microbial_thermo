"""Figures.

These modules may import from the numerical core; the core must never import
from here. That keeps ``microbial_thermo`` usable headless and keeps plotting
dependencies out of the calculation path.
"""

from .explorer import plot_energy_explorer
from .frost import frost_diagram, plot_frost
from .halfreaction import plot_half_reactions
from .interactive_tower import interactive_tower, plot_interactive_tower
from .ladder import plot_affinity_ladder
from .latimer import latimer_diagram, plot_latimer
from .light import light_budget, light_budget_table, plot_light_budget
from .pourbaix import plot_pourbaix, pourbaix_field, water_stability
from .syntrophy import SyntrophyWindow, plot_syntrophy_window, syntrophy_window
from .tower import plot_redox_tower

__all__ = [
    "SyntrophyWindow",
    "frost_diagram",
    "interactive_tower",
    "latimer_diagram",
    "light_budget",
    "light_budget_table",
    "plot_affinity_ladder",
    "plot_energy_explorer",
    "plot_frost",
    "plot_half_reactions",
    "plot_interactive_tower",
    "plot_latimer",
    "plot_pourbaix",
    "pourbaix_field",
    "plot_light_budget",
    "plot_redox_tower",
    "plot_syntrophy_window",
    "syntrophy_window",
    "water_stability",
]
