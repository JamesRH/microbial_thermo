"""microbial_thermo -- thermodynamics for microbial physiology and biogeochemistry.

Typical use::

    import microbial_thermo as mt

    rxn = mt.Reaction.from_couples(
        donor=("H2(g)", "H+"),          # (reduced, oxidized)
        acceptor=("HS-", "SO4-2"),
        conditions=mt.Conditions(temperature_c=25.0, pH=7.0),
    )
    print(rxn.summary())

The thermodynamic backend is selected once, at first use, and shared thereafter.
pyGCC is the only backend currently implemented.
"""

from __future__ import annotations

from .balance import (
    HalfReaction,
    balance_equation,
    balance_half_reaction,
    combine_half_reactions,
)
from .conditions import Conditions
from .exceptions import (
    AmbiguousReactionError,
    BalancingError,
    MicrobialThermoError,
    MissingDataError,
    OutOfRangeError,
    SpeciesNotFoundError,
    ThermodynamicConsistencyError,
)
from .library import default_library, energy_table, metabolism
from .oxidation import atom_oxidation_states, mean_oxidation_state, nosc
from .reaction import Couple, HalfReactionResult, Reaction, half_reaction
from .speciation import (
    default_families,
    dominant,
    fractions,
    pKa,
    pKa_ladder,
    speciation_table,
)
from .species import Species, default_registry, resolve
from .units import Quantity, ureg

__version__ = "0.1.0"

_BACKEND = None
_BACKEND_NAME = "pygcc"


def configure(backend: str = "pygcc", **kwargs) -> object:
    """Select and configure the thermodynamic backend.

    Called implicitly on first use with the defaults, so most scripts never
    need it. Pass ``database=`` to point pyGCC at a different data file.
    """
    global _BACKEND, _BACKEND_NAME
    name = backend.lower()
    if name != "pygcc":
        raise ValueError(f"unknown backend {backend!r}; pyGCC is the only backend implemented")
    from .backends import PygccBackend

    _BACKEND = PygccBackend(**kwargs)
    _BACKEND_NAME = name
    return _BACKEND


def get_backend() -> object:
    """The active backend, configuring the default on first use."""
    global _BACKEND
    if _BACKEND is None:
        configure(_BACKEND_NAME)
    return _BACKEND


def reset_backend() -> None:
    """Drop the configured backend, clearing its caches."""
    global _BACKEND
    _BACKEND = None


def versions() -> dict[str, str]:
    """Versions of this library and its computational dependencies."""
    import importlib.metadata as md

    out = {"microbial_thermo": __version__}
    for package in (
        "pygcc",
        "numpy",
        "pandas",
        "scipy",
        "chempy",
        "sympy",
        "pint",
        "matplotlib",
        "plotly",
    ):
        try:
            out[package] = md.version(package)
        except md.PackageNotFoundError:
            out[package] = "not installed"
    return out


__all__ = [
    "Conditions",
    "half_reaction",
    "default_library",
    "energy_table",
    "metabolism",
    "default_families",
    "dominant",
    "fractions",
    "pKa",
    "pKa_ladder",
    "speciation_table",
    "Couple",
    "HalfReaction",
    "HalfReactionResult",
    "Reaction",
    "Species",
    "balance_equation",
    "balance_half_reaction",
    "combine_half_reactions",
    "configure",
    "default_registry",
    "get_backend",
    "reset_backend",
    "resolve",
    "versions",
    "mean_oxidation_state",
    "nosc",
    "atom_oxidation_states",
    "Quantity",
    "ureg",
    "__version__",
    "MicrobialThermoError",
    "SpeciesNotFoundError",
    "BalancingError",
    "AmbiguousReactionError",
    "ThermodynamicConsistencyError",
    "MissingDataError",
    "OutOfRangeError",
]
