"""Reference redox couples and their potentials.

Kept in the numerical core rather than in ``figures`` so the potentials can be
tabulated without importing a plotting library.

Potentials are computed from the backend at the requested pH and temperature
rather than read from a stored table, so a tower drawn at pH 6 and 60 C is
actually correct rather than a pH 7, 25 C table with a misleading axis label.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

import yaml

from .conditions import Conditions
from .exceptions import MicrobialThermoError
from .reaction import Couple, HalfReactionResult

_DATA_FILE = "tower_couples.yaml"


@dataclass(frozen=True)
class TowerEntry:
    """One couple on the tower, with its computed potential."""

    label: str
    couple: Couple
    group: str
    potential_v: float
    n_electrons: float


@lru_cache(maxsize=1)
def _couple_definitions() -> tuple:
    text = resources.files("microbial_thermo.data").joinpath(_DATA_FILE).read_text()
    return tuple(yaml.safe_load(text)["couples"])


def reference_couples(
    conditions: Conditions | None = None,
    backend=None,
    primed: bool = True,
    skip_unavailable: bool = True,
) -> list[TowerEntry]:
    """Compute potentials for the built-in reference couples.

    ``primed`` selects E-standard-prime (at the conditions' pH) over
    E-standard. Couples whose species are missing from the loaded database are
    skipped when ``skip_unavailable``, so a partial database still yields a
    usable tower.
    """
    from . import get_backend

    backend = backend or get_backend()
    conditions = conditions or Conditions()

    entries: list[TowerEntry] = []
    for definition in _couple_definitions():
        try:
            couple = Couple.make(
                definition["reduced"],
                definition["oxidized"],
                definition.get("key_element"),
            )
            half = couple.half_reaction()
            result = HalfReactionResult(
                couple=couple, half=half, conditions=conditions, backend=backend
            )
            potential = (result.E_standard_prime if primed else result.E_standard).to("V").magnitude
        except MicrobialThermoError:
            if skip_unavailable:
                continue
            raise
        entries.append(
            TowerEntry(
                label=definition["label"],
                couple=couple,
                group=definition.get("group", "other"),
                potential_v=float(potential),
                n_electrons=float(half.n_electrons),
            )
        )
    entries.sort(key=lambda e: e.potential_v)
    return entries


def tower_table(conditions: Conditions | None = None, backend=None, primed: bool = True):
    """The reference couples as a pandas DataFrame, ordered by potential."""
    import pandas as pd

    entries = reference_couples(conditions, backend, primed=primed)
    column = "E0' (V)" if primed else "E0 (V)"
    return pd.DataFrame(
        {
            "couple": [e.label for e in entries],
            column: [e.potential_v for e in entries],
            "n electrons": [e.n_electrons for e in entries],
            "group": [e.group for e in entries],
        }
    )
