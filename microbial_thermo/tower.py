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


#: Default grid for the interactive tower. pH is sampled finely because it is
#: nearly free -- formation energies are cached per temperature, and pH enters
#: only through the proton activity term. Temperature is sampled coarsely
#: because each new value costs a full round of backend calls.
DEFAULT_GRID_PH = tuple(x / 2 for x in range(6, 23))  # 3.0 to 11.0 by 0.5
DEFAULT_GRID_TEMPERATURES = (4.0, 15.0, 25.0, 37.0, 55.0, 80.0)


@dataclass(frozen=True)
class TowerGrid:
    """Reference-couple potentials precomputed over pH and temperature.

    The interactive tower cannot call the backend per slider movement --
    a single temperature change costs a couple of seconds -- so the whole
    surface is computed once and then indexed. The activity-ratio dimension is
    not stored: it is a closed-form Nernst shift, applied on lookup.

    ``potentials`` is indexed ``[pH, temperature, couple]`` and is NaN where a
    couple could not be computed at that grid point.
    """

    ph_values: tuple[float, ...]
    temperature_values: tuple[float, ...]
    labels: tuple[str, ...]
    groups: tuple[str, ...]
    n_electrons: tuple[float, ...]
    potentials: object  # numpy array [pH, temperature, couple]
    primed: bool = True

    def nearest_index(self, axis: str, value: float) -> int:
        """Index of the closest grid point. Sliders are restricted to grid
        values, so this is only for callers passing arbitrary numbers."""
        import numpy as np

        values = self.ph_values if axis == "pH" else self.temperature_values
        return int(np.argmin(np.abs(np.asarray(values) - value)))

    def at(self, pH: float, temperature_c: float, ratio: float = 1.0):
        """Potentials in volts at the nearest grid point, Nernst-shifted.

        ``ratio`` is the oxidised-to-reduced activity ratio, applied to every
        couple alike. The shift is ``(RT/nF) ln(ratio)``, so it is larger for
        a one-electron couple than for an eight-electron one -- which is the
        reason the slider is worth having rather than a curiosity.
        """
        import numpy as np

        from .units import FARADAY, R, celsius_to_kelvin, ureg

        base = self.potentials[self.nearest_index("pH", pH), self.nearest_index("T", temperature_c)]
        if ratio == 1.0:
            return base
        kelvin = celsius_to_kelvin(self.temperature_values[self.nearest_index("T", temperature_c)])
        # RT/F, in volts. The kelvin must carry its unit or it will not cancel
        # against the gas constant.
        slope = (R * (kelvin * ureg.kelvin) / FARADAY.to("C/mol")).to("V").magnitude
        return base + slope * np.log(ratio) / np.asarray(self.n_electrons)

    def entries_at(self, pH: float, temperature_c: float, ratio: float = 1.0):
        """The same thing as :class:`TowerEntry` objects, ordered by potential
        and with unavailable couples dropped."""
        import numpy as np

        potentials = self.at(pH, temperature_c, ratio)
        rows = [
            (label, group, n, float(v))
            for label, group, n, v in zip(
                self.labels, self.groups, self.n_electrons, potentials, strict=True
            )
            if np.isfinite(v)
        ]
        rows.sort(key=lambda r: r[3])
        return rows

    def to_frame(self, pH: float, temperature_c: float, ratio: float = 1.0):
        import pandas as pd

        rows = self.entries_at(pH, temperature_c, ratio)
        column = "E0' (V)" if self.primed else "E0 (V)"
        return pd.DataFrame(
            {
                "couple": [r[0] for r in rows],
                column: [r[3] for r in rows],
                "n electrons": [r[2] for r in rows],
                "group": [r[1] for r in rows],
            }
        )


def tower_grid(
    ph_values=DEFAULT_GRID_PH,
    temperature_values=DEFAULT_GRID_TEMPERATURES,
    backend=None,
    primed: bool = True,
    progress=None,
) -> TowerGrid:
    """Precompute reference potentials across pH and temperature.

    Costs roughly a couple of seconds per temperature and almost nothing per
    pH. Pass ``progress`` a callable taking ``(done, total)`` to report back
    during the wait.

    Couples are collected across the whole grid rather than at one corner, so
    a couple available at 25 C but not at 80 C occupies a row of NaN instead
    of silently shortening the tower.
    """
    import numpy as np

    ph_values = tuple(float(x) for x in ph_values)
    temperature_values = tuple(float(x) for x in temperature_values)

    collected: dict[str, dict] = {}
    cells: dict[tuple[int, int], dict[str, float]] = {}
    total = len(ph_values) * len(temperature_values)
    done = 0
    for j, temperature in enumerate(temperature_values):
        for i, pH in enumerate(ph_values):
            conditions = Conditions(temperature_c=temperature, pH=pH)
            entries = reference_couples(conditions, backend, primed=primed)
            cells[(i, j)] = {e.label: e.potential_v for e in entries}
            for entry in entries:
                collected.setdefault(
                    entry.label,
                    {"group": entry.group, "n_electrons": entry.n_electrons},
                )
            done += 1
            if progress is not None:
                progress(done, total)

    labels = tuple(collected)
    potentials = np.full((len(ph_values), len(temperature_values), len(labels)), np.nan)
    for (i, j), found in cells.items():
        for k, label in enumerate(labels):
            if label in found:
                potentials[i, j, k] = found[label]

    return TowerGrid(
        ph_values=ph_values,
        temperature_values=temperature_values,
        labels=labels,
        groups=tuple(collected[x]["group"] for x in labels),
        n_electrons=tuple(collected[x]["n_electrons"] for x in labels),
        potentials=potentials,
        primed=primed,
    )
