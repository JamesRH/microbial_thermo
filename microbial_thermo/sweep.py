"""Sweeping a reaction over an environmental variable.

pyGCC costs roughly 28 ms per species per call and gains nothing from
vectorisation, so nothing interactive may call it per frame. Instead a sweep
is computed once on a grid here, and the figures interpolate or simply plot
the precomputed points. Baking the grid into the traces is also what lets an
exported HTML file keep working with no Python process behind it.

Lives in the core rather than in ``figures`` so sweeps can be tabulated
without importing a plotting library.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .conditions import Conditions
from .units import KJ_PER_MOL_STR, as_magnitude

#: Variables a reaction can be swept over.
PH = "pH"
TEMPERATURE = "temperature_c"
CONCENTRATION = "concentration"
PARTIAL_PRESSURE = "partial_pressure"


@dataclass
class SweepAxis:
    """One swept variable and the grid to evaluate it on."""

    kind: str
    label: str
    values: np.ndarray
    species: str | None = None
    log_scale: bool = False
    unit: str = ""

    @property
    def key(self) -> str:
        return f"{self.kind}:{self.species}" if self.species else self.kind

    def apply(self, conditions: Conditions, value: float) -> Conditions:
        """Return ``conditions`` with this variable set to ``value``."""
        if self.kind == PH:
            return conditions.replace(pH=float(value))
        if self.kind == TEMPERATURE:
            return conditions.replace(temperature_c=float(value))
        if self.kind == CONCENTRATION:
            updated = dict(conditions.concentrations)
            updated[self.species] = float(value)
            return conditions.replace(concentrations=updated)
        if self.kind == PARTIAL_PRESSURE:
            updated = dict(conditions.partial_pressures)
            updated[self.species] = float(value)
            return conditions.replace(partial_pressures=updated)
        raise ValueError(f"unknown sweep kind {self.kind!r}")


@dataclass
class SweepResult:
    """Free energies along one swept axis."""

    axis: SweepAxis
    delta_g: np.ndarray
    delta_g_per_electron: np.ndarray
    n_electrons: float
    failures: list[tuple[float, str]] = field(default_factory=list)

    def to_frame(self):
        import pandas as pd

        return pd.DataFrame(
            {
                self.axis.label: self.axis.values,
                "dG (kJ/mol)": self.delta_g,
                "dG per e- (kJ/mol)": self.delta_g_per_electron,
            }
        )


def ph_axis(low: float = 3.0, high: float = 11.0, n: int = 33) -> SweepAxis:
    return SweepAxis(kind=PH, label="pH", values=np.linspace(low, high, n))


def temperature_axis(low: float = 0.01, high: float = 100.0, n: int = 31) -> SweepAxis:
    """Temperature grid. The lower bound avoids exactly 0 C, where the
    IAPWS-95 water equation of state fails."""
    return SweepAxis(
        kind=TEMPERATURE,
        label="temperature (°C)",
        values=np.linspace(max(low, 0.01), min(high, 100.0), n),
        unit="°C",
    )


def partial_pressure_axis(
    species: str, low: float = 1e-10, high: float = 1.0, n: int = 25
) -> SweepAxis:
    return SweepAxis(
        kind=PARTIAL_PRESSURE,
        label=f"p{species} (bar)",
        values=np.logspace(np.log10(low), np.log10(high), n),
        species=species,
        log_scale=True,
        unit="bar",
    )


def concentration_axis(
    species: str, low: float = 1e-9, high: float = 1e-1, n: int = 25
) -> SweepAxis:
    return SweepAxis(
        kind=CONCENTRATION,
        label=f"[{species}] (molal)",
        values=np.logspace(np.log10(low), np.log10(high), n),
        species=species,
        log_scale=True,
        unit="molal",
    )


def sweep(reaction, axis: SweepAxis, skip_failures: bool = True) -> SweepResult:
    """Evaluate ``reaction``'s free energy across ``axis``.

    Rebuilds the reaction at each grid point so that temperature sweeps pick
    up the temperature dependence of every formation energy, not just the
    RT ln Q factor.
    """
    from .reaction import Reaction

    energies, per_electron, failures = [], [], []
    for value in axis.values:
        conditions = axis.apply(reaction.conditions, value)
        try:
            point = Reaction.from_couples(
                donor=reaction.donor,
                acceptor=reaction.acceptor,
                conditions=conditions,
                backend=reaction.backend,
                n_electrons=reaction.n_electrons,
            )
            energy = as_magnitude(point.delta_G, KJ_PER_MOL_STR)
        except Exception as exc:
            if not skip_failures:
                raise
            failures.append((float(value), f"{type(exc).__name__}: {exc}"))
            energy = float("nan")
        energies.append(energy)
        per_electron.append(energy / float(reaction.n_electrons))

    return SweepResult(
        axis=axis,
        delta_g=np.asarray(energies, dtype=float),
        delta_g_per_electron=np.asarray(per_electron, dtype=float),
        n_electrons=float(reaction.n_electrons),
        failures=failures,
    )


def default_axes(reaction) -> list[SweepAxis]:
    """A sensible set of axes for ``reaction``: pH, temperature, and every
    gas or aqueous species that actually participates."""
    axes = [ph_axis(), temperature_axis()]
    for species in sorted(reaction.coefficients, key=lambda s: s.backend):
        if species.backend in ("H2O", "H+"):
            continue  # pH and water activity are handled separately
        if species.is_gas:
            axes.append(partial_pressure_axis(species.backend))
        elif species.is_aqueous:
            axes.append(concentration_axis(species.backend))
    return axes
