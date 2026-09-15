"""The backend contract.

A backend supplies standard-state formation free energies and solvent
properties. It knows nothing about reactions, half reactions, or figures.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

from ..units import Quantity


@dataclass(frozen=True)
class SpeciesRecord:
    """What a backend knows about one species.

    ``source`` records where the number came from so that every downstream
    result can carry provenance.
    """

    name: str
    backend_name: str
    source: str
    verified: bool = True


class ThermoBackend(ABC):
    """Standard-state thermodynamic data at a temperature and pressure."""

    #: Human-readable backend identity, recorded in provenance.
    name: str = "abstract"

    @property
    @abstractmethod
    def version(self) -> str:
        """Version string of the underlying data package."""

    @abstractmethod
    def resolve(self, name: str) -> SpeciesRecord:
        """Map a user-facing species name onto a backend species.

        Raises :class:`~microbial_thermo.exceptions.SpeciesNotFoundError` with
        near-miss suggestions when the name cannot be resolved.
        """

    @abstractmethod
    def delta_Gf(self, name: str, temperature_c: float, pressure_bar=None) -> Quantity:
        """Standard-state Gibbs free energy of formation, in kJ/mol."""

    @abstractmethod
    def solvent_properties(self, temperature_c: float, pressure_bar=None) -> dict:
        """Debye-Huckel A and B and the b-dot parameter at T, P."""

    @abstractmethod
    def available_species(self) -> list[str]:
        """Every species name the backend can supply."""

    def temperature_range_c(self) -> tuple[float, float]:
        """Supported temperature window, in Celsius."""
        return (0.01, 100.0)
