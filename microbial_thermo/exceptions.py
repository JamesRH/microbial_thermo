"""Exception hierarchy.

Errors are specific so that a teaching context can explain *what* went wrong
rather than surfacing a bare KeyError from deep inside a database lookup.
"""

from __future__ import annotations


class MicrobialThermoError(Exception):
    """Base class for every error raised by this library."""


class SpeciesNotFoundError(MicrobialThermoError):
    """A species name could not be resolved against the backend database."""

    def __init__(self, name: str, suggestions: list[str] | None = None):
        self.name = name
        self.suggestions = suggestions or []
        msg = f"species {name!r} not found in the thermodynamic database"
        if self.suggestions:
            msg += f"; did you mean: {', '.join(self.suggestions)}?"
        super().__init__(msg)


class BalancingError(MicrobialThermoError):
    """A reaction could not be balanced for atoms and charge."""


class AmbiguousReactionError(BalancingError):
    """The reaction has more than one independent balanced solution.

    Carries the solution basis so the caller can choose, rather than the library
    silently picking one. Disproportionation is the usual cause.
    """

    def __init__(self, message: str, basis=None):
        super().__init__(message)
        self.basis = basis


class ThermodynamicConsistencyError(MicrobialThermoError):
    """The two independent free-energy paths disagreed.

    Raised when summing formation energies over the full reaction disagrees with
    -nF*deltaE computed from the two half reactions. Almost always indicates a
    sign-convention slip, a miscounted electron, or a bad half-reaction split.
    """


class OutOfRangeError(MicrobialThermoError):
    """A requested temperature or pressure is outside the supported window."""


class MissingDataError(MicrobialThermoError):
    """Required thermodynamic data is unavailable for a species."""
