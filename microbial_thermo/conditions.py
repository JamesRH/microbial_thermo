"""Environmental conditions for a thermodynamic calculation."""

from __future__ import annotations

from dataclasses import dataclass, field

from .units import STANDARD_BIOCHEMICAL_PH

#: Activity models available for aqueous species.
ACTIVITY_MODELS = ("bdot", "ideal", "unit")


@dataclass
class Conditions:
    """Temperature, pH, and composition at which a reaction is evaluated.

    ``concentrations`` are molalities keyed by any name the species registry
    resolves; ``partial_pressures`` are in bar for gaseous species. Anything not
    listed is held at unit activity, which reproduces the standard state.

    ``total_concentrations`` are keyed by an acid-base family instead -- pass
    ``{"sulfide": 1e-6}`` or ``{"DIC": 2.2e-3}`` -- and the pH-dependent
    fraction belonging to each member is worked out for you. This is usually
    what you want, because measurements are reported as totals while reactions
    are written with one specific form. A per-species entry in
    ``concentrations`` takes precedence over a family total.

    ``activity_model``:

    ``"bdot"``
        Extended Debye-Huckel with pyGCC's solvent parameters. Requires an
        ionic strength.
    ``"ideal"``
        Activity equals concentration; activity coefficients are 1.
    ``"unit"``
        Every activity is 1 regardless of concentration, i.e. the standard
        state. Used internally for the standard-state quantities.
    """

    temperature_c: float = 25.0
    pH: float = STANDARD_BIOCHEMICAL_PH
    pressure_bar: float | None = None
    ionic_strength: float = 0.0
    concentrations: dict[str, float] = field(default_factory=dict)
    total_concentrations: dict[str, float] = field(default_factory=dict)
    partial_pressures: dict[str, float] = field(default_factory=dict)
    activity_model: str = "ideal"
    water_activity: float = 1.0

    def __post_init__(self):
        if self.activity_model not in ACTIVITY_MODELS:
            raise ValueError(
                f"activity_model must be one of {ACTIVITY_MODELS}, got {self.activity_model!r}"
            )
        if self.activity_model == "bdot" and self.ionic_strength <= 0:
            raise ValueError(
                "the 'bdot' activity model needs a positive ionic_strength; "
                "pass ionic_strength= or use activity_model='ideal'"
            )

    @property
    def proton_activity(self) -> float:
        return 10.0**-self.pH

    @property
    def uses_speciation(self) -> bool:
        return bool(self.total_concentrations)

    def replace(self, **changes) -> Conditions:
        """Return a copy with ``changes`` applied."""
        from dataclasses import replace as _replace

        return _replace(self, **changes)

    def standard(self) -> Conditions:
        """The same T and P, but every activity at unity ([H+] = 1 M)."""
        return Conditions(
            temperature_c=self.temperature_c,
            pH=0.0,
            pressure_bar=self.pressure_bar,
            activity_model="unit",
        )

    def standard_prime(self, pH: float | None = None) -> Conditions:
        """The same T and P, unit activities, but protons at the given pH."""
        return Conditions(
            temperature_c=self.temperature_c,
            pH=self.pH if pH is None else pH,
            pressure_bar=self.pressure_bar,
            activity_model="unit",
        )
