"""Redox couples, half reactions, and whole reactions with their energetics.

The central objects are :class:`Couple` (a reduced/oxidized pair),
:class:`HalfReactionResult` (a balanced half reaction plus its potentials), and
:class:`Reaction` (a donor couple driving an acceptor couple).

Sign conventions, stated once and enforced everywhere:

* Half reactions are stored as **reductions**, per IUPAC.
* :math:`E = -\\Delta G_\\mathrm{half} / (nF)`.
* For the whole reaction, :math:`\\Delta E = E_\\mathrm{acceptor} - E_\\mathrm{donor}`
  and :math:`\\Delta G = -nF\\Delta E`.
* :math:`\\Delta G_f^\\circ(\\mathrm{H}^+) = 0` and :math:`\\Delta G_f^\\circ(e^-) = 0`,
  which places everything on the standard hydrogen electrode.

Every :class:`Reaction` verifies its free energy along two independent paths and
raises if they disagree. Redox sign errors are the single most common mistake in
this domain -- textbooks included -- so the check runs on every calculation
rather than only in the test suite.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from .balance import (
    ELECTRON,
    HalfReaction,
    balance_half_reaction,
    combine_half_reactions,
)
from .conditions import Conditions
from .exceptions import ThermodynamicConsistencyError
from .species import Species, default_registry
from .units import (
    FARADAY,
    KJ_PER_MOL_STR,
    Quantity,
    R,
    as_magnitude,
    celsius_to_kelvin,
    ureg,
)

#: Tolerance for the two-path free-energy cross-check, in kJ/mol.
CONSISTENCY_TOLERANCE_KJ = 1e-6


@dataclass(frozen=True)
class Couple:
    """A redox couple, always written reduced/oxidized."""

    reduced: Species
    oxidized: Species
    key_element: str | None = None

    @classmethod
    def make(cls, reduced, oxidized, key_element=None, registry=None) -> Couple:
        registry = registry or default_registry()
        return cls(
            reduced=registry.resolve(reduced),
            oxidized=registry.resolve(oxidized),
            key_element=key_element,
        )

    def half_reaction(self, registry=None) -> HalfReaction:
        return balance_half_reaction(
            self.reduced, self.oxidized, self.key_element, registry=registry
        )

    @property
    def label(self) -> str:
        return f"{self.oxidized.label}/{self.reduced.label}"

    def __str__(self) -> str:
        return f"{self.oxidized.backend}/{self.reduced.backend}"


def _activity(species: Species, conditions: Conditions, backend) -> float:
    """Activity of ``species`` under ``conditions``.

    Protons always follow the pH, in every activity model -- pH is the whole
    point of the primed convention, so it is never treated as "unit activity".
    Water takes ``conditions.water_activity``.
    """
    if species.backend == "H+":
        return conditions.proton_activity
    if species.backend == "H2O":
        return conditions.water_activity
    if species.is_solid:
        return 1.0

    registry = default_registry()

    def lookup(table):
        for name, value in table.items():
            resolved = registry.get(name)
            if resolved is not None and resolved.backend == species.backend:
                return value
        return None

    if species.is_gas:
        pressure = lookup(conditions.partial_pressures)
        return 1.0 if pressure is None else float(pressure)

    if conditions.activity_model == "unit":
        return 1.0

    molality = lookup(conditions.concentrations)
    if molality is None:
        molality = _molality_from_family_total(species, conditions, backend)
        if molality is None:
            return 1.0
    else:
        # An explicit per-species concentration near a pKa is exactly where
        # naming one form of a family goes wrong, so say so.
        from .speciation import warn_if_ambiguous

        warn_if_ambiguous(
            species,
            conditions.pH,
            conditions.temperature_c,
            conditions.pressure_bar,
            backend,
        )
    molality = float(molality)

    if conditions.activity_model == "ideal":
        return molality

    gamma = backend.activity_coefficient(
        charge=species.charge,
        ionic_strength=conditions.ionic_strength,
        temperature_c=conditions.temperature_c,
        pressure_bar=conditions.pressure_bar,
        species_name=species.backend,
    )
    return gamma * molality


def _molality_from_family_total(species: Species, conditions: Conditions, backend):
    """Molality of ``species`` implied by a family total, or None.

    Turns a measured total -- total sulfide, DIC, total ammonia -- into the
    molality of the one form the reaction is written with, using the
    pH-dependent distribution across the family.
    """
    if not conditions.total_concentrations:
        return None

    from .speciation import family_by_name, family_for, fraction_of

    family = family_for(species)
    if family is None:
        return None
    for key, total in conditions.total_concentrations.items():
        if family_by_name(key) is family:
            fraction = fraction_of(
                species,
                conditions.pH,
                conditions.temperature_c,
                conditions.pressure_bar,
                backend,
            )
            return float(total) * fraction
    return None


def standard_gibbs(coefficients: dict, conditions: Conditions, backend) -> Quantity:
    """Sum of ``nu_i * dGf_i`` over a signed coefficient set, in kJ/mol."""
    total = 0.0
    for species, nu in coefficients.items():
        if species == ELECTRON:
            continue  # dGf(e-) = 0 by the SHE convention
        value = backend.delta_Gf(species.backend, conditions.temperature_c, conditions.pressure_bar)
        total += float(nu) * as_magnitude(value, KJ_PER_MOL_STR)
    return Quantity(total, KJ_PER_MOL_STR)


def log_reaction_quotient(coefficients: dict, conditions: Conditions, backend) -> float:
    """``ln Q`` for a signed coefficient set."""
    total = 0.0
    for species, nu in coefficients.items():
        if species == ELECTRON:
            continue
        activity = _activity(species, conditions, backend)
        if activity <= 0:
            raise ValueError(
                f"activity of {species.backend} is {activity}; concentrations must be positive"
            )
        total += float(nu) * math.log(activity)
    return total


def gibbs(coefficients: dict, conditions: Conditions, backend) -> Quantity:
    """``dG = dG_standard + RT ln Q``, in kJ/mol."""
    standard = standard_gibbs(coefficients, conditions, backend)
    rt = (R * (celsius_to_kelvin(conditions.temperature_c) * ureg.kelvin)).to(KJ_PER_MOL_STR)
    correction = rt.magnitude * log_reaction_quotient(coefficients, conditions, backend)
    return Quantity(standard.magnitude + correction, KJ_PER_MOL_STR)


def potential(half: HalfReaction, conditions: Conditions, backend) -> Quantity:
    """``E = -dG_half / (nF)`` for a half reaction, in volts."""
    if half.n_electrons == 0:
        raise ValueError("half reaction transfers no electrons; E is undefined")
    delta_g = gibbs(half.coefficients, conditions, backend)
    n = float(half.n_electrons)
    volts = -delta_g.to(KJ_PER_MOL_STR) / (n * FARADAY.to("C/mol"))
    return volts.to("V")


@dataclass
class HalfReactionResult:
    """A balanced half reaction together with its potentials."""

    couple: Couple
    half: HalfReaction
    conditions: Conditions
    backend: object
    role: str = ""  # "donor" or "acceptor"

    @property
    def n_electrons(self) -> Fraction:
        return self.half.n_electrons

    @property
    def E_standard(self) -> Quantity:
        """E at unit activity of every species, including [H+] = 1 M."""
        return potential(self.half, self.conditions.standard(), self.backend)

    @property
    def E_standard_prime(self) -> Quantity:
        """E at unit activity but protons at the conditions' pH."""
        return potential(self.half, self.conditions.standard_prime(), self.backend)

    @property
    def E(self) -> Quantity:
        """E under the full conditions, including concentrations."""
        return potential(self.half, self.conditions, self.backend)

    def oxidation_states(self):
        return self.half.oxidation_states()

    def format(self, direction: str | None = None) -> str:
        if direction is None:
            direction = "oxidation" if self.role == "donor" else "reduction"
        return self.half.format(direction)

    def __str__(self) -> str:
        return self.format()


@dataclass
class Reaction:
    """A donor couple driving an acceptor couple.

    Build with :meth:`from_couples`. The reaction is normalised to
    ``n_electrons`` transferred, defaulting to an electron pair.
    """

    donor: Couple
    acceptor: Couple
    conditions: Conditions
    backend: object
    n_electrons: int = 2
    coefficients: dict = None
    donor_half: HalfReactionResult = None
    acceptor_half: HalfReactionResult = None

    @classmethod
    def from_couples(
        cls,
        donor,
        acceptor,
        conditions: Conditions | None = None,
        backend=None,
        n_electrons: int = 2,
        normalize_to=None,
        registry=None,
    ) -> Reaction:
        """Construct from two couples, each given as ``(reduced, oxidized)``.

        ``normalize_to`` overrides the electron-pair default: pass ``"donor"``,
        ``"acceptor"``, ``"electron"``, ``"electron_pair"``, or any species name
        to scale the reaction to one mole of that component.
        """
        from . import get_backend

        registry = registry or default_registry()
        backend = backend or get_backend()
        conditions = conditions or Conditions()

        donor_couple = (
            donor if isinstance(donor, Couple) else Couple.make(*donor, registry=registry)
        )
        acceptor_couple = (
            acceptor if isinstance(acceptor, Couple) else Couple.make(*acceptor, registry=registry)
        )

        donor_half = donor_couple.half_reaction(registry)
        acceptor_half = acceptor_couple.half_reaction(registry)

        n = _resolve_normalization(normalize_to, n_electrons, donor_half, acceptor_half, registry)
        coefficients = combine_half_reactions(donor_half, acceptor_half, n_electrons=n)

        reaction = cls(
            donor=donor_couple,
            acceptor=acceptor_couple,
            conditions=conditions,
            backend=backend,
            n_electrons=n,
            coefficients=coefficients,
            donor_half=HalfReactionResult(
                couple=donor_couple,
                half=donor_half.normalized_to_electrons(n),
                conditions=conditions,
                backend=backend,
                role="donor",
            ),
            acceptor_half=HalfReactionResult(
                couple=acceptor_couple,
                half=acceptor_half.normalized_to_electrons(n),
                conditions=conditions,
                backend=backend,
                role="acceptor",
            ),
        )
        reaction.verify_consistency()
        return reaction

    def renormalized(self, n_electrons: int = 2, normalize_to=None) -> Reaction:
        """Return the same chemistry scaled to a different electron count.

        Potentials and free energy per electron are intensive and do not
        change; only the stoichiometry and the total free energy scale.
        """
        return Reaction.from_couples(
            donor=self.donor,
            acceptor=self.acceptor,
            conditions=self.conditions,
            backend=self.backend,
            n_electrons=n_electrons,
            normalize_to=normalize_to,
        )

    # --- energetics ------------------------------------------------------------

    @property
    def delta_G_standard(self) -> Quantity:
        """dG at unit activity of everything, including [H+] = 1 M."""
        return standard_gibbs(self.coefficients, self.conditions.standard(), self.backend)

    @property
    def delta_G_standard_prime(self) -> Quantity:
        """dG at unit activity but protons at the conditions' pH."""
        return gibbs(self.coefficients, self.conditions.standard_prime(), self.backend)

    @property
    def delta_G(self) -> Quantity:
        """dG under the full conditions."""
        return gibbs(self.coefficients, self.conditions, self.backend)

    @property
    def delta_E(self) -> Quantity:
        """E_acceptor - E_donor under the full conditions."""
        return self.acceptor_half.E - self.donor_half.E

    @property
    def delta_E_standard_prime(self) -> Quantity:
        return self.acceptor_half.E_standard_prime - self.donor_half.E_standard_prime

    @property
    def delta_G_per_electron(self) -> Quantity:
        """dG normalised per mole of electrons transferred."""
        return self.delta_G / float(self.n_electrons)

    # --- the two-path cross-check ---------------------------------------------

    def verify_consistency(self, tolerance_kj: float = CONSISTENCY_TOLERANCE_KJ) -> None:
        """Check dG from formation energies against dG from the two potentials.

        Path A sums ``nu_i * dGf_i`` over the combined reaction. Path B computes
        ``-nF(E_acceptor - E_donor)`` from the half reactions. They are
        algebraically equivalent, so any disagreement means the half reactions
        do not actually sum to the combined reaction -- a bad split, a
        miscounted electron, or a sign error.
        """
        path_a = as_magnitude(self.delta_G, KJ_PER_MOL_STR)
        n = float(self.n_electrons)
        path_b_q = -(n * FARADAY.to("C/mol") * self.delta_E.to("V")).to(KJ_PER_MOL_STR)
        path_b = as_magnitude(path_b_q, KJ_PER_MOL_STR)
        if abs(path_a - path_b) > tolerance_kj:
            raise ThermodynamicConsistencyError(
                f"free energy disagrees between the two paths for {self}: "
                f"sum(nu*dGf) = {path_a:.6f} kJ/mol but -nF*dE = {path_b:.6f} "
                f"kJ/mol (difference {path_a - path_b:.2e}). The half reactions "
                "do not sum to the combined reaction."
            )

    # --- presentation ----------------------------------------------------------

    def format(self) -> str:
        left = [(s, -v) for s, v in self.coefficients.items() if v < 0]
        right = [(s, v) for s, v in self.coefficients.items() if v > 0]
        from .balance import _format_side

        return f"{_format_side(left)} -> {_format_side(right)}"

    def __str__(self) -> str:
        return self.format()

    def show_work(self):
        """A step-by-step derivation of every number in this reaction.

        Renders as Markdown with LaTeX in a notebook, as plain text elsewhere.
        """
        from .showwork import derivation

        return derivation(self)

    def summary(self) -> str:
        """A short multi-line report, useful at a notebook prompt."""
        lines = [
            f"Reaction:  {self.format()}",
            f"  donor     {self.donor}  (oxidative)  "
            f"E0' = {self.donor_half.E_standard_prime.to('V').magnitude:+.3f} V",
            f"  acceptor  {self.acceptor}  (reductive)  "
            f"E0' = {self.acceptor_half.E_standard_prime.to('V').magnitude:+.3f} V",
            f"  n electrons        {self.n_electrons}",
            f"  dE0'               {self.delta_E_standard_prime.to('V').magnitude:+.3f} V",
            f"  dG0'               {self.delta_G_standard_prime.magnitude:+.1f} kJ/mol",
            f"  dG (conditions)    {self.delta_G.magnitude:+.1f} kJ/mol",
            f"  dG per electron    {self.delta_G_per_electron.magnitude:+.1f} kJ/mol e-",
        ]
        return "\n".join(lines)


def _resolve_normalization(normalize_to, n_electrons, donor_half, acceptor_half, registry) -> int:
    """Map a ``normalize_to`` request onto an electron count."""
    if normalize_to is None:
        return int(n_electrons)
    if normalize_to == "electron_pair":
        return 2
    if normalize_to == "electron":
        return 1
    if normalize_to == "donor":
        return int(abs(donor_half.n_electrons))
    if normalize_to == "acceptor":
        return int(abs(acceptor_half.n_electrons))
    if normalize_to in ("integer", "integers", "whole"):
        return _smallest_integer_electron_count(donor_half, acceptor_half)

    # Any other value names a species: scale so one mole of it participates.
    target = registry.resolve(normalize_to)
    for half in (donor_half, acceptor_half):
        for species, value in half.coefficients.items():
            if species.backend == target.backend and value != 0:
                n = abs(half.n_electrons / value)
                return int(n) if float(n).is_integer() else n
    raise ValueError(f"cannot normalise to {normalize_to!r}: it appears in neither half reaction")


def _smallest_integer_electron_count(donor_half, acceptor_half) -> int:
    """Smallest electron count that clears every fractional coefficient.

    Normalising to an electron pair is the default because it makes couples
    directly comparable, but it often leaves coefficients in fifths or
    sevenths, which reads badly on a teaching figure. This finds the smallest
    n for which both half reactions come out in whole numbers -- the form a
    textbook would print.
    """
    from math import gcd

    denominator = 1
    for half in (donor_half, acceptor_half):
        n = abs(half.n_electrons)
        if n == 0:
            continue
        for value in half.coefficients.values():
            # Coefficient per electron; its denominator must divide n.
            per_electron = Fraction(value) / n
            denominator = (
                denominator * per_electron.denominator // gcd(denominator, per_electron.denominator)
            )
    return int(denominator)
