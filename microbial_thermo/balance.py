"""Balancing reactions and half reactions.

Half reactions are balanced by the standard biogeochemical procedure: balance
the redox-active element first, then oxygen with water, then hydrogen with
protons, then charge with electrons. Doing this explicitly rather than through a
general nullspace keeps the result deterministic and keeps each step available
for show-your-work output.

Everything here operates on :class:`~microbial_thermo.species.Species` objects
rather than bare strings, because a pyGCC database name ("Methane(aq)") is not a
parseable formula. The species carries both identities.

Every coefficient is a :class:`~fractions.Fraction`; floats are never used for
stoichiometry, so normalising to an electron pair gives exact halves rather
than 0.4999999.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction
from math import gcd

from .exceptions import AmbiguousReactionError, BalancingError
from .oxidation import mean_oxidation_state
from .species import Species, default_registry

#: The free electron. Not a database species anywhere, but a first-class
#: participant in half reactions.
ELECTRON = Species(backend="e-", formula="e-", display="e^-", phase="aq")

#: Elements balanced by the auxiliary species rather than treated as redox-active.
_AUXILIARY_ELEMENTS = frozenset({"H", "O"})

_ARROW = re.compile(r"\s*(?:->|=>|-->|=|→|⟶)\s*")

#: A separating "+" must have whitespace on both sides; a bare "+" is a charge.
_PLUS_SEPARATOR = re.compile(r"\s+\+\s+")


def _water(registry=None) -> Species:
    return (registry or default_registry()).resolve("H2O")


def _proton(registry=None) -> Species:
    return (registry or default_registry()).resolve("H+")


def _lcm(a: int, b: int) -> int:
    return abs(a * b) // gcd(a, b) if a and b else abs(a or b)


def clear_denominators(coefficients: dict) -> dict:
    """Scale a coefficient set to the smallest integers."""
    denominator = 1
    for value in coefficients.values():
        denominator = _lcm(denominator, value.denominator)
    scaled = {k: v * denominator for k, v in coefficients.items()}
    numerators = [abs(v.numerator) for v in scaled.values() if v != 0]
    if not numerators:
        return scaled
    common = 0
    for n in numerators:
        common = gcd(common, n)
    if common > 1:
        scaled = {k: v / common for k, v in scaled.items()}
    return scaled


def _format_coefficient(value: Fraction) -> str:
    if value == 1:
        return ""
    if value.denominator == 1:
        return f"{value.numerator} "
    return f"{value.numerator}/{value.denominator} "


def _format_side(terms) -> str:
    return " + ".join(f"{_format_coefficient(v)}{s.backend}" for s, v in terms) or "0"


@dataclass(frozen=True)
class HalfReaction:
    """A half reaction stored canonically as a *reduction*.

    ``coefficients`` is signed: negative for reactants, positive for products,
    so ``sum(nu_i * dGf_i)`` is a direct dot product. Storing it as a reduction
    regardless of display direction keeps every potential comparable, which is
    what makes a redox tower meaningful.
    """

    oxidized: Species
    reduced: Species
    key_element: str
    coefficients: dict
    n_electrons: Fraction

    @property
    def species(self) -> list[Species]:
        return [s for s in self.coefficients if s != ELECTRON]

    def oxidation_states(self) -> tuple[Fraction, Fraction]:
        """Mean oxidation state of the key element: (oxidized form, reduced form)."""
        return (
            mean_oxidation_state(self.key_element, self.oxidized.formula),
            mean_oxidation_state(self.key_element, self.reduced.formula),
        )

    def scaled(self, factor) -> HalfReaction:
        factor = Fraction(factor)
        return HalfReaction(
            oxidized=self.oxidized,
            reduced=self.reduced,
            key_element=self.key_element,
            coefficients={k: v * factor for k, v in self.coefficients.items()},
            n_electrons=self.n_electrons * factor,
        )

    def normalized_to_electrons(self, n: int = 2) -> HalfReaction:
        """Rescale so exactly ``n`` electrons are transferred."""
        if self.n_electrons == 0:
            raise BalancingError("cannot normalise a half reaction with no electrons")
        return self.scaled(Fraction(n) / self.n_electrons)

    def normalized_to(self, species) -> HalfReaction:
        """Rescale so one mole of ``species`` participates."""
        target = default_registry().resolve(species)
        for entry, value in self.coefficients.items():
            if entry.backend == target.backend and value != 0:
                return self.scaled(Fraction(1) / abs(value))
        raise BalancingError(f"{target.backend} does not appear in half reaction {self.format()}")

    def reversed(self) -> HalfReaction:
        return HalfReaction(
            oxidized=self.oxidized,
            reduced=self.reduced,
            key_element=self.key_element,
            coefficients={k: -v for k, v in self.coefficients.items()},
            n_electrons=-self.n_electrons,
        )

    def format(self, direction: str = "reduction") -> str:
        coefficients = (
            self.coefficients
            if direction == "reduction"
            else {k: -v for k, v in self.coefficients.items()}
        )
        left = [(s, -v) for s, v in coefficients.items() if v < 0]
        right = [(s, v) for s, v in coefficients.items() if v > 0]
        return f"{_format_side(left)} -> {_format_side(right)}"

    def __str__(self) -> str:
        return self.format()


def find_redox_element(reduced: Species, oxidized: Species) -> str:
    """Identify the element that changes oxidation state between two forms.

    Prefers an element other than hydrogen or oxygen, since those are the
    auxiliary balancing species.
    """
    red, ox = reduced.parsed, oxidized.parsed
    shared = set(red.elements) & set(ox.elements)
    if not shared:
        raise BalancingError(
            f"{reduced.backend!r} and {oxidized.backend!r} share no element, "
            "so they are not a redox couple"
        )

    candidates = sorted(shared - _AUXILIARY_ELEMENTS) or sorted(shared)
    changed = []
    for element in candidates:
        try:
            if mean_oxidation_state(element, reduced.formula) != mean_oxidation_state(
                element, oxidized.formula
            ):
                changed.append(element)
        except Exception:
            continue

    if len(changed) == 1:
        return changed[0]
    if len(changed) > 1:
        raise AmbiguousReactionError(
            f"more than one element changes oxidation state between "
            f"{reduced.backend!r} and {oxidized.backend!r}: {', '.join(changed)}. "
            "Pass key_element explicitly.",
            basis=changed,
        )
    if len(candidates) == 1:
        return candidates[0]
    raise BalancingError(
        f"could not identify the redox-active element between "
        f"{reduced.backend!r} and {oxidized.backend!r}; pass key_element explicitly"
    )


def balance_half_reaction(
    reduced, oxidized, key_element: str | None = None, registry=None
) -> HalfReaction:
    """Balance a redox couple as a reduction: ``oxidized + n e- -> reduced``.

    Balances the key element, then oxygen with water, hydrogen with protons,
    and charge with electrons.
    """
    registry = registry or default_registry()
    reduced = registry.resolve(reduced)
    oxidized = registry.resolve(oxidized)
    element = key_element or find_redox_element(reduced, oxidized)

    red, ox = reduced.parsed, oxidized.parsed
    n_red = red.element_count(element)
    n_ox = ox.element_count(element)
    if n_red == 0 or n_ox == 0:
        raise BalancingError(
            f"key element {element} missing from {reduced.backend!r} or {oxidized.backend!r}"
        )

    # One mole of the reduced form fixes the scale.
    a_red = Fraction(1)
    a_ox = Fraction(n_red, n_ox)

    # Written as a reduction:
    #     a_ox Oxidized + h H+ + n e-  ->  a_red Reduced + w H2O
    # Negative w or h simply moves that species to the other side. When the
    # couple *is* water or the proton, these terms come out at zero, which is
    # why the O2/H2O and H+/H2 couples need no special case.
    water = a_ox * ox.element_count("O") - a_red * red.element_count("O")
    proton = a_red * red.element_count("H") + 2 * water - a_ox * ox.element_count("H")
    electrons = a_ox * ox.charge + proton - a_red * red.charge

    coefficients: dict = {}
    _add(coefficients, oxidized, -a_ox)
    _add(coefficients, reduced, a_red)
    _add(coefficients, _water(registry), water)
    _add(coefficients, _proton(registry), -proton)
    _add(coefficients, ELECTRON, -electrons)

    half = HalfReaction(
        oxidized=oxidized,
        reduced=reduced,
        key_element=element,
        coefficients={k: v for k, v in coefficients.items() if v != 0},
        n_electrons=electrons,
    )
    verify_conservation(half.coefficients, label=half.format())
    return half


def _add(target: dict, species: Species, value: Fraction) -> None:
    target[species] = target.get(species, Fraction(0)) + value


def verify_conservation(coefficients: dict, label: str = "reaction") -> None:
    """Assert atom and charge conservation. A failure here is a library bug."""
    elements: dict[str, Fraction] = {}
    charge = Fraction(0)
    for species, nu in coefficients.items():
        parsed = species.parsed
        charge += nu * parsed.charge
        for symbol, count in parsed.elements.items():
            elements[symbol] = elements.get(symbol, Fraction(0)) + nu * count
    bad = {k: v for k, v in elements.items() if v != 0}
    if bad:
        raise BalancingError(f"{label} does not conserve atoms: {bad}")
    if charge != 0:
        raise BalancingError(f"{label} does not conserve charge (net {charge})")


def combine_half_reactions(
    donor_half: HalfReaction, acceptor_half: HalfReaction, n_electrons: int = 2
) -> dict:
    """Combine two half reactions into a full reaction.

    The donor couple is reversed into its oxidative direction, the acceptor
    stays reductive, and both are scaled to the same electron count so the
    electrons cancel exactly.
    """
    donor = donor_half.normalized_to_electrons(n_electrons).reversed()
    acceptor = acceptor_half.normalized_to_electrons(n_electrons)

    combined: dict = {}
    for source in (donor, acceptor):
        for species, value in source.coefficients.items():
            _add(combined, species, value)

    residual_electrons = combined.pop(ELECTRON, Fraction(0))
    if residual_electrons != 0:
        raise BalancingError(
            f"electrons did not cancel when combining half reactions "
            f"(residual {residual_electrons}); this is a library bug"
        )
    combined = {k: v for k, v in combined.items() if v != 0}
    verify_conservation(combined, label="combined reaction")
    return combined


def parse_equation(equation: str) -> tuple[list[str], list[str]]:
    """Split ``A + B -> C + D`` into reactant and product name lists."""
    parts = _ARROW.split(equation)
    if len(parts) != 2:
        raise BalancingError(
            f"could not find a single reaction arrow in {equation!r}; "
            "use '->' between reactants and products"
        )
    left = _split_side(parts[0], equation)
    right = _split_side(parts[1], equation)
    if not left or not right:
        raise BalancingError(f"reaction {equation!r} has an empty side")
    return left, right


def _split_side(side: str, equation: str) -> list[str]:
    """Split one side of an equation on its separating plus signs.

    The separator must be surrounded by whitespace, because a bare ``+`` is
    also the charge marker: splitting ``HS- + H+`` on every ``+`` would yield a
    species called ``H``.
    """
    tokens = [t.strip() for t in _PLUS_SEPARATOR.split(side) if t.strip()]
    if len(tokens) == 1 and "+" in tokens[0][:-1]:
        raise BalancingError(
            f"could not separate the species in {side.strip()!r} (from "
            f"{equation!r}); put spaces around the separating '+', as in "
            "'SO4-2 + H+' -- an unspaced '+' is read as a charge"
        )
    return tokens


def balance_equation(equation: str, registry=None) -> dict:
    """Balance a full reaction string into signed, Species-keyed coefficients.

    Solves the homogeneous system ``A x = 0`` where the rows of ``A`` are the
    conserved quantities (each element, plus charge) and the columns are the
    species. The nullspace is computed in exact rational arithmetic by sympy,
    so coefficients are never floating point.

    The nullspace dimension tells us how determined the problem is:

    * 0 -- no non-trivial solution; the reaction cannot be balanced.
    * 1 -- a unique answer up to scale, which is the normal case.
    * >1 -- genuinely ambiguous, disproportionation being the usual cause. We
      raise :class:`AmbiguousReactionError` carrying the basis rather than
      silently returning one of infinitely many answers.

    We deliberately do not use ``chempy.balance_stoichiometry`` here: its
    ``underdetermined=None`` mode solves an integer program through ``pulp``,
    which requires an external CBC solver binary that is not part of the conda
    environment, and it offers no way to inspect the nullspace dimension.
    """
    from sympy import Matrix, Rational

    registry = registry or default_registry()
    reactant_names, product_names = parse_equation(equation)
    reactants = [registry.resolve(n) for n in reactant_names]
    products = [registry.resolve(n) for n in product_names]
    species = reactants + products
    if len(set(species)) != len(species):
        raise BalancingError(
            f"a species appears more than once in {equation!r}; "
            "each species may appear on only one side"
        )

    conserved = sorted({e for s in species for e in s.parsed.elements})
    rows = [[Rational(s.parsed.element_count(e)) for s in species] for e in conserved]
    rows.append([Rational(s.parsed.charge) for s in species])

    nullspace = Matrix(rows).nullspace()
    if not nullspace:
        raise BalancingError(
            f"{equation!r} cannot be balanced: no combination of these species "
            "conserves both atoms and charge. Check for a missing species such "
            "as H2O, H+, or an electron acceptor."
        )
    if len(nullspace) > 1:
        raise AmbiguousReactionError(
            f"reaction {equation!r} has {len(nullspace)} independent balanced "
            "solutions, so the stoichiometry is not determined by conservation "
            "alone (disproportionation is the usual cause). Specify the donor "
            "and acceptor couples explicitly instead.",
            basis=nullspace,
        )

    vector = list(nullspace[0])
    # Orient the solution so the declared reactants come out negative.
    if vector[0] > 0:
        vector = [-v for v in vector]

    coefficients: dict = {}
    for entry, value in zip(species, vector, strict=True):
        _add(coefficients, entry, Fraction(int(value.p), int(value.q)))
    coefficients = clear_denominators({k: v for k, v in coefficients.items() if v != 0})

    misplaced = [s.backend for s in reactants if coefficients.get(s, Fraction(0)) > 0]
    misplaced += [s.backend for s in products if coefficients.get(s, Fraction(0)) < 0]
    if misplaced:
        raise BalancingError(
            f"{equation!r} only balances if {', '.join(sorted(set(misplaced)))} "
            "switch sides; check which species are reactants and which are products"
        )

    verify_conservation(coefficients, label=f"balanced {equation!r}")
    return coefficients
