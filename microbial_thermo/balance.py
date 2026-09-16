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

from .exceptions import (
    AmbiguousReactionError,
    BalancingError,
    MicrobialThermoError,
)
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

#: "+" immediately followed by a letter, i.e. a separator missing its spaces.
_UNSPACED_SEPARATOR = re.compile(r"\+[A-Za-z]")


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

    #: ``((Species, coefficient), ...)`` for each side. Usually one entry, but
    #: an incomplete oxidation has several -- syntrophic propionate oxidation
    #: yields acetate *and* bicarbonate, in a ratio conservation cannot fix.
    oxidized_side: tuple
    reduced_side: tuple
    key_element: str
    coefficients: dict
    n_electrons: Fraction

    @property
    def species(self) -> list[Species]:
        return [s for s in self.coefficients if s != ELECTRON]

    @property
    def is_simple(self) -> bool:
        """True when each side names a single species, the usual case."""
        return len(self.oxidized_side) == 1 and len(self.reduced_side) == 1

    @property
    def oxidized(self) -> Species:
        """The single oxidized species; raises for a multi-product side."""
        return _sole(self.oxidized_side, "oxidized")

    @property
    def reduced(self) -> Species:
        """The single reduced species; raises for a multi-product side."""
        return _sole(self.reduced_side, "reduced")

    def oxidation_states(self) -> tuple[Fraction, Fraction]:
        """Mean oxidation state of the key element: (oxidized side, reduced side).

        Averaged across the whole side, weighted by how many atoms of the key
        element each species carries. For propionate oxidising to acetate plus
        bicarbonate, the oxidized state is the mean over both products.
        """
        return (
            side_oxidation_state(self.oxidized_side, self.key_element),
            side_oxidation_state(self.reduced_side, self.key_element),
        )

    def scaled(self, factor) -> HalfReaction:
        factor = Fraction(factor)
        return HalfReaction(
            oxidized_side=self.oxidized_side,
            reduced_side=self.reduced_side,
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
            oxidized_side=self.oxidized_side,
            reduced_side=self.reduced_side,
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


def normalize_side(spec, registry=None) -> tuple:
    """Normalise one side of a couple into ``((Species, coefficient), ...)``.

    Accepts a species name, a :class:`Species`, a sequence of names (taken in
    equal proportion), a mapping of name to coefficient, or a sequence of
    ``(name, coefficient)`` pairs.

    Proportions *within* a side are the caller's to state, because conservation
    cannot supply them: propionate could in principle oxidise to acetate plus
    bicarbonate, to one and a half acetate, or to three bicarbonate, and which
    one happens is biochemistry, not stoichiometry. Conservation then fixes the
    scale *between* the two sides.
    """
    registry = registry or default_registry()

    if isinstance(spec, (Species, str)):
        return ((registry.resolve(spec), Fraction(1)),)

    if isinstance(spec, dict):
        items = spec.items()
    else:
        try:
            items = list(spec)
        except TypeError as exc:
            raise BalancingError(f"could not read {spec!r} as one side of a redox couple") from exc
        # A bare sequence of names means equal proportions.
        if items and all(isinstance(i, (str, Species)) for i in items):
            items = [(entry, 1) for entry in items]

    side = []
    for entry in items:
        try:
            name, coefficient = entry
        except (TypeError, ValueError) as exc:
            raise BalancingError(
                f"could not read {entry!r} as a (species, coefficient) pair"
            ) from exc
        value = Fraction(coefficient)
        if value <= 0:
            raise BalancingError(f"coefficient for {name!r} must be positive, got {coefficient}")
        side.append((registry.resolve(name), value))

    if not side:
        raise BalancingError("a redox couple side cannot be empty")
    return tuple(side)


def _sole(side: tuple, which: str) -> Species:
    if len(side) != 1:
        raise BalancingError(
            f"this couple's {which} side names "
            f"{', '.join(s.backend for s, _ in side)}; use .{which}_side for a "
            "multi-product couple"
        )
    return side[0][0]


def side_element_count(side: tuple, element: str) -> Fraction:
    """Total atoms of ``element`` across a side."""
    return sum(
        (coefficient * species.parsed.element_count(element) for species, coefficient in side),
        Fraction(0),
    )


def side_charge(side: tuple) -> Fraction:
    return sum((coefficient * species.charge for species, coefficient in side), Fraction(0))


def side_oxidation_state(side: tuple, element: str) -> Fraction:
    """Mean oxidation state of ``element`` across a side, atom-weighted."""
    total_atoms = side_element_count(side, element)
    if total_atoms == 0:
        raise BalancingError(f"no {element} on the side {', '.join(s.backend for s, _ in side)}")
    weighted = sum(
        (
            coefficient
            * species.parsed.element_count(element)
            * mean_oxidation_state(element, species.formula)
            for species, coefficient in side
            if species.parsed.element_count(element)
        ),
        Fraction(0),
    )
    return weighted / total_atoms


def _side_elements(side: tuple) -> set:
    elements: set = set()
    for species, _ in side:
        elements |= set(species.parsed.elements)
    return elements


def _describe(side: tuple) -> str:
    return " + ".join(f"{_format_coefficient(c)}{s.backend}" for s, c in side).strip()


def find_redox_element(reduced_side: tuple, oxidized_side: tuple) -> str:
    """Identify the element that changes oxidation state between two sides.

    Prefers an element other than hydrogen or oxygen, since those are the
    auxiliary balancing species.
    """
    shared = _side_elements(reduced_side) & _side_elements(oxidized_side)
    if not shared:
        raise BalancingError(
            f"{_describe(reduced_side)!r} and {_describe(oxidized_side)!r} share "
            "no element, so they are not a redox couple"
        )

    candidates = sorted(shared - _AUXILIARY_ELEMENTS) or sorted(shared)
    changed = []
    for element in candidates:
        try:
            if side_oxidation_state(reduced_side, element) != side_oxidation_state(
                oxidized_side, element
            ):
                changed.append(element)
        except Exception:
            continue

    if len(changed) == 1:
        return changed[0]
    if len(changed) > 1:
        raise AmbiguousReactionError(
            f"more than one element changes oxidation state between "
            f"{_describe(reduced_side)!r} and {_describe(oxidized_side)!r}: "
            f"{', '.join(changed)}. Pass key_element explicitly.",
            basis=changed,
        )
    if len(candidates) == 1:
        return candidates[0]
    raise BalancingError(
        f"could not identify the redox-active element between "
        f"{_describe(reduced_side)!r} and {_describe(oxidized_side)!r}; "
        "pass key_element explicitly"
    )


def balance_half_reaction(
    reduced, oxidized, key_element: str | None = None, registry=None
) -> HalfReaction:
    """Balance a redox couple as a reduction: ``oxidized + n e- -> reduced``.

    Balances the key element, then oxygen with water, hydrogen with protons,
    and charge with electrons.

    Either side may name several species -- see :func:`normalize_side` -- which
    is what allows incomplete oxidations such as propionate to acetate plus
    bicarbonate. Proportions within a side come from the caller; the scale
    between the sides is set by conserving the key element.
    """
    registry = registry or default_registry()
    reduced_side = normalize_side(reduced, registry)
    oxidized_side = normalize_side(oxidized, registry)
    element = key_element or find_redox_element(reduced_side, oxidized_side)

    key_in_reduced = side_element_count(reduced_side, element)
    key_in_oxidized = side_element_count(oxidized_side, element)
    if key_in_reduced == 0 or key_in_oxidized == 0:
        raise BalancingError(
            f"key element {element} missing from {_describe(reduced_side)!r} or "
            f"{_describe(oxidized_side)!r}"
        )

    # The reduced side as given fixes the scale; the oxidized side is scaled to
    # conserve the key element across the arrow.
    scale = key_in_reduced / key_in_oxidized
    oxidized_side = tuple((s, c * scale) for s, c in oxidized_side)

    # Written as a reduction:
    #     oxidized side + h H+ + n e-  ->  reduced side + w H2O
    # Negative w or h simply moves that species to the other side. When the
    # couple *is* water or the proton, these terms come out at zero, which is
    # why the O2/H2O and H+/H2 couples need no special case.
    water = side_element_count(oxidized_side, "O") - side_element_count(reduced_side, "O")
    proton = (
        side_element_count(reduced_side, "H") + 2 * water - side_element_count(oxidized_side, "H")
    )
    electrons = side_charge(oxidized_side) + proton - side_charge(reduced_side)

    coefficients: dict = {}
    for species, coefficient in oxidized_side:
        _add(coefficients, species, -coefficient)
    for species, coefficient in reduced_side:
        _add(coefficients, species, coefficient)
    _add(coefficients, _water(registry), water)
    _add(coefficients, _proton(registry), -proton)
    _add(coefficients, ELECTRON, -electrons)

    half = HalfReaction(
        oxidized_side=oxidized_side,
        reduced_side=reduced_side,
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
    # A "+" followed by a letter is a separator someone forgot to space out.
    # A "+" followed by a digit, or ending the token, is a charge: "Fe+3", "NH4+".
    if len(tokens) == 1 and _UNSPACED_SEPARATOR.search(tokens[0]):
        raise BalancingError(
            f"could not separate the species in {side.strip()!r} (from "
            f"{equation!r}); put spaces around the separating '+', as in "
            "'SO4-2 + H+' -- an unspaced '+' is read as a charge"
        )
    return tokens


def _readable_basis(species, nullspace) -> list:
    """Turn sympy nullspace vectors into labelled coefficient sets."""
    out = []
    for vector in nullspace:
        values = [Fraction(int(v.p), int(v.q)) for v in vector]
        if any(v != 0 for v in values):
            first = next(v for v in values if v != 0)
            if first < 0:
                values = [-v for v in values]
        entry = {sp: v for sp, v in zip(species, values, strict=True) if v != 0}
        out.append(clear_denominators(entry))
    return out


def _describe_solution(coefficients: dict) -> str:
    return " , ".join(f"{v}·{s.backend}" for s, v in coefficients.items())


def balance_equation(equation: str, registry=None, fix=None) -> dict:
    """Balance a full reaction string into signed, Species-keyed coefficients.

    Solves the homogeneous system ``A x = 0`` where the rows of ``A`` are the
    conserved quantities (each element, plus charge) and the columns are the
    species. The nullspace is computed in exact rational arithmetic by sympy,
    so coefficients are never floating point.

    The nullspace dimension says how determined the problem is:

    * 0 -- no non-trivial solution; the reaction cannot be balanced.
    * 1 -- a unique answer up to scale, the normal case.
    * >1 -- genuinely underdetermined. Conservation alone does not pick an
      answer, and nor should we.

    For that last case, pass ``fix`` to supply the missing information: a
    mapping of species name to how many moles of it participate, as a positive
    number, e.g. ``fix={"O2": 1}``. Each entry is one extra constraint. Without
    it the raised :class:`AmbiguousReactionError` carries the solution basis in
    readable form, so you can see what the choices actually are.

    ``chempy.balance_stoichiometry`` is deliberately not used here: its
    ``underdetermined=None`` mode solves an integer program through ``pulp``,
    needing an external CBC binary that is not in the environment, and it
    offers no way to inspect the nullspace.
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

    if len(nullspace) > 1 and fix is None:
        options = _readable_basis(species, nullspace)
        listed = "; ".join(f"({i}) {_describe_solution(o)}" for i, o in enumerate(options))
        raise AmbiguousReactionError(
            f"reaction {equation!r} has {len(nullspace)} independent balanced "
            "solutions, so conservation alone does not determine the "
            f"stoichiometry. Any combination of these balances: {listed}. "
            "Say which you mean by passing fix={'SPECIES': moles} -- one entry "
            "per degree of freedom, which sets both the ratio and the scale, so "
            f"{len(nullspace)} here -- or "
            "specify the donor and acceptor couples with "
            "Reaction.from_couples instead.",
            basis=options,
        )

    if len(nullspace) > 1:
        coefficients = _solve_with_constraints(
            species, reactants, rows, fix, registry, equation, len(nullspace)
        )
    else:
        vector = list(nullspace[0])
        if vector[0] > 0:
            vector = [-v for v in vector]
        coefficients = {}
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


def _solve_with_constraints(species, reactants, rows, fix, registry, equation, dimension) -> dict:
    """Re-solve an underdetermined system with the caller's extra constraints.

    ``fix`` gives moles as a positive number; the sign comes from which side of
    the arrow the species was written on, so the caller need not think about it.
    """
    from sympy import Rational, linsolve, symbols

    unknowns = symbols(f"x0:{len(species)}")
    equations = [sum(row[i] * unknowns[i] for i in range(len(species))) for row in rows]

    reactant_backends = {s.backend for s in reactants}
    for name, moles in fix.items():
        target = registry.resolve(name)
        try:
            index = next(i for i, s in enumerate(species) if s.backend == target.backend)
        except StopIteration:
            raise BalancingError(
                f"cannot fix {name!r}: it does not appear in {equation!r}"
            ) from None
        value = Rational(moles)
        if target.backend in reactant_backends:
            value = -value
        equations.append(unknowns[index] - value)

    solution = linsolve(equations, unknowns)
    if not solution:
        raise BalancingError(
            f"no balanced solution for {equation!r} with fix={fix!r}; the "
            "constraint is inconsistent with conservation"
        )

    values = list(next(iter(solution)))
    free = sorted({str(sym) for value in values for sym in value.free_symbols})
    if free:
        raise AmbiguousReactionError(
            f"{equation!r} is still underdetermined with fix={fix!r}: "
            f"{len(free)} degree(s) of freedom remain. A nullspace of dimension "
            f"{dimension} needs {dimension} constraints in all -- they set the "
            "ratio between the independent solutions and the overall scale.",
            basis=None,
        )

    coefficients = {}
    for entry, value in zip(species, values, strict=True):
        fraction = Fraction(int(value.p), int(value.q))
        if fraction != 0:
            _add(coefficients, entry, fraction)
    return coefficients


#: Species that never form a couple of their own -- they are the auxiliaries
#: used to balance oxygen, hydrogen and charge.
_NEVER_A_COUPLE = frozenset({"H+", "H2O", "OH-"})

#: An element in its elemental form couples to a conventional partner that the
#: written equation usually leaves implicit: H2 oxidises to H+, O2 reduces to
#: water.
_IMPLICIT_PARTNERS = {"H2": "H+", "O2": "H2O"}


def infer_couples(equation: str, registry=None) -> tuple:
    """Work out the donor and acceptor couples from a written equation.

    Returns ``(donor_couple, acceptor_couple)`` as ``(reduced, oxidized)``
    name pairs, ready for :meth:`Reaction.from_couples`.

    The method: find every element whose oxidation state differs between a
    reactant and a product, which gives the couples directly. Any redox-active
    species left over -- typically H2 or O2, whose partner the writer left
    implicit -- is paired with that conventional partner. Hydrogen and oxygen
    are considered last, since they are also the auxiliary balancing species
    and would otherwise produce spurious couples: in
    ``NO3- + H2 -> NH2OH`` the hydrogens of hydroxylamine are not an oxidation
    product of H2, they are just hydrogens.

    Raises when the result is not a single donor and a single acceptor, since
    guessing there would be worse than asking.
    """
    registry = registry or default_registry()
    reactant_names, product_names = parse_equation(equation)
    reactants = [registry.resolve(n) for n in reactant_names]
    products = [registry.resolve(n) for n in product_names]

    candidates = [s for s in reactants + products if s.backend not in _NEVER_A_COUPLE]
    elements = sorted({e for s in candidates for e in s.parsed.elements} - _AUXILIARY_ELEMENTS)

    found: list[tuple] = []  # (reactant, product, element)
    assigned: set = set()
    # Species whose oxidation state cannot be read from the formula alone --
    # vanadyl sulfate carries both V and S, neither of them a spectator. These
    # are invisible to the search below, so remember them for the error message.
    unassignable: dict = {}
    for element in elements:
        left = [
            s
            for s in reactants
            if s.parsed.element_count(element) and s.backend not in _NEVER_A_COUPLE
        ]
        right = [
            s
            for s in products
            if s.parsed.element_count(element) and s.backend not in _NEVER_A_COUPLE
        ]
        for r in left:
            for p in right:
                try:
                    if mean_oxidation_state(element, r.formula) == mean_oxidation_state(
                        element, p.formula
                    ):
                        continue
                except MicrobialThermoError:
                    for candidate in (r, p):
                        try:
                            mean_oxidation_state(element, candidate.formula)
                        except MicrobialThermoError:
                            unassignable[candidate.backend] = element
                    continue
                except Exception:
                    continue
                found.append((r, p, element))
                assigned.add(r.backend)
                assigned.add(p.backend)

    # Anything still unassigned should be an elemental species whose partner the
    # writer left out.
    for species in candidates:
        if species.backend in assigned:
            continue
        partner_name = _IMPLICIT_PARTNERS.get(species.formula)
        if partner_name is None:
            continue
        partner = registry.resolve(partner_name)
        element = next(iter(species.parsed.elements))
        on_left = any(s.backend == species.backend for s in reactants)
        if on_left:
            found.append((species, partner, element))
        else:
            found.append((partner, species, element))
        assigned.add(species.backend)

    donors, acceptors = [], []
    for reactant, product, element in found:
        reactant_state = mean_oxidation_state(element, reactant.formula)
        product_state = mean_oxidation_state(element, product.formula)
        if reactant_state < product_state:
            # The reactant loses electrons: it is the donor, written oxidatively.
            donors.append((reactant, product, element))
        else:
            acceptors.append((product, reactant, element))

    if len(donors) != 1 or len(acceptors) != 1:
        message = (
            f"could not read a single donor and a single acceptor from "
            f"{equation!r}: found {len(donors)} donor(s) "
            f"[{', '.join(f'{r.backend}->{p.backend}' for r, p, _ in donors)}] and "
            f"{len(acceptors)} acceptor(s) "
            f"[{', '.join(f'{o.backend}->{r.backend}' for r, o, _ in acceptors)}]."
        )
        if unassignable:
            listed = ", ".join(
                f"{name} ({element})" for name, element in sorted(unassignable.items())
            )
            message += (
                f" The oxidation state of {listed} cannot be read from the formula "
                "alone, because more than one of its elements has no conventional "
                "state, so it could not take part in a couple here."
            )
        message += (
            " Build it with Reaction.from_couples(donor=..., acceptor=...) instead, "
            "naming the couples yourself."
        )
        raise AmbiguousReactionError(message, basis=(donors, acceptors))

    donor_reduced, donor_oxidized, _ = donors[0]
    acceptor_reduced, acceptor_oxidized, _ = acceptors[0]
    return (
        (donor_reduced.backend, donor_oxidized.backend),
        (acceptor_reduced.backend, acceptor_oxidized.backend),
    )
