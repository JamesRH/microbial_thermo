"""Formal oxidation states and the nominal oxidation state of carbon (NOSC).

Two independent routes are provided:

``mean_oxidation_state``
    Formula-only. Solves for the mean state of one element by assigning
    conventional states to the spectator elements. This is what the
    half-reaction figure displays, since the specification calls for the
    *average* state when a species holds several atoms of the element.

``atom_oxidation_states``
    Per-atom, from a SMILES via RDKit's bond graph, using electronegativity
    partitioning: each bond's electrons are assigned to the more electronegative
    partner and homonuclear bonds are split evenly.

RDKit's own ``rdMolDescriptors.CalcOxidationNumbers`` is documented as
experimental and was observed to violate charge conservation on neutral ethanol
(states summed to +6 rather than 0), so it is deliberately not used. Our
implementation enforces the conservation invariant in code rather than trusting
it.
"""

from __future__ import annotations

from fractions import Fraction

from .exceptions import MicrobialThermoError
from .formula import parse_formula

#: Pauling electronegativities for the elements this library encounters.
#: Hard-coded rather than pulled from a dependency so the values used in a
#: teaching tool are auditable in one place.
PAULING_ELECTRONEGATIVITY: dict[str, float] = {
    "H": 2.20,
    "Li": 0.98,
    "Be": 1.57,
    "B": 2.04,
    "C": 2.55,
    "N": 3.04,
    "O": 3.44,
    "F": 3.98,
    "Na": 0.93,
    "Mg": 1.31,
    "Al": 1.61,
    "Si": 1.90,
    "P": 2.19,
    "S": 2.58,
    "Cl": 3.16,
    "K": 0.82,
    "Ca": 1.00,
    "V": 1.63,
    "Cr": 1.66,
    "Mn": 1.55,
    "Fe": 1.83,
    "Co": 1.88,
    "Ni": 1.91,
    "Cu": 1.90,
    "Zn": 1.65,
    "As": 2.18,
    "Se": 2.55,
    "Br": 2.96,
    "Mo": 2.16,
    "I": 2.66,
    "W": 2.36,
    "Hg": 2.00,
    "U": 1.38,
}

#: Conventional oxidation states for spectator elements, applied in this order.
#: Order matters: the first matching rule wins, so halogens are resolved before
#: oxygen in species such as perchlorate.
_SPECTATOR_STATES: dict[str, int] = {
    "F": -1,
    "Li": 1,
    "Na": 1,
    "K": 1,
    "Rb": 1,
    "Cs": 1,
    "Be": 2,
    "Mg": 2,
    "Ca": 2,
    "Sr": 2,
    "Ba": 2,
    "H": 1,
    "O": -2,
    "Cl": -1,
    "Br": -1,
    "I": -1,
}

#: Species where the conventional rules give the wrong answer and the true
#: state is fixed by definition. Keyed by formula body (charge stripped).
_EXCEPTIONS: dict[str, dict[str, Fraction]] = {
    "H2O2": {"O": Fraction(-1), "H": Fraction(1)},
    "O2": {"O": Fraction(0)},
    "H2": {"H": Fraction(0)},
    "N2": {"N": Fraction(0)},
    "S2O3": {"S": Fraction(2)},  # thiosulfate: mean of S(-1) and S(+5)
}


def mean_oxidation_state(element: str, formula: str) -> Fraction:
    """Mean formal oxidation state of ``element`` in ``formula``.

    Returns an exact :class:`~fractions.Fraction` so that species such as
    thiosulfate or magnetite report ``+2`` and ``+8/3`` rather than a rounded
    float.
    """
    parsed = parse_formula(formula)
    if element not in parsed.elements:
        raise MicrobialThermoError(f"{formula!r} contains no {element}")

    body = _formula_body(formula)
    if body in _EXCEPTIONS and element in _EXCEPTIONS[body]:
        return _EXCEPTIONS[body][element]

    count = parsed.elements[element]

    # An elemental species carries its charge on the only element present.
    if len(parsed.elements) == 1:
        return Fraction(parsed.charge, count)

    known = Fraction(0)
    unresolved: list[str] = []
    for symbol, n in parsed.elements.items():
        if symbol == element:
            continue
        if symbol in _SPECTATOR_STATES:
            known += Fraction(_SPECTATOR_STATES[symbol] * n)
        else:
            unresolved.append(symbol)

    if unresolved:
        raise MicrobialThermoError(
            f"cannot assign a mean oxidation state for {element} in {formula!r}: "
            f"no conventional state for {', '.join(unresolved)}. "
            "Use atom_oxidation_states() with a SMILES instead."
        )

    return (Fraction(parsed.charge) - known) / count


def _formula_body(formula: str) -> str:
    from .formula import strip_charge

    body = strip_charge(formula)
    # Drop phase/state annotations such as "(aq)", "(g)", "(s)".
    for suffix in ("(aq)", "(g)", "(s)", "(l)"):
        body = body.replace(suffix, "")
    return body


def nosc(formula: str) -> Fraction:
    """Nominal Oxidation State of Carbon (LaRowe & Van Cappellen, 2011).

    For :math:`C_a H_b N_c O_d P_e S_f^{\\,Z}`:

    .. math::
        NOSC = 4 - \\frac{4a + b - 3c - 2d + 5e - 2f - Z}{a}

    Reference values used as test fixtures: CO2 -> +4, CH4 -> -4,
    acetate -> 0, glucose -> 0.
    """
    parsed = parse_formula(formula)
    a = parsed.element_count("C")
    if a == 0:
        raise MicrobialThermoError(f"NOSC is undefined for {formula!r}: no carbon")
    b = parsed.element_count("H")
    c = parsed.element_count("N")
    d = parsed.element_count("O")
    e = parsed.element_count("P")
    f = parsed.element_count("S")
    z = parsed.charge
    return Fraction(4) - Fraction(4 * a + b - 3 * c - 2 * d + 5 * e - 2 * f - z, a)


def atom_oxidation_states(smiles: str) -> list[tuple[str, int]]:
    """Per-atom formal oxidation states from a SMILES, via electronegativity.

    For each atom the oxidation number is its formal charge plus, for every
    bond, the bond order signed by which partner is more electronegative
    (homonuclear bonds contribute zero).

    Returns ``[(element_symbol, oxidation_number), ...]`` in RDKit atom order,
    explicit hydrogens included.

    Raises if the states do not sum to the molecular charge, which is the
    defining invariant of an oxidation-state assignment.
    """
    try:
        from rdkit import Chem
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise MicrobialThermoError(
            "atom_oxidation_states() requires rdkit; install the 'figures' extra"
        ) from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise MicrobialThermoError(f"RDKit could not parse SMILES {smiles!r}")
    mol = Chem.AddHs(mol)

    states: list[tuple[str, int]] = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        en_self = _electronegativity(symbol)
        oxidation = atom.GetFormalCharge()
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            en_other = _electronegativity(other.GetSymbol())
            order = int(round(bond.GetBondTypeAsDouble()))
            if en_other > en_self:
                oxidation += order
            elif en_other < en_self:
                oxidation -= order
        states.append((symbol, oxidation))

    total = sum(state for _, state in states)
    expected = Chem.GetFormalCharge(mol)
    if total != expected:
        raise MicrobialThermoError(
            f"oxidation states for {smiles!r} sum to {total}, expected {expected}; "
            "this indicates a bug in the assignment, not a chemistry edge case"
        )
    return states


def _electronegativity(symbol: str) -> float:
    try:
        return PAULING_ELECTRONEGATIVITY[symbol]
    except KeyError as exc:
        raise MicrobialThermoError(
            f"no Pauling electronegativity tabulated for {symbol!r}; "
            "add it to PAULING_ELECTRONEGATIVITY"
        ) from exc


def format_oxidation_state(state: Fraction) -> str:
    """Render an oxidation state for display: +6, -2, 0, +8/3."""
    if state == 0:
        return "0"
    sign = "+" if state > 0 else "-"
    magnitude = abs(state)
    if magnitude.denominator == 1:
        return f"{sign}{magnitude.numerator}"
    return f"{sign}{magnitude.numerator}/{magnitude.denominator}"
