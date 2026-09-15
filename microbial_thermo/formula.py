"""Chemical formula parsing: string -> element counts + charge.

Built on :mod:`chempy`, which parses charge as a first-class part of the
composition (stored under composition key ``0``) and uses exact integer counts.

``periodictable`` was evaluated and rejected here: it cannot parse charge
notation at all, which makes it unusable for redox work.

Accepted charge notations::

    SO4-2   SO4--   Fe+2   Fe++   NH4+   e-   H+

The doubled-sign forms are the legacy SUPCRT/GWB spellings that pyGCC's own
databases use, so both are supported and normalised to the same composition.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from functools import lru_cache

from chempy.util import periodic

from .exceptions import MicrobialThermoError

#: The free electron. Not a database species anywhere, but a first-class
#: participant in half reactions, so it gets an explicit composition.
ELECTRON = "e-"

_REPEATED_SIGN = re.compile(r"^(?P<body>.*?)(?P<signs>[+-]{2,})$")
_TRAILING_CHARGE = re.compile(r"^(?P<body>.*?)(?P<sign>[+-])(?P<count>\d*)$")


def _symbol_for_atomic_number(z: int) -> str:
    return periodic.symbols[z - 1]


def normalise_charge_notation(formula: str) -> str:
    """Convert legacy doubled-sign notation to the ``X-2`` form chempy expects.

    ``SO4--`` -> ``SO4-2``, ``Fe+++`` -> ``Fe+3``. Single signs and explicit
    counts are returned unchanged.
    """
    formula = formula.strip()
    m = _REPEATED_SIGN.match(formula)
    if m:
        signs = m.group("signs")
        # Mixed signs are meaningless; leave alone so chempy raises a clear error.
        if len(set(signs)) == 1:
            return f"{m.group('body')}{signs[0]}{len(signs)}"
    return formula


@dataclass(frozen=True)
class ParsedFormula:
    """Element counts and net charge for one chemical formula."""

    formula: str
    elements: dict[str, int] = field(default_factory=dict)
    charge: int = 0

    @property
    def is_electron(self) -> bool:
        return not self.elements and self.charge == -1

    def element_count(self, symbol: str) -> int:
        return self.elements.get(symbol, 0)

    def contains(self, symbol: str) -> bool:
        return self.elements.get(symbol, 0) > 0

    def __str__(self) -> str:
        return self.formula


@lru_cache(maxsize=2048)
def parse_formula(formula: str) -> ParsedFormula:
    """Parse ``formula`` into element counts and net charge.

    Results are cached; formulas are parsed constantly during balancing.
    """
    original = formula.strip()
    if not original:
        raise MicrobialThermoError("empty chemical formula")

    if original in ("e-", "e−"):
        return ParsedFormula(formula="e-", elements={}, charge=-1)

    candidate = normalise_charge_notation(original)

    from chempy import Substance  # imported lazily; chempy is slow to import

    try:
        with warnings.catch_warnings():
            # chempy still calls pyparsing's old camelCase API.
            warnings.simplefilter("ignore", DeprecationWarning)
            warnings.filterwarnings("ignore", message=".*deprecated.*")
            substance = Substance.from_formula(candidate)
    except Exception as exc:  # chempy raises a variety of parse errors
        raise MicrobialThermoError(f"could not parse chemical formula {original!r}: {exc}") from exc

    composition = dict(substance.composition)
    charge = int(composition.pop(0, 0))
    elements = {_symbol_for_atomic_number(z): int(n) for z, n in sorted(composition.items())}
    return ParsedFormula(formula=original, elements=elements, charge=charge)


def formula_charge(formula: str) -> int:
    return parse_formula(formula).charge


def element_counts(formula: str) -> dict[str, int]:
    return dict(parse_formula(formula).elements)


def all_elements(formulas) -> list[str]:
    """Sorted union of every element appearing across ``formulas``."""
    seen: set[str] = set()
    for f in formulas:
        seen.update(parse_formula(f).elements)
    return sorted(seen)


def strip_charge(formula: str) -> str:
    """Return the formula body without its trailing charge annotation."""
    candidate = normalise_charge_notation(formula)
    m = _TRAILING_CHARGE.match(candidate)
    if m and m.group("body"):
        return m.group("body")
    return candidate


def format_charge(charge: int) -> str:
    """Render a charge as a superscript-ready string: 2 -> '2+', -1 -> '-'."""
    if charge == 0:
        return ""
    sign = "+" if charge > 0 else "-"
    magnitude = abs(charge)
    return sign if magnitude == 1 else f"{magnitude}{sign}"
