"""A renderer-agnostic layout model for chemical equations.

The half-reaction figure has to stack two equations, align them at their
arrows, and centre an oxidation-state label over a particular species. Doing
that with hand-tuned coordinates is fragile and untestable, so layout is built
here as a data structure first: a list of tokens, each with a width and an
anchor, computed from measured text extents.

A renderer then draws the tokens. Only the renderer knows about matplotlib,
which keeps the layout logic testable without drawing anything.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from .balance import ELECTRON, HalfReaction
from .oxidation import format_oxidation_state, mean_oxidation_state

#: Token kinds. The arrow is the alignment anchor for stacked equations.
COEFFICIENT = "coefficient"
SPECIES = "species"
OPERATOR = "operator"
ARROW = "arrow"


@dataclass
class Token:
    """One drawable piece of an equation."""

    text: str
    kind: str
    species: object = None
    oxidation_state: str | None = None
    width: float = 0.0
    x: float = 0.0

    @property
    def is_arrow(self) -> bool:
        return self.kind == ARROW

    @property
    def center(self) -> float:
        return self.x + self.width / 2.0


@dataclass
class EquationLayout:
    """A laid-out equation: tokens with widths and x positions."""

    tokens: list[Token] = field(default_factory=list)
    arrow_index: int = -1

    @property
    def arrow(self) -> Token | None:
        if self.arrow_index < 0:
            return None
        return self.tokens[self.arrow_index]

    @property
    def left_width(self) -> float:
        """Total width of everything before the arrow."""
        return sum(t.width for t in self.tokens[: self.arrow_index])

    @property
    def right_width(self) -> float:
        return sum(t.width for t in self.tokens[self.arrow_index + 1 :])

    @property
    def total_width(self) -> float:
        return sum(t.width for t in self.tokens)

    def species_tokens(self) -> list[Token]:
        return [t for t in self.tokens if t.kind == SPECIES]


def format_coefficient_latex(value: Fraction) -> str:
    """Render a stoichiometric coefficient, omitting an implicit 1."""
    if value == 1:
        return ""
    if value.denominator == 1:
        return str(value.numerator)
    return rf"$\frac{{{value.numerator}}}{{{value.denominator}}}$"


def species_latex(species) -> str:
    """Wrap a species label as mathtext."""
    label = getattr(species, "display", "") or getattr(species, "formula", str(species))
    return f"${label}$"


def build_equation_tokens(
    half: HalfReaction,
    direction: str,
    annotate_element: str | None = None,
) -> EquationLayout:
    """Turn a half reaction into an ordered token list.

    ``direction`` is ``"oxidation"`` or ``"reduction"``. When
    ``annotate_element`` is given, species containing that element carry the
    mean oxidation state of it, which the renderer draws above or below.
    """
    coefficients = (
        half.coefficients
        if direction == "reduction"
        else {k: -v for k, v in half.coefficients.items()}
    )
    left = [(s, -v) for s, v in coefficients.items() if v < 0]
    right = [(s, v) for s, v in coefficients.items() if v > 0]

    # Put the redox-active species first on each side, then the electron last,
    # so the equation reads the way it would be written by hand.
    def sort_key(item):
        species, _ = item
        if species == ELECTRON:
            return (2, species.backend)
        if annotate_element and annotate_element in species.parsed.elements:
            return (0, species.backend)
        return (1, species.backend)

    left.sort(key=sort_key)
    right.sort(key=sort_key)

    layout = EquationLayout()

    def emit_side(terms):
        for index, (species, value) in enumerate(terms):
            if index:
                layout.tokens.append(Token(text="+", kind=OPERATOR))
            coefficient = format_coefficient_latex(value)
            if coefficient:
                layout.tokens.append(Token(text=coefficient, kind=COEFFICIENT))
            state = None
            if (
                annotate_element
                and species != ELECTRON
                and annotate_element in species.parsed.elements
            ):
                state = format_oxidation_state(
                    mean_oxidation_state(annotate_element, species.formula)
                )
            layout.tokens.append(
                Token(
                    text=species_latex(species),
                    kind=SPECIES,
                    species=species,
                    oxidation_state=state,
                )
            )

    emit_side(left)
    layout.arrow_index = len(layout.tokens)
    layout.tokens.append(Token(text=r"$\longrightarrow$", kind=ARROW))
    emit_side(right)
    return layout


def assign_positions(layout: EquationLayout, arrow_x: float, gap: float = 0.0) -> EquationLayout:
    """Place tokens so the arrow's left edge sits at ``arrow_x``.

    Tokens before the arrow are laid out right-to-left from the arrow; tokens
    after it run left-to-right. This is what makes two stacked equations line
    up on their arrows regardless of how wide their left-hand sides are.
    """
    cursor = arrow_x
    for token in reversed(layout.tokens[: layout.arrow_index]):
        cursor -= token.width + gap
        token.x = cursor

    cursor = arrow_x
    for token in layout.tokens[layout.arrow_index :]:
        token.x = cursor
        cursor += token.width + gap
    return layout


def align_at_arrows(layouts, gap: float = 0.0) -> float:
    """Align several equations on a common arrow position.

    Returns the chosen arrow x, which is set by the widest left-hand side so
    that no equation overflows to the left.
    """
    arrow_x = max(layout.left_width + gap * max(layout.arrow_index, 0) for layout in layouts)
    for layout in layouts:
        assign_positions(layout, arrow_x, gap=gap)
    return arrow_x
