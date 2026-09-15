"""The stacked half-reaction diagram.

Two equations, one above the other, aligned at their arrows: the electron
donor written in the oxidative direction on top, the electron acceptor written
in the reductive direction below. A curved arrow between them carries the
electron count. Oxidation states of the redox-active element sit above the top
equation and below the bottom one, centred on the species they describe.
E-standard and E-standard-prime are reported at the right.

Text is always measured with the Agg renderer, even when the output is SVG:
the SVG backend's text metrics can differ, and mixing the two would put the
oxidation-state labels off-centre in one format but not the other.
"""

from __future__ import annotations

from pathlib import Path

from ..typeset import ARROW, SPECIES, align_at_arrows, build_equation_tokens
from .style import PALETTE, SIZES

#: Vertical positions in axes coordinates.
_Y_TOP = 0.66
_Y_BOTTOM = 0.30

#: Horizontal space reserved for the role labels and the potential columns.
_LEFT_MARGIN = 0.14
_RIGHT_MARGIN = 0.20


#: Gap between adjacent tokens, in axes-fraction units.
_TOKEN_GAP = 0.008


def _measure(fig, ax, tokens, fontsize):
    """Fill in each token's width, in axes-fraction units.

    Measured against the Agg renderer regardless of the eventual output
    format, so SVG and PNG place identically.
    """
    renderer = fig.canvas.get_renderer()
    axes_width = ax.get_window_extent(renderer=renderer).width
    for token in tokens:
        artist = ax.text(0, -5, token.text, fontsize=fontsize)
        extent = artist.get_window_extent(renderer=renderer)
        token.width = extent.width / axes_width
        artist.remove()


def plot_half_reactions(
    reaction,
    n_electrons: int | None = 2,
    figsize=(11.0, 4.0),
    title: str | None = None,
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
):
    """Draw the two half reactions of ``reaction``, aligned at their arrows.

    Normalised to an electron pair by default, which is the convention the
    figure is specified around: both halves are balanced so exactly two
    electrons transfer. Pass ``n_electrons=None`` to draw the reaction exactly
    as it was built, or another integer to normalise to that instead.

    ``save`` gives a path without an extension; one file is written per entry
    in ``formats``. Returns the matplotlib figure.
    """
    import matplotlib.pyplot as plt

    if n_electrons is not None and reaction.n_electrons != n_electrons:
        reaction = reaction.renormalized(n_electrons)

    donor = reaction.donor_half
    acceptor = reaction.acceptor_half

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.canvas.draw()  # provision a renderer before measuring anything

    top = build_equation_tokens(
        donor.half, direction="oxidation", annotate_element=donor.half.key_element
    )
    bottom = build_equation_tokens(
        acceptor.half, direction="reduction", annotate_element=acceptor.half.key_element
    )

    for layout in (top, bottom):
        _measure(fig, ax, layout.tokens, SIZES["equation"])

    # Align on a common arrow position, then fit both equations into the space
    # between the role labels and the potentials column.
    #
    # Measure the *actual* extent after alignment rather than summing token
    # widths: the sum omits the inter-token gaps, which underestimates the true
    # width and lets a long equation overrun the potentials.
    align_at_arrows([top, bottom], gap=_TOKEN_GAP)
    tokens = [t for layout in (top, bottom) for t in layout.tokens]
    leftmost = min(t.x for t in tokens)
    extent = max(t.x + t.width for t in tokens) - leftmost

    available = 1.0 - _LEFT_MARGIN - _RIGHT_MARGIN
    scale = min(1.0, available / extent) if extent > 0 else 1.0

    for token in tokens:
        token.x = _LEFT_MARGIN + (token.x - leftmost) * scale
        token.width *= scale

    fontsize = SIZES["equation"] * scale

    _draw_equation(ax, top, _Y_TOP, PALETTE["oxidation"], "above", fontsize)
    _draw_equation(ax, bottom, _Y_BOTTOM, PALETTE["reduction"], "below", fontsize)

    _draw_role_labels(ax)
    _draw_electron_arrow(ax, top, bottom, reaction)
    _draw_potentials(ax, donor, acceptor)

    if title is None:
        title = _default_title(reaction)
    ax.set_title(title, fontsize=SIZES["title"], pad=14)

    fig.tight_layout()
    if save is not None:
        _save(fig, save, formats, dpi)
    return fig


def _draw_equation(ax, layout, y, colour, state_side, fontsize):
    """Draw one equation's tokens and its oxidation-state annotations."""
    offset = 0.11 if state_side == "above" else -0.11
    for token in layout.tokens:
        weight = "bold" if token.kind == ARROW else "normal"
        ax.text(
            token.x,
            y,
            token.text,
            fontsize=fontsize,
            color=colour if token.kind == ARROW else PALETTE["annotation"],
            ha="left",
            va="center",
            fontweight=weight,
        )
        if token.kind == SPECIES and token.oxidation_state:
            ax.text(
                token.center,
                y + offset,
                token.oxidation_state,
                fontsize=SIZES["oxidation_state"],
                color=colour,
                ha="center",
                va="center",
                fontweight="bold",
            )


def _draw_role_labels(ax):
    ax.text(
        0.005,
        _Y_TOP,
        "Oxidation",
        fontsize=SIZES["role_label"],
        color=PALETTE["oxidation"],
        ha="left",
        va="center",
        fontweight="bold",
    )
    ax.text(
        0.005,
        _Y_TOP - 0.075,
        "(electron donor)",
        fontsize=SIZES["annotation"],
        color=PALETTE["muted"],
        ha="left",
        va="center",
    )
    ax.text(
        0.005,
        _Y_BOTTOM,
        "Reduction",
        fontsize=SIZES["role_label"],
        color=PALETTE["reduction"],
        ha="left",
        va="center",
        fontweight="bold",
    )
    ax.text(
        0.005,
        _Y_BOTTOM - 0.075,
        "(electron acceptor)",
        fontsize=SIZES["annotation"],
        color=PALETTE["muted"],
        ha="left",
        va="center",
    )


def _draw_electron_arrow(ax, top, bottom, reaction):
    """Curved arrow from the donor equation to the acceptor equation."""
    x = max(top.arrow.center, bottom.arrow.center)
    ax.annotate(
        "",
        xy=(x, _Y_BOTTOM + 0.07),
        xytext=(x, _Y_TOP - 0.07),
        arrowprops={
            "arrowstyle": "-|>",
            "color": PALETTE["electron"],
            "linewidth": 1.6,
            "connectionstyle": "arc3,rad=-0.35",
            "shrinkA": 0,
            "shrinkB": 0,
        },
    )
    n = reaction.n_electrons
    label = f"{_pretty(n)} e$^-$"
    ax.text(
        x + 0.045,
        (_Y_TOP + _Y_BOTTOM) / 2.0,
        label,
        fontsize=SIZES["oxidation_state"],
        color=PALETTE["electron"],
        ha="left",
        va="center",
        fontweight="bold",
    )


def _draw_potentials(ax, donor, acceptor):
    """The E-standard-prime and E-standard columns on the right."""
    x = 1.0 - _RIGHT_MARGIN + 0.035
    ax.text(
        x,
        _Y_TOP + 0.20,
        r"$E^{\circ\prime}$" + "\n" + r"$E^{\circ}$",
        fontsize=SIZES["annotation"],
        color=PALETTE["muted"],
        ha="left",
        va="center",
        linespacing=1.5,
    )
    for result, y, colour in (
        (donor, _Y_TOP, PALETTE["oxidation"]),
        (acceptor, _Y_BOTTOM, PALETTE["reduction"]),
    ):
        primed = result.E_standard_prime.to("V").magnitude
        standard = result.E_standard.to("V").magnitude
        ax.text(
            x,
            y,
            f"{_signed(primed)} V\n{_signed(standard)} V",
            fontsize=SIZES["potential"],
            color=colour,
            ha="left",
            va="center",
            linespacing=1.5,
        )


def _default_title(reaction) -> str:
    from ..showwork import MATHTEXT_FRAC, latex_equation

    delta_g = reaction.delta_G_standard_prime.magnitude
    delta_e = reaction.delta_E_standard_prime.to("V").magnitude
    ph = reaction.conditions.pH
    temperature = reaction.conditions.temperature_c
    return (
        f"${latex_equation(reaction.coefficients, frac=MATHTEXT_FRAC)}$\n"
        f"$\\Delta E^{{\\circ\\prime}}$ = {delta_e:+.3f} V, "
        f"$\\Delta G^{{\\circ\\prime}}$ = {delta_g:+.1f} kJ/mol "
        f"({_pretty(reaction.n_electrons)} e$^-$, pH {ph:g}, {temperature:g} °C)"
    )


def _signed(value: float, places: int = 3) -> str:
    """Format with an explicit sign, without producing "-0.000"."""
    if abs(value) < 0.5 * 10**-places:
        return f"{0.0:.{places}f}"
    return f"{value:+.{places}f}"


def _pretty(value) -> str:
    """Render a Fraction as an integer when it is one."""
    try:
        if value.denominator == 1:
            return str(value.numerator)
        return f"{value.numerator}/{value.denominator}"
    except AttributeError:
        return str(value)


def _save(fig, save, formats, dpi):
    base = Path(save)
    base.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for suffix in formats:
        path = base.with_suffix(f".{suffix}")
        fig.savefig(path, format=suffix, dpi=dpi, bbox_inches="tight")
        written.append(path)
    return written
