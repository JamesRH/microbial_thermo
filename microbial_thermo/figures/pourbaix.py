"""Eh-pH (Pourbaix) diagrams: which form of an element wins, and where.

A predominance diagram over the two variables that actually decide aqueous
speciation. Every species of an element is formed from a common basis --
the element itself, water, protons and electrons -- and at each point in
(Eh, pH) the one with the lowest free energy per mole of element is drawn.

The construction, since it is short enough to state:

    n_E * E + b * H2O + c * H+ + d * e-  ->  S

with ``b`` set by the oxygen in S, ``c`` by its hydrogen after the water is
accounted for, and ``d`` by charge. Then at a given pH and Eh

    dG(S) = dGf(S) + RT ln(a) - n_E*dGf(E) - b*dGf(H2O)
            + c * RT ln(10) * pH + d * F * Eh

all divided by ``n_E``, because the comparison is per mole of the element.
Protons and electrons carry no formation energy on the SHE convention, which
is what makes those last two terms the whole pH and Eh dependence.

**Why this is built rather than borrowed.** pymatgen has a mature
``PourbaixDiagram``, but its pipeline wants a Materials Project API key and
builds entries from DFT solid energies plus experimental ion energies. Mixing
those with our HKF aqueous species would put two provenances inside one
diagram -- the same objection that rules out PHREEQC for speciation. The
geometry here is a few dozen lines; the data consistency is the part worth
protecting.

The dissolved activity matters and is not a detail. A Pourbaix diagram is
drawn at a chosen total dissolved activity, conventionally 1e-6, and the
solid fields grow as you lower it. The value is always printed on the figure
for that reason.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..units import FARADAY, KJ_PER_MOL_STR, R, as_magnitude, celsius_to_kelvin, ureg
from .style import PALETTE, SIZES

#: Conventional dissolved activity for a Pourbaix diagram.
DEFAULT_ACTIVITY = 1e-6

#: Basis species and activities for elements other than the one being drawn.
#: An iron diagram has a siderite field only because there is carbonate about,
#: and how much is a modelling input rather than a detail -- so these are
#: stated, and printed on the figure when they are used.
#: The reference form of each element -- its zero point, in whatever species
#: the element actually exists as. Mostly the monatomic solid, but nitrogen is
#: N2 and carbon is graphite, and monatomic "N" is in no database at all. The
#: atom count is read from the formula, so a diatomic reference is handled by
#: taking a half of it.
ELEMENT_REFERENCE = {
    "N": "N2(aq)",
    "C": "Graphite",
    "S": "Sulfur(s)",
    "O": "O2(aq)",
    "H": "H2(aq)",
}

FIXED_DEFAULTS = {
    "S": ("SO4-2", 1e-3),
    "C": ("HCO3-", 2e-3),
    "P": ("HPO4-2", 1e-5),
}

#: Elements we can draw straight away, with the species to include. Kept
#: explicit rather than "everything containing the element": a diagram with
#: forty fields teaches nothing, and the choice of which phases to allow is a
#: modelling decision the caller should see.
DEFAULT_SPECIES = {
    "As": ["As", "As(OH)3(aq)", "H2AsO3-", "H3AsO4(aq)", "H2AsO4-", "HAsO4--", "AsO4---"],
    "Fe": ["Fe", "Fe+2", "Fe+3", "Goethite", "Hematite", "Magnetite", "Siderite", "Pyrite"],
    "Mn": ["Mn", "Mn+2", "Mn+3", "Pyrolusite", "Manganite", "Hausmannite", "Bixbyite"],
    "S": ["Sulfur(s)", "H2S(aq)", "HS-", "SO4--", "HSO4-", "SO3--", "HSO3-"],
    "Se": ["Se", "SeO4--", "HSeO4-", "SeO3--", "HSeO3-"],
    "N": ["N2(aq)", "NO3-", "NO2-", "NH4+", "NH3(aq)"],
    "C": ["Graphite", "CO2(aq)", "HCO3-", "CO3--", "Methane(aq)", "Acetate", "Formate(aq)"],
    "Cr": ["Cr++", "Cr+++", "CrO4--", "HCrO4-"],
    "U": ["U++++", "UO2++", "Uraninite"],
    "Cu": ["Cu", "Cu+", "Cu++", "Covellite"],
}


@dataclass
class PourbaixField:
    """Which species predominates across an (Eh, pH) grid."""

    element: str
    ph: np.ndarray
    eh: np.ndarray
    winner: np.ndarray  # index into `species`, shape (len(eh), len(ph))
    species: tuple
    activity: float
    temperature_c: float
    fixed: tuple = ()

    @property
    def present(self) -> list:
        """Species that actually win somewhere. Often far fewer than offered."""
        return [self.species[i] for i in sorted(set(self.winner.ravel().tolist()))]

    def at(self, pH: float, eh: float) -> str:
        i = int(np.argmin(np.abs(self.eh - eh)))
        j = int(np.argmin(np.abs(self.ph - pH)))
        return self.species[int(self.winner[i, j])]


def _basis_coefficients(species, element: str, auxiliaries=()):
    """Decompose ``species`` onto the basis.

    Returns ``(n_E, b, c, d, extras)`` for

        n_E E + b H2O + c H+ + d e- + sum(k_X B_X)  ->  species

    where each ``B_X`` is an auxiliary basis species carrying one extra
    element at a fixed activity -- sulfate for sulfur, bicarbonate for carbon.

    Without those, a species like siderite or pyrite has elements nothing in
    the basis accounts for, and its energy comes out meaningless (and, as it
    happens, spuriously low, so it wins everywhere). That was the first
    version of this function and it drew an iron diagram made entirely of
    siderite and pyrite.
    """
    parsed = species.parsed
    n_e = parsed.element_count(element)
    if n_e == 0:
        raise ValueError(f"{species.backend} contains no {element}")

    extras = []
    for aux_element, aux_species in auxiliaries:
        needed = parsed.element_count(aux_element)
        if not needed:
            continue
        per = aux_species.parsed.element_count(aux_element)
        extras.append((aux_element, aux_species, needed / per))

    unexplained = set(parsed.elements) - {element, "H", "O"} - {name for name, _, _ in extras}
    if unexplained:
        raise ValueError(
            f"{species.backend} contains {', '.join(sorted(unexplained))}, which "
            "the basis does not account for. Give an auxiliary species for each "
            "via fixed={'S': 'SO4-2', ...}, or leave the species out."
        )

    b = parsed.element_count("O") - sum(k * aux.parsed.element_count("O") for _, aux, k in extras)
    c = (
        parsed.element_count("H")
        - 2 * b
        - sum(k * aux.parsed.element_count("H") for _, aux, k in extras)
    )
    d = c - parsed.charge + sum(k * aux.parsed.charge for _, aux, k in extras)
    return n_e, b, c, d, extras


def _reference_for(element: str, entries, gibbs):
    """A zero point for the element: species, atoms of the element, energy.

    **The choice does not affect the diagram.** Every species picks up the
    same ``-G(reference)/atoms`` per mole of element, so the term is a
    constant offset and cancels in the comparison. It is kept because it makes
    the energies mean something -- formation from the element -- and dropped
    to a fallback without ceremony when no elemental form exists.

    Which is often. Monatomic nitrogen is in no database, uranium and chromium
    metal are in none of ours, and carbon's reference is graphite. So the
    fallback is simply the first species offered.
    """
    from ..species import resolve

    candidates = [ELEMENT_REFERENCE.get(element, element)]
    candidates += [entry.backend for entry in entries]
    for name in candidates:
        try:
            species = resolve(name)
            atoms = species.parsed.element_count(element)
            if atoms:
                return species, atoms, gibbs(species.backend)
        except Exception:
            continue
    raise ValueError(f"no usable reference form for {element!r}")


def pourbaix_field(
    element: str,
    species=None,
    temperature_c: float = 25.0,
    activity: float = DEFAULT_ACTIVITY,
    ph_range=(0.0, 14.0),
    eh_range=(-1.0, 1.4),
    points: int = 240,
    fixed=None,
    backend=None,
) -> PourbaixField:
    """Work out the predominant species of ``element`` across (Eh, pH).

    ``fixed`` supplies a basis species and activity for every element other
    than the one being drawn, as ``{"S": ("SO4-2", 1e-3)}``. Iron's diagram
    only has a siderite or pyrite field because carbonate and sulfide are
    present at some concentration, so that concentration is an input, not a
    detail -- and a species whose extra elements are not covered is refused
    rather than silently mis-weighted.
    """
    from .. import get_backend
    from ..species import resolve

    backend = backend or get_backend()
    names = species or DEFAULT_SPECIES.get(element)
    if not names:
        raise ValueError(f"no default species list for {element!r}; pass species= explicitly")
    entries = [resolve(n) for n in names]
    fixed = dict(FIXED_DEFAULTS if fixed is None else fixed)
    auxiliaries = [
        (el, resolve(spec if isinstance(spec, str) else spec[0]))
        for el, spec in fixed.items()
        if el != element
    ]
    aux_activity = {
        el: (DEFAULT_ACTIVITY if isinstance(spec, str) else spec[1]) for el, spec in fixed.items()
    }

    kelvin = celsius_to_kelvin(temperature_c)
    rt = as_magnitude((R * (kelvin * ureg.kelvin)).to(KJ_PER_MOL_STR), KJ_PER_MOL_STR)
    rt_ln10 = rt * np.log(10.0)
    faraday = FARADAY.to("C/mol").magnitude / 1000.0  # kJ/(mol V)

    def gibbs(name):
        return backend.delta_Gf(name, temperature_c).to(KJ_PER_MOL_STR).magnitude

    water = gibbs("H2O")
    reference, reference_atoms, reference_g = _reference_for(element, entries, gibbs)

    ph = np.linspace(*ph_range, points)
    eh = np.linspace(*eh_range, points)
    grid_ph, grid_eh = np.meshgrid(ph, eh)

    energies, kept = [], []
    for entry in entries:
        try:
            n_e, b, c, d, extras = _basis_coefficients(entry, element, auxiliaries)
        except ValueError:
            continue
        base = gibbs(entry.backend) - (n_e / reference_atoms) * reference_g - b * water
        for aux_element, aux, k in extras:
            base -= k * (
                gibbs(aux.backend)
                + (
                    rt * np.log(aux_activity.get(aux_element, DEFAULT_ACTIVITY))
                    if aux.is_aqueous
                    else 0.0
                )
            )
        if entry.is_aqueous:
            base += rt * np.log(activity)
        energies.append((base + c * rt_ln10 * grid_ph + d * faraday * grid_eh) / n_e)
        kept.append(entry)
    if not energies:
        raise ValueError(f"no usable species for {element!r}")
    entries = kept

    stacked = np.stack(energies)
    return PourbaixField(
        element=element,
        ph=ph,
        eh=eh,
        winner=np.argmin(stacked, axis=0),
        species=tuple(e.backend for e in entries),
        activity=activity,
        temperature_c=temperature_c,
        fixed=tuple(
            (el, aux.backend, aux_activity.get(el, DEFAULT_ACTIVITY))
            for el, aux in auxiliaries
            if any(el in e.parsed.elements for e in entries)
        ),
    )


def water_stability(temperature_c: float = 25.0, ph=None, backend=None):
    """The Eh bounds within which liquid water is stable.

    Computed from this library's own O2/H2O and H+/H2 couples rather than
    quoted as 1.229 and 0, so the lines move correctly with temperature.
    Returns ``(pH, upper, lower)``.
    """
    from .. import get_backend
    from ..conditions import Conditions
    from ..reaction import Couple, HalfReactionResult

    backend = backend or get_backend()
    ph = np.linspace(0.0, 14.0, 57) if ph is None else np.asarray(ph)

    def line(reduced, oxidized):
        values = []
        for value in ph:
            conditions = Conditions(temperature_c=temperature_c, pH=float(value))
            couple = Couple.make(reduced, oxidized)
            result = HalfReactionResult(
                couple=couple,
                half=couple.half_reaction(),
                conditions=conditions,
                backend=backend,
            )
            values.append(result.E_standard_prime.to("V").magnitude)
        return np.asarray(values)

    # O2(g) and H2(g) at 1 bar, which is what the published 1.229 V refers to.
    # Against O2(aq) at unit activity the upper line sits 43 mV higher -- a
    # different reference state, not an error, but the wrong one to draw here.
    return ph, line("H2O", "O2(g)"), line("H2(g)", "H+")


def plot_pourbaix(
    element: str,
    field: PourbaixField | None = None,
    show_water: bool = True,
    figsize=(8.5, 6.5),
    save: str | Path | None = None,
    formats=("svg", "png"),
    dpi: int = 200,
    **kwargs,
):
    """Draw the predominance diagram. Returns ``(figure, field)``."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    if field is None:
        field = pourbaix_field(element, **kwargs)

    present = sorted(set(field.winner.ravel().tolist()))
    remap = {old: new for new, old in enumerate(present)}
    drawn = np.vectorize(remap.get)(field.winner)
    labels = [field.species[i] for i in present]

    palette = plt.cm.tab20(np.linspace(0, 1, max(len(labels), 2)))
    figure, ax = plt.subplots(figsize=figsize)
    ax.pcolormesh(
        field.ph,
        field.eh,
        drawn,
        cmap=ListedColormap(palette[: len(labels)]),
        shading="auto",
        alpha=0.75,
    )
    ax.contour(
        field.ph,
        field.eh,
        drawn,
        levels=np.arange(len(labels)) + 0.5,
        colors=PALETTE["annotation"],
        linewidths=0.8,
    )

    # Label each field at its centroid, which is where it is widest.
    for new, name in enumerate(labels):
        mask = drawn == new
        if mask.sum() < drawn.size * 0.004:
            continue
        rows, cols = np.nonzero(mask)
        ax.text(
            field.ph[int(np.median(cols))],
            field.eh[int(np.median(rows))],
            _pretty(name),
            ha="center",
            va="center",
            fontsize=SIZES["annotation"] + 1,
            color=PALETTE["annotation"],
            bbox={"facecolor": "white", "alpha": 0.55, "edgecolor": "none", "pad": 1.5},
        )

    if show_water:
        ph, upper, lower = water_stability(field.temperature_c)
        for values, label in ((upper, "O$_2$/H$_2$O"), (lower, "H$^+$/H$_2$")):
            ax.plot(ph, values, color=PALETTE["reduction"], linestyle="--", linewidth=1.3)
            ax.annotate(
                label,
                xy=(ph[-1], values[-1]),
                xytext=(-4, 4),
                textcoords="offset points",
                ha="right",
                fontsize=SIZES["annotation"] - 1,
                color=PALETTE["reduction"],
            )

    ax.set_xlabel("pH", fontsize=SIZES["potential"])
    ax.set_ylabel("Eh (V vs SHE)", fontsize=SIZES["potential"])
    ax.set_xlim(field.ph[0], field.ph[-1])
    ax.set_ylim(field.eh[0], field.eh[-1])
    subtitle = f"dissolved activity {field.activity:g} · {field.temperature_c:g} °C"
    if field.fixed:
        held = ", ".join(f"{_pretty(name)} {value:g}" for _, name, value in field.fixed)
        subtitle += f"\nheld fixed: {held}"
    ax.set_title(f"{element} predominance\n{subtitle}", fontsize=SIZES["title"], pad=10)
    figure.tight_layout()

    if save is not None:
        base = Path(save)
        if base.parent != Path(""):
            base.parent.mkdir(parents=True, exist_ok=True)
        for suffix in formats:
            figure.savefig(
                base.with_suffix(f".{suffix}"), format=suffix, dpi=dpi, bbox_inches="tight"
            )
    return figure, field


def _pretty(name: str) -> str:
    """Backend names are legacy SUPCRT; make them readable on a figure."""
    from ..species import resolve

    try:
        display = resolve(name).display
    except Exception:
        display = ""
    return f"${display}$" if display else name
