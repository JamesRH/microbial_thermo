"""The common basis behind Pourbaix, Latimer and Frost diagrams.

All three answer questions about one element's redox forms, and all three
need the same thing first: every species written as a formation from the
element itself, water, protons and electrons.

    n_E E + b H2O + c H+ + d e- (+ auxiliaries)  ->  S

Once that decomposition exists, everything else is arithmetic on it. Per mole
of the element,

    G(S; pH, Eh) = base + (c/n_E) RT ln10 pH + (d/n_E) F Eh

where ``base`` gathers the formation energies, the activity term and the
auxiliary basis species. The three diagrams differ only in what they do with
that expression:

* **Pourbaix** evaluates it on a (pH, Eh) grid and takes the lowest at each
  point.
* **Latimer** sets Eh to zero, compares pairs, and reports the potential at
  which two forms are in equilibrium.
* **Frost** sets Eh to zero and plots ``base / F`` -- the volt equivalent --
  against the oxidation state.

Two facts about the element reference that are easy to get backwards, and
which are the reason this module exists rather than three copies of it:

**For Pourbaix the reference cancels; for Frost it does not.** Every species
picks up the same ``-G(reference)/atoms`` per mole of element, so it shifts
every curve equally and moves no boundary. A Frost diagram plots the value
itself, so the reference *is* the zero of the y axis and must be the element
in its standard state. Drawing nitrogen against ``N2(aq)`` rather than
``N2(g)`` puts every point 0.094 V low -- enough to turn the textbook -0.82 V
for ammonium into -0.92 V, which is how this was found.

**The oxidation state falls out of ``d``, and needs no separate convention.**
``d`` is the electrons consumed in forming the species from the element, so
the oxidation state per atom is ``-d/n_E``. Magnetite gives -(-8)/3 = +8/3,
ammonium -3, nitrate +5. Nothing here treats negative states specially, which
is what makes the volt equivalent of a hydride come out right.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..exceptions import OutOfRangeError
from ..units import FARADAY, KJ_PER_MOL_STR, R, as_magnitude, celsius_to_kelvin, ureg

#: Faraday constant in kJ/(mol V), so a potential times electrons is an energy.
FARADAY_KJ = FARADAY.to("C/mol").magnitude / 1000.0

#: Conventional dissolved activity for a Pourbaix diagram.
DEFAULT_ACTIVITY = 1e-6

#: The element's zero point, in whatever species it actually exists as.
#:
#: These are standard states: the gas at 1 bar for nitrogen, oxygen and
#: hydrogen, graphite for carbon, rhombic sulfur, and the metal otherwise.
#: Monatomic "N" is in no database at all, which is why nitrogen needs an
#: entry here; carbon's would otherwise resolve to nothing.
ELEMENT_REFERENCE = {
    "N": "N2(g)",
    "C": "Graphite",
    "S": "Sulfur(s)",
    "O": "O2(g)",
    "H": "H2(g)",
}

#: Basis species and activities for elements other than the one being drawn.
#: An iron diagram has a siderite field only because there is carbonate about,
#: and how much is a modelling input rather than a detail -- so these are
#: stated, and printed on the figure when they are used.
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


def basis_coefficients(species, element: str, auxiliaries=()):
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


def reference_for(element: str, entries, gibbs, preferred=None):
    """A zero point for the element: species, atoms of the element, energy.

    For a predominance diagram **the choice does not matter** -- every species
    picks up the same ``-G(reference)/atoms`` per mole of element, so the term
    is a constant offset that cancels in the comparison. For a Frost diagram
    it matters completely, because the value plotted *is* that difference.

    So the fallback below, which simply takes the first species offered, is
    harmless on an Eh-pH diagram and wrong on a Frost one. Callers who plot
    the energy itself should check ``ElementSeries.reference_is_element``.
    """
    from ..species import resolve

    candidates = [preferred] if preferred else []
    candidates += [ELEMENT_REFERENCE.get(element, element)]
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


@dataclass(frozen=True)
class SpeciesEnergy:
    """One species of an element, written per mole of that element.

    ``base`` is everything independent of pH and Eh; the two coefficients
    carry the rest. Splitting it this way is what lets one decomposition
    serve a scalar Latimer potential and a 240x240 Pourbaix grid.
    """

    species: object
    n_element: int
    base: float  # kJ per mole of element, at pH 0 and Eh 0
    protons: float  # c / n_E
    electrons: float  # d / n_E
    auxiliaries: tuple = ()  # (element, backend) held at fixed activity

    @property
    def backend(self) -> str:
        return self.species.backend

    @property
    def oxidation_state(self) -> float:
        """The element's oxidation state, straight from the electron count.

        Negative states need no special handling: ``d`` counts electrons
        *consumed* in forming the species from the element, so a hydride has
        a positive ``d`` and a negative state without any extra convention.
        """
        # 0.0 rather than -0.0: the sign flip makes an elemental form print as
        # "-0" everywhere it is formatted, which looks like a bug.
        return 0.0 if self.electrons == 0 else -self.electrons

    @property
    def state_is_conventional(self) -> bool:
        """Whether :attr:`oxidation_state` is the element's real state.

        The electron count is measured against the basis, so it is an
        oxidation state only when every *other* element in the species sits
        at the same state in the species as in the basis species holding it.
        Hydrogen (+1) and oxygen (-2) always do. An auxiliary element need
        not, and pyrite is the case that matters: its sulfur is S(-I) while
        the sulfur basis is sulfate, so the count comes out at -12 for an
        iron that is plainly Fe(II).

        So this refuses any species carrying an auxiliary element at all --
        including siderite, whose carbonate carbon *is* C(+IV) in both and
        whose +2 would have been right. Telling those two apart needs an
        oxidation-state assignment independent of this decomposition, and
        the one this library has refuses metal sulfides outright, which is
        exactly the case in question. A mechanical rule that is never wrong
        beats a clever one that is sometimes wrong and silent about it.

        Nothing is lost by it. Siderite is Fe(II), so it can only ever be an
        alternative form on a rung Fe(2+) already occupies, and its energy
        depends on a carbonate activity that a single-element ladder has
        nowhere to state. The Eh-pH diagram keeps it, and that is where a
        siderite field belongs.
        """
        return not self.auxiliaries


@dataclass(frozen=True)
class ElementSeries:
    """Every usable redox form of one element, decomposed onto the basis."""

    element: str
    entries: tuple
    reference: str
    reference_atoms: int
    reference_gibbs: float
    temperature_c: float
    activity: float
    pH: float
    rt: float
    fixed: tuple = ()
    skipped: tuple = ()

    @property
    def rt_ln10(self) -> float:
        return self.rt * np.log(10.0)

    @property
    def reference_is_element(self) -> bool:
        """Whether the zero point is a genuine elemental form.

        False when the element has no elemental species available and the
        fallback picked an arbitrary one. A Frost diagram drawn on a series
        where this is False has a meaningless y-axis zero, and says so.

        The test is neutral *and* single-element, both parts earned: neither
        chromium nor uranium metal is in any of our databases, so the
        fallback lands on Cr(2+) and U(4+), whose only element is the one in
        question and which are emphatically not zero points.
        """
        from ..species import resolve

        if self.reference_gibbs is None:
            return False
        try:
            parsed = resolve(self.reference).parsed
        except Exception:
            return False
        return set(parsed.elements) == {self.element} and parsed.charge == 0

    def energy(self, entry: SpeciesEnergy, pH=None, eh=0.0):
        """Free energy of ``entry`` per mole of the element, kJ/mol.

        ``pH`` and ``eh`` may be arrays, which is how the Pourbaix grid is
        built; they broadcast together.
        """
        ph_value = self.pH if pH is None else pH
        return (
            entry.base
            + entry.protons * self.rt_ln10 * np.asarray(ph_value)
            + entry.electrons * FARADAY_KJ * np.asarray(eh)
        )

    def energies(self, pH=None, eh=0.0):
        """:meth:`energy` for every entry, aligned with :attr:`entries`."""
        return np.array([self.energy(entry, pH, eh) for entry in self.entries])

    def rescaled(self, activity=None, pH=None) -> ElementSeries:
        """The same series at a different dissolved activity or pH.

        Both are closed form -- an activity change is ``RT ln(a2/a1)`` per mole
        of the element on the aqueous entries only, and pH is already a
        coefficient -- so neither touches the backend. That is what makes a
        slider possible: the expensive axis is temperature, and only
        temperature.
        """
        import dataclasses

        entries = self.entries
        if activity is not None and activity != self.activity:
            shift = self.rt * np.log(activity / self.activity)
            entries = tuple(
                dataclasses.replace(entry, base=entry.base + shift / entry.n_element)
                if entry.species.is_aqueous
                else entry
                for entry in entries
            )
        return dataclasses.replace(
            self,
            entries=entries,
            activity=self.activity if activity is None else activity,
            pH=self.pH if pH is None else pH,
        )

    def __iter__(self):
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)


def element_series(
    element: str,
    species=None,
    temperature_c: float = 25.0,
    activity: float = 1.0,
    pH: float = 0.0,
    fixed=None,
    backend=None,
    reference=None,
    skip_out_of_range: bool = False,
) -> ElementSeries:
    """Decompose every listed form of ``element`` onto the common basis.

    ``activity`` is the dissolved activity assumed for every aqueous species.
    Unit activity is the standard-state convention a Latimer or Frost diagram
    is quoted at; a Pourbaix diagram conventionally uses 1e-6.

    ``fixed`` supplies a basis species and activity for every element other
    than the one being drawn, as ``{"S": ("SO4-2", 1e-3)}``. A species whose
    extra elements are not covered is recorded in ``skipped`` with the reason
    rather than silently mis-weighted.

    **An auxiliary activity of exactly zero means the element is absent**, and
    every species needing it is dropped. Arithmetically ``log(0)`` gives an
    energy of ``+inf`` and the same species lose everywhere, so the picture
    was already right; but an infinity in the array is a bug waiting for the
    first caller who averages or interpolates one, and "there is no sulfide
    here" is a statement the code should make in words. A negative activity is
    refused outright, and the element being drawn may not be set to zero at
    all -- that would mean there is none of it to draw.

    ``reference`` overrides the element's zero point. It changes nothing on
    an Eh-pH diagram and everything on a Frost one, which is the whole reason
    it is exposed.

    ``skip_out_of_range`` turns a species the backend cannot evaluate at this
    temperature into a ``skipped`` entry instead of an error. It is off by
    default, because losing a phase silently changes the diagram and the
    backend's own message names the species and says what to do. It is on
    for :func:`element_grid`, which spans temperatures and has to decide
    something.
    """
    from .. import get_backend
    from ..species import resolve

    backend = backend or get_backend()
    names = species or DEFAULT_SPECIES.get(element)
    if not names:
        raise ValueError(f"no default species list for {element!r}; pass species= explicitly")

    entries = [resolve(n) if isinstance(n, str) else n for n in names]
    fixed = dict(FIXED_DEFAULTS if fixed is None else fixed)
    auxiliaries = [
        (el, resolve(spec if isinstance(spec, str) else spec[0]))
        for el, spec in fixed.items()
        if el != element
    ]
    aux_activity = {
        el: (DEFAULT_ACTIVITY if isinstance(spec, str) else spec[1]) for el, spec in fixed.items()
    }

    if activity <= 0:
        raise ValueError(
            f"activity must be positive, not {activity!r}; it is the dissolved "
            f"activity of {element} itself, and zero would mean there is none to draw"
        )
    for name, value in aux_activity.items():
        if value < 0:
            raise ValueError(
                f"the fixed activity for {name} is negative ({value!r}); use 0 to "
                "mean the element is absent"
            )

    #: Auxiliary elements declared absent. A species that needs one cannot
    #: form, so it is dropped rather than evaluated at log(0).
    absent = {name for name, value in aux_activity.items() if value == 0}

    kelvin = celsius_to_kelvin(temperature_c)
    rt = as_magnitude((R * (kelvin * ureg.kelvin)).to(KJ_PER_MOL_STR), KJ_PER_MOL_STR)

    def gibbs(name):
        return backend.delta_Gf(name, temperature_c).to(KJ_PER_MOL_STR).magnitude

    water = gibbs("H2O")
    reference, reference_atoms, reference_g = reference_for(element, entries, gibbs, reference)

    built, skipped = [], []
    for entry in entries:
        try:
            n_e, b, c, d, extras = basis_coefficients(entry, element, auxiliaries)
        except ValueError as exc:
            skipped.append((entry.backend, str(exc)))
            continue
        missing = sorted({el for el, _, _ in extras} & absent)
        if missing:
            skipped.append(
                (
                    entry.backend,
                    f"needs {', '.join(missing)}, which is fixed at zero activity and so is absent",
                )
            )
            continue
        try:
            entry_gibbs = gibbs(entry.backend)
        except OutOfRangeError as exc:
            if not skip_out_of_range:
                raise
            skipped.append((entry.backend, str(exc)))
            continue
        base = entry_gibbs - (n_e / reference_atoms) * reference_g - b * water
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
        built.append(
            SpeciesEnergy(
                species=entry,
                n_element=n_e,
                base=base / n_e,
                protons=c / n_e,
                electrons=d / n_e,
                auxiliaries=tuple((el, aux.backend) for el, aux, _ in extras),
            )
        )

    if not built:
        raise ValueError(f"no usable species for {element!r}")

    return ElementSeries(
        element=element,
        entries=tuple(built),
        reference=reference.backend,
        reference_atoms=reference_atoms,
        reference_gibbs=reference_g,
        temperature_c=temperature_c,
        activity=activity,
        pH=pH,
        rt=rt,
        fixed=tuple(
            (el, aux.backend, aux_activity.get(el, DEFAULT_ACTIVITY))
            for el, aux in auxiliaries
            if any(el in e.species.parsed.elements for e in built)
        ),
        skipped=tuple(skipped),
    )


def pretty(name: str) -> str:
    """Backend names are legacy SUPCRT; make them readable on a figure."""
    from ..species import resolve

    try:
        display = resolve(name).display
    except Exception:
        display = ""
    return f"${display}$" if display else name


def ladder_entries(series):
    """Split a series into entries usable as oxidation-state rungs, and not.

    A Latimer or Frost diagram plots against the oxidation state, so it may
    only use species whose electron count *is* one --
    :attr:`SpeciesEnergy.state_is_conventional`. Returns
    ``(kept, [(backend, reason), ...])``.

    This is what keeps pyrite off the iron ladder. Before it existed, iron's
    Latimer diagram ran ``Fe(0) --+0.462 V--> pyrite`` with twelve electrons,
    which is a real reaction and a nonsensical rung: the electrons go to the
    sulfur, not the iron. Siderite goes with it, for the reason given on
    :attr:`SpeciesEnergy.state_is_conventional`.
    """
    kept, dropped = [], []
    for entry in series:
        if entry.state_is_conventional:
            kept.append(entry)
        else:
            held = ", ".join(f"{el} (as {name})" for el, name in entry.auxiliaries)
            dropped.append(
                (
                    entry.backend,
                    f"its electron count is measured against {held}, so it is not "
                    f"{series.element}'s oxidation state",
                )
            )
    return kept, dropped


#: Temperatures a precomputed element grid covers by default. Chosen for what
#: a course actually asks about: ice-covered water, a temperate lake, the
#: bench, a warm aquifer, a compost heap and a hot spring.
DEFAULT_TEMPERATURES = (2.0, 10.0, 25.0, 40.0, 60.0, 80.0)


@dataclass(frozen=True)
class ElementGrid:
    """One element's series precomputed across temperature.

    Temperature is the only axis that costs backend calls, so it is the only
    one gridded. pH and activity are applied in closed form on lookup, which
    keeps them continuous rather than quantised to slider steps.
    """

    element: str
    temperatures: tuple
    series: tuple
    activity: float
    dropped: tuple = ()  # (backend, why) -- not available at every temperature

    def at(self, temperature_c: float, pH: float = 7.0, activity=None) -> ElementSeries:
        index = min(
            range(len(self.temperatures)),
            key=lambda i: abs(self.temperatures[i] - temperature_c),
        )
        return self.series[index].rescaled(activity=activity, pH=pH)

    @property
    def species(self) -> tuple:
        return tuple(entry.backend for entry in self.series[0])


def element_grid(
    element: str,
    species=None,
    temperatures=DEFAULT_TEMPERATURES,
    activity: float = 1.0,
    fixed=None,
    backend=None,
    reference=None,
) -> ElementGrid:
    """Build an :class:`ElementGrid`, one series per temperature.

    Costs one pass over the backend per temperature.

    **A species has to be available at every temperature or at none**, because
    a slider that changes which species exist cannot be indexed and would
    silently redraw a different diagram. Manganese is the case that forced
    this: manganite carries a single log K at 25 C, so it is in the static
    25 C diagram and cannot be in an interactive one. Such species are
    dropped from the whole grid and listed in ``dropped``, which the figures
    print, rather than appearing and disappearing as the slider moves.
    """
    built = []
    for temperature in temperatures:
        built.append(
            element_series(
                element,
                species=species,
                temperature_c=temperature,
                activity=activity,
                fixed=fixed,
                backend=backend,
                reference=reference,
                skip_out_of_range=True,
            )
        )

    common = set.intersection(*({entry.backend for entry in series} for series in built))
    reasons = {}
    for series in built:
        for name, why in series.skipped:
            reasons.setdefault(name, why)
    dropped = tuple(
        (name, reasons.get(name, "not available at every temperature"))
        for name in sorted({e.backend for s in built for e in s} - common)
    )
    if not common:
        raise ValueError(
            f"no species of {element} is available across {list(temperatures)} °C: "
            + "; ".join(f"{name}: {why}" for name, why in dropped)
        )

    import dataclasses

    built = [
        dataclasses.replace(
            series, entries=tuple(e for e in series if e.backend in common), skipped=dropped
        )
        for series in built
    ]
    return ElementGrid(
        element=element,
        temperatures=tuple(float(t) for t in temperatures),
        series=tuple(built),
        activity=activity,
        dropped=dropped,
    )
