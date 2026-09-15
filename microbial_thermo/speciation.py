"""Acid-base speciation: which form of a species exists at a given pH.

You name a species -- ``H2S``, ``sulfide``, ``HS-`` -- and this module works out
what fraction of the family each member holds at the working pH and
temperature, so a measured *total* can be turned into the activity of the
particular species a reaction is written with.

This matters because analytical measurements are reported as totals: total
sulfide, dissolved inorganic carbon, total ammonia. A reaction, by contrast, is
written with one specific form. At pH 7 sulfide is close to an even split
between H2S and HS-, so treating a measured total as though it were all one
form is wrong by a factor of about two, and choosing "the dominant species" is
very nearly a coin flip.

p*K*a values are **computed from the backend at the working temperature**
rather than read from a table, so a calculation at 60 C uses 60 C p*K*a values.
At 25 C the computed values reproduce the literature closely: H2S 6.99,
carbonate 6.34 and 10.33, ammonium 9.24, phosphate 2.17 / 7.21 / 12.32,
acetate 4.76.

The deprotonation step is balanced with water where the two forms differ in
oxygen as well as hydrogen. That is essential for the carbonate system, whose
first step is :math:`CO_2(aq) + H_2O \\rightleftharpoons HCO_3^- + H^+` rather
than a bare proton loss; treating it as one would put p*K*a1 tens of units out.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

import yaml

from .exceptions import MicrobialThermoError
from .species import Species, default_registry
from .units import KJ_PER_MOL_STR, R, as_magnitude, celsius_to_kelvin, ureg

_DATA_FILE = "acid_base.yaml"

#: Within this many pH units of a pKa, both forms are present in comparable
#: amounts and the choice of a single species is unreliable.
NEAR_PKA_UNITS = 1.0


@dataclass(frozen=True)
class AcidBaseFamily:
    """A set of species differing only by protonation, most protonated first."""

    name: str
    members: tuple
    aliases: tuple = ()

    def __len__(self) -> int:
        return len(self.members)

    def __iter__(self):
        return iter(self.members)

    @property
    def most_protonated(self) -> Species:
        return self.members[0]

    @property
    def least_protonated(self) -> Species:
        return self.members[-1]

    def contains(self, species) -> bool:
        backend = getattr(species, "backend", species)
        return any(m.backend == backend for m in self.members)


@lru_cache(maxsize=1)
def default_families() -> tuple:
    """Families shipped with the package, skipping any absent species."""
    text = resources.files("microbial_thermo.data").joinpath(_DATA_FILE).read_text()
    document = yaml.safe_load(text)
    registry = default_registry()

    families = []
    for item in document["families"]:
        members = []
        for name in item["members"]:
            resolved = registry.get(name)
            if resolved is not None:
                members.append(resolved)
        if len(members) >= 2:
            families.append(
                AcidBaseFamily(
                    name=item["name"],
                    members=tuple(members),
                    aliases=tuple(str(a) for a in item.get("aliases", [])),
                )
            )
    return tuple(families)


def family_for(species, families=None) -> AcidBaseFamily | None:
    """The family a species belongs to, or None."""
    registry = default_registry()
    resolved = registry.get(species) if not isinstance(species, Species) else species
    if resolved is None:
        return None
    for family in families or default_families():
        if family.contains(resolved):
            return family
    return None


def family_by_name(name: str, families=None) -> AcidBaseFamily | None:
    """Look a family up by its name or any alias, or by one of its members."""
    key = str(name).strip().lower()
    for family in families or default_families():
        if family.name.lower() == key:
            return family
        if any(alias.lower() == key for alias in family.aliases):
            return family
    return family_for(name, families)


def deprotonation_water(acid: Species, base: Species) -> int:
    """Water coefficient for ``acid -> base + H+ + w H2O``.

    Negative means water sits on the left, which is the carbonate case:
    ``CO2(aq) + H2O -> HCO3- + H+``.
    """
    return acid.parsed.element_count("O") - base.parsed.element_count("O")


def _check_step(acid: Species, base: Species) -> int:
    """Validate that two members really differ by exactly one proton."""
    if base.charge - acid.charge != -1:
        raise MicrobialThermoError(
            f"{acid.backend} and {base.backend} do not differ by one proton "
            f"(charges {acid.charge} and {base.charge})"
        )
    water = deprotonation_water(acid, base)
    expected = base.parsed.element_count("H") + 1 + 2 * water
    if acid.parsed.element_count("H") != expected:
        raise MicrobialThermoError(
            f"the step {acid.backend} -> {base.backend} does not balance for "
            "hydrogen; check the family definition"
        )
    return water


def pKa(acid, base, temperature_c: float = 25.0, pressure_bar=None, backend=None) -> float:
    """p*K*a of the step ``acid -> base + H+``, computed from the backend.

    Balanced with water where the two forms differ in oxygen.
    """
    from . import get_backend

    backend = backend or get_backend()
    registry = default_registry()
    acid = registry.resolve(acid) if not isinstance(acid, Species) else acid
    base = registry.resolve(base) if not isinstance(base, Species) else base
    water_coefficient = _check_step(acid, base)

    def gibbs(name: str) -> float:
        return as_magnitude(backend.delta_Gf(name, temperature_c, pressure_bar), KJ_PER_MOL_STR)

    delta_g = gibbs(base.backend) - gibbs(acid.backend)
    if water_coefficient:
        delta_g += water_coefficient * gibbs("H2O")
    # dGf(H+) is zero by convention and so contributes nothing.

    rt = as_magnitude(
        (R * (celsius_to_kelvin(temperature_c) * ureg.kelvin)).to(KJ_PER_MOL_STR),
        KJ_PER_MOL_STR,
    )
    return delta_g / (rt * math.log(10.0))


def pKa_ladder(family, temperature_c: float = 25.0, pressure_bar=None, backend=None) -> list[float]:
    """Successive p*K*a values for a family, one per consecutive pair."""
    family = _as_family(family)
    return [
        pKa(acid, base, temperature_c, pressure_bar, backend)
        for acid, base in zip(family.members, family.members[1:], strict=False)
    ]


def fractions(
    family, pH: float, temperature_c: float = 25.0, pressure_bar=None, backend=None
) -> dict:
    """Fraction of the family held by each member at ``pH``.

    Uses the standard polyprotic distribution: relative abundances are built up
    as :math:`r_j = r_{j-1} K_{a,j} / [H^+]`, then normalised. Returns a dict
    keyed by :class:`~microbial_thermo.species.Species`, summing to 1.
    """
    family = _as_family(family)
    ladder = pKa_ladder(family, temperature_c, pressure_bar, backend)

    # Work in log space: at low pH the raw ratios overflow for a triprotic acid.
    log_relative = [0.0]
    for value in ladder:
        log_relative.append(log_relative[-1] + (pH - value))
    largest = max(log_relative)
    relative = [10.0 ** (entry - largest) for entry in log_relative]
    total = sum(relative)
    return {member: value / total for member, value in zip(family.members, relative, strict=True)}


def fraction_of(
    species, pH: float, temperature_c: float = 25.0, pressure_bar=None, backend=None
) -> float:
    """Fraction of its family that ``species`` holds at ``pH``.

    Returns 1.0 for a species belonging to no family, so callers can apply this
    unconditionally.
    """
    family = family_for(species)
    if family is None:
        return 1.0
    registry = default_registry()
    resolved = registry.resolve(species) if not isinstance(species, Species) else species
    distribution = fractions(family, pH, temperature_c, pressure_bar, backend)
    for member, value in distribution.items():
        if member.backend == resolved.backend:
            return value
    return 1.0


def dominant(family, pH: float, temperature_c: float = 25.0, pressure_bar=None, backend=None):
    """The most abundant member at ``pH``, with its fraction."""
    distribution = fractions(family, pH, temperature_c, pressure_bar, backend)
    member = max(distribution, key=distribution.get)
    return member, distribution[member]


def near_pKa(
    family,
    pH: float,
    temperature_c: float = 25.0,
    pressure_bar=None,
    backend=None,
    tolerance: float = NEAR_PKA_UNITS,
) -> bool:
    """True when ``pH`` is close enough to a p*K*a that both forms matter."""
    ladder = pKa_ladder(family, temperature_c, pressure_bar, backend)
    return any(abs(pH - value) <= tolerance for value in ladder)


def warn_if_ambiguous(
    species, pH: float, temperature_c: float = 25.0, pressure_bar=None, backend=None
) -> None:
    """Warn when a single species is being used near a p*K*a."""
    family = family_for(species)
    if family is None:
        return
    if not near_pKa(family, pH, temperature_c, pressure_bar, backend):
        return
    member, value = dominant(family, pH, temperature_c, pressure_bar, backend)
    ladder = pKa_ladder(family, temperature_c, pressure_bar, backend)
    warnings.warn(
        f"pH {pH:g} is within {NEAR_PKA_UNITS:g} unit of a pKa for the "
        f"{family.name} family (pKa {', '.join(f'{v:.2f}' for v in ladder)}); "
        f"{member.backend} holds only {value:.0%} of the total. Give a total "
        f"concentration via total_concentrations={{'{family.name}': ...}} so "
        "the speciation is handled for you.",
        UserWarning,
        stacklevel=3,
    )


def speciation_table(
    family, pH: float, temperature_c: float = 25.0, pressure_bar=None, backend=None
):
    """Family distribution as a pandas DataFrame, for teaching and inspection."""
    import pandas as pd

    family = _as_family(family)
    distribution = fractions(family, pH, temperature_c, pressure_bar, backend)
    ladder = pKa_ladder(family, temperature_c, pressure_bar, backend)
    return pd.DataFrame(
        {
            "species": [m.backend for m in family.members],
            "fraction": [distribution[m] for m in family.members],
            "percent": [100.0 * distribution[m] for m in family.members],
            "pKa to next": [*[f"{v:.2f}" for v in ladder], ""],
        }
    )


def _as_family(family) -> AcidBaseFamily:
    if isinstance(family, AcidBaseFamily):
        return family
    resolved = family_by_name(str(family))
    if resolved is None:
        raise MicrobialThermoError(
            f"no acid-base family called {family!r}; known families: "
            + ", ".join(f.name for f in default_families())
        )
    return resolved
