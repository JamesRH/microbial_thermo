"""Species registry: user-facing names <-> pyGCC database names <-> formulas.

pyGCC inherits SUPCRT/GWB naming, where several species are named rather than
written as formulas (``Methane(aq)``, ``Acetate``, ``Lactate(aq)``) and charge
uses doubled signs (``SO4--``). Those strings are not parseable formulas, so
balancing and oxidation-state code cannot consume them directly.

This module keeps the three identities separate and explicit:

``backend``
    The exact string pyGCC expects.
``formula``
    A parseable chemical formula, used for balancing and oxidation states.
``display``
    A LaTeX-ish label for figures.
"""

from __future__ import annotations

import difflib
from collections.abc import Iterable
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources

import yaml

from .exceptions import SpeciesNotFoundError
from .formula import ParsedFormula, parse_formula

_DATA_FILE = "species.yaml"


@dataclass(frozen=True)
class Species:
    """One chemical species, in all three of its identities."""

    backend: str
    formula: str
    display: str = ""
    phase: str = "aq"
    smiles: str | None = None
    aliases: tuple[str, ...] = field(default_factory=tuple)
    #: Declares this entry the intended answer when a name it claims is also
    #: claimed by another phase -- a bare "H2" being the case that matters.
    #: Without it, which entry wins would depend on file ordering.
    canonical: bool = False

    @property
    def parsed(self) -> ParsedFormula:
        return parse_formula(self.formula)

    @property
    def charge(self) -> int:
        return self.parsed.charge

    @property
    def is_gas(self) -> bool:
        return self.phase == "g"

    @property
    def is_solid(self) -> bool:
        return self.phase == "s"

    @property
    def is_aqueous(self) -> bool:
        return self.phase == "aq"

    @property
    def label(self) -> str:
        return self.display or self.formula

    def __str__(self) -> str:
        return self.backend


class SpeciesRegistry:
    """Resolves any reasonable user string onto a :class:`Species`."""

    def __init__(self, species: Iterable[Species]):
        self._species: list[Species] = list(species)
        self._index: dict[str, Species] = {}
        #: Keys claimed by more than one species, and every claimant. A bare
        #: "H2" is claimed by both the dissolved and the gaseous entry.
        self._contested: dict[str, list[Species]] = {}

        for entry in self._species:
            for key in self._keys_for(entry):
                claimants = self._contested.setdefault(key, [])
                # One species can claim a key twice -- its formula and an
                # identical alias -- which is not an ambiguity.
                if entry.backend not in {e.backend for e in claimants}:
                    claimants.append(entry)

        for key, claimants in self._contested.items():
            if len(claimants) == 1:
                self._index[key] = claimants[0]
                continue
            declared = [e for e in claimants if e.canonical]
            if len(declared) == 1:
                self._index[key] = declared[0]
            else:
                # Nothing declared, so fall back on order -- but record it so
                # ambiguities() can report the guess rather than hiding it.
                self._index[key] = claimants[0]

        self._contested = {k: v for k, v in self._contested.items() if len(v) > 1}

    @staticmethod
    def _keys_for(entry: Species) -> list[str]:
        keys = [entry.backend, entry.formula, *entry.aliases]
        return [_normalise(k) for k in keys if k]

    @property
    def undeclared_ambiguities(self) -> dict:
        """Contested names where no entry claims to be canonical.

        These resolve by file ordering, which is not a decision anyone made.
        """
        return {
            key: list(claimants)
            for key, claimants in self._contested.items()
            if not any(e.canonical for e in claimants)
        }

    def ambiguities(self) -> dict:
        """Every name claimed by more than one species, with what it resolves to.

        ``{"h2": (resolved, [all claimants], declared)}`` -- ``declared`` says
        whether the winner was chosen deliberately or fell out of file order.
        """
        out = {}
        for key, claimants in self._contested.items():
            resolved = self._index[key]
            out[key] = (resolved, list(claimants), bool(resolved.canonical))
        return out

    def resolve(self, name, phase: str | None = None) -> Species:
        """Look up ``name``; raises with suggestions when it is unknown.

        ``phase`` picks between entries sharing a name -- ``resolve("H2",
        phase="g")`` for the gas rather than the dissolved form.
        """
        if isinstance(name, Species):
            return name
        key = _normalise(name)
        if key not in self._index:
            raise SpeciesNotFoundError(name, suggestions=self.suggest(name))

        if phase is None:
            return self._index[key]

        for entry in self._contested.get(key, [self._index[key]]):
            if entry.phase == phase:
                return entry
        available = sorted({e.phase for e in self._contested.get(key, [self._index[key]])})
        raise SpeciesNotFoundError(f"{name} (phase {phase!r}; available: {', '.join(available)})")

    def get(self, name: str, default=None):
        try:
            return self.resolve(name)
        except SpeciesNotFoundError:
            return default

    def suggest(self, name: str, n: int = 5) -> list[str]:
        matches = difflib.get_close_matches(_normalise(name), list(self._index), n=n, cutoff=0.6)
        # Report canonical backend names rather than whichever alias matched.
        seen, out = set(), []
        for m in matches:
            backend = self._index[m].backend
            if backend not in seen:
                seen.add(backend)
                out.append(backend)
        return out

    def __contains__(self, name: str) -> bool:
        return _normalise(name) in self._index

    def __iter__(self):
        return iter(self._species)

    def __len__(self) -> int:
        return len(self._species)

    def all_species(self) -> list[Species]:
        return list(self._species)


def _normalise(name: str) -> str:
    """Case- and space-insensitive lookup key."""
    return str(name).strip().lower().replace(" ", "_")


@lru_cache(maxsize=1)
def default_registry() -> SpeciesRegistry:
    """The registry shipped with the package."""
    text = resources.files("microbial_thermo.data").joinpath(_DATA_FILE).read_text()
    document = yaml.safe_load(text)
    entries = [
        Species(
            backend=item["backend"],
            formula=item["formula"],
            display=item.get("display", ""),
            phase=item.get("phase", "aq"),
            smiles=item.get("smiles"),
            aliases=tuple(str(a) for a in item.get("aliases", [])),
            canonical=bool(item.get("canonical", False)),
        )
        for item in document["species"]
    ]

    # Species no pyGCC database carries, supplied by the hand-entered table.
    from .supplemental import supplemental_species

    entries.extend(
        Species(
            backend=entry.backend,
            formula=entry.formula,
            display=entry.display,
            phase=entry.phase,
            smiles=entry.smiles,
            aliases=entry.aliases,
        )
        for entry in supplemental_species().values()
    )
    return SpeciesRegistry(entries)


def resolve(name, phase: str | None = None) -> Species:
    """Resolve ``name`` against the default registry.

    ``phase`` disambiguates a name two entries share: ``resolve("H2", "g")``.
    """
    return default_registry().resolve(name, phase=phase)
