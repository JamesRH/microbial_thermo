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
        for entry in self._species:
            for key in self._keys_for(entry):
                # First registration wins, so the canonical backend name is not
                # displaced by an alias belonging to another entry.
                self._index.setdefault(key, entry)

    @staticmethod
    def _keys_for(entry: Species) -> list[str]:
        keys = [entry.backend, entry.formula, *entry.aliases]
        return [_normalise(k) for k in keys if k]

    def resolve(self, name: str) -> Species:
        """Look up ``name``; raises with suggestions when it is unknown."""
        if isinstance(name, Species):
            return name
        key = _normalise(name)
        if key in self._index:
            return self._index[key]
        raise SpeciesNotFoundError(name, suggestions=self.suggest(name))

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
        )
        for item in document["species"]
    ]
    return SpeciesRegistry(entries)


def resolve(name: str) -> Species:
    """Resolve ``name`` against the default registry."""
    return default_registry().resolve(name)
