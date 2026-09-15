"""The curated library of microbial metabolisms.

A named catalogue of the reactions a course actually covers, so a student can
write ``library.reaction("anammox")`` instead of looking up which couples to
pair. It also gives the test suite a large corpus: every entry is balanced and
put through the two-path cross-check.

Nothing thermodynamic is stored here. Potentials and free energies are computed
from the backend at whatever conditions are asked for, so the catalogue stays
correct at any pH and temperature rather than being a table pinned to pH 7.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources

import yaml

from .conditions import Conditions
from .exceptions import MicrobialThermoError

_DATA_FILE = "reactions.yaml"


@dataclass(frozen=True)
class Metabolism:
    """One named metabolism: a donor couple, an acceptor couple, and context."""

    name: str
    label: str
    group: str
    donor: tuple
    acceptor: tuple
    key_element: str | None = None
    organisms: str = ""
    note: str = ""

    def build(
        self,
        conditions: Conditions | None = None,
        backend=None,
        n_electrons: int = 2,
        normalize_to=None,
    ):
        """Construct the :class:`~microbial_thermo.reaction.Reaction`."""
        from .reaction import Couple, Reaction

        donor = Couple.make(self.donor[0], self.donor[1], self._key_element_for(self.donor))
        acceptor = Couple.make(
            self.acceptor[0], self.acceptor[1], self._key_element_for(self.acceptor)
        )
        return Reaction.from_couples(
            donor=donor,
            acceptor=acceptor,
            conditions=conditions,
            backend=backend,
            n_electrons=n_electrons,
            normalize_to=normalize_to,
        )

    def _key_element_for(self, couple) -> str | None:
        """Apply ``key_element`` only to the couple it actually describes.

        A metabolism usually needs the hint for just one of its two couples --
        carbon, for an organic donor driving an inorganic acceptor. Passing it
        to the other couple would assert an element that is not in it.
        """
        if self.key_element is None:
            return None
        from .balance import normalize_side

        # A side may name several species, so normalise before asking whether
        # the element is present throughout.
        in_both = all(
            any(self.key_element in species.parsed.elements for species, _ in normalize_side(side))
            for side in couple
        )
        return self.key_element if in_both else None

    def __str__(self) -> str:
        return f"{self.name} ({self.label})"


@dataclass
class MetabolismLibrary:
    """Lookup over the curated metabolisms."""

    entries: tuple = field(default_factory=tuple)

    def __iter__(self):
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)

    def __contains__(self, name: str) -> bool:
        return any(e.name == name for e in self.entries)

    def get(self, name: str) -> Metabolism:
        for entry in self.entries:
            if entry.name == name:
                return entry
        import difflib

        suggestions = difflib.get_close_matches(
            name, [e.name for e in self.entries], n=3, cutoff=0.4
        )
        message = f"no metabolism called {name!r}"
        if suggestions:
            message += f"; did you mean: {', '.join(suggestions)}?"
        raise MicrobialThermoError(message)

    def groups(self) -> list[str]:
        seen: list[str] = []
        for entry in self.entries:
            if entry.group not in seen:
                seen.append(entry.group)
        return seen

    def in_group(self, group: str) -> list[Metabolism]:
        return [e for e in self.entries if e.group == group]

    def names(self) -> list[str]:
        return [e.name for e in self.entries]


@lru_cache(maxsize=1)
def default_library() -> MetabolismLibrary:
    """The catalogue shipped with the package."""
    text = resources.files("microbial_thermo.data").joinpath(_DATA_FILE).read_text()
    document = yaml.safe_load(text)
    entries = tuple(
        Metabolism(
            name=item["name"],
            label=item.get("label", item["name"]),
            group=item.get("group", "other"),
            donor=tuple(item["donor"]),
            acceptor=tuple(item["acceptor"]),
            key_element=item.get("key_element"),
            organisms=item.get("organisms", ""),
            note=item.get("note", ""),
        )
        for item in document["metabolisms"]
    )
    return MetabolismLibrary(entries=entries)


def metabolism(name: str) -> Metabolism:
    """Look up one metabolism by name."""
    return default_library().get(name)


def reaction(
    name: str,
    conditions: Conditions | None = None,
    backend=None,
    n_electrons: int = 2,
    normalize_to=None,
):
    """Build a named metabolism straight into a reaction."""
    return (
        default_library()
        .get(name)
        .build(
            conditions=conditions,
            backend=backend,
            n_electrons=n_electrons,
            normalize_to=normalize_to,
        )
    )


def energy_table(
    conditions: Conditions | None = None,
    backend=None,
    n_electrons: int = 2,
    group: str | None = None,
    sort: bool = True,
):
    """Free energies of every curated metabolism, as a pandas DataFrame.

    Normalised per the given electron count so the rows are comparable. This is
    the data behind an affinity ladder: at one measured condition it shows which
    metabolisms pay and which do not.

    Entries that cannot be evaluated under these conditions -- an isothermal
    mineral away from 25 C, say -- are reported with NaN rather than dropped, so
    the omission is visible.
    """
    import pandas as pd

    library = default_library()
    entries = library.in_group(group) if group else list(library)

    rows = []
    for entry in entries:
        record = {
            "name": entry.name,
            "label": entry.label,
            "group": entry.group,
            "dG (kJ/mol)": float("nan"),
            "dG per e- (kJ/mol)": float("nan"),
            "dE (V)": float("nan"),
            "reaction": "",
            "problem": "",
        }
        try:
            built = entry.build(conditions, backend, n_electrons=n_electrons)
            record["dG (kJ/mol)"] = built.delta_G.magnitude
            record["dG per e- (kJ/mol)"] = built.delta_G_per_electron.magnitude
            record["dE (V)"] = built.delta_E.to("V").magnitude
            record["reaction"] = built.format()
        except Exception as exc:
            record["problem"] = f"{type(exc).__name__}: {exc}"
        rows.append(record)

    frame = pd.DataFrame(rows)
    if sort:
        frame = frame.sort_values("dG per e- (kJ/mol)", na_position="last")
    return frame.reset_index(drop=True)
