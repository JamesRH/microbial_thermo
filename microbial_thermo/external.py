"""Thermodynamic data imported wholesale from a published database.

Distinct from :mod:`~microbial_thermo.supplemental`, which holds numbers typed
in by hand for species no database carries. These arrive from a source that
already cites its own primaries, so they come with provenance rather than
needing it added.

They are still not the primary database. Every source here is consulted after
``speq23``, after the GWB minerals and after ``supcrtbl``, so importing one can
only ever *add* a species -- it can never move a number that already resolved.
That ordering is what keeps a single provenance per result, and it is the
reason a source can be added without revalidating everything downstream.

A source only belongs here if its standard state already matches: aqueous,
1 molal, 25 C, dGf(H+) = 0, and the same revised-HKF equation of state. OBIGT
qualifies because it is SUPCRT-lineage, which is checkable rather than
asserted -- acetate's HKF parameters are identical in OBIGT and in speq23 to
the last digit.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

import yaml

_PACKAGE = "microbial_thermo.data.external"
_MANIFEST = "sources.yaml"

#: The order of an aqueous HKF record as pyGCC expects it, after the formula
#: and reference string: G, H, S, then the eight equation-of-state parameters.
#: The first three are named for their unit, so a file says what it holds.
_ENTRY_COLUMNS = {
    "cal": ("G_cal", "H_cal", "S_cal"),
    "J": ("G_J", "H_J", "S_J"),
}
_HKF_COLUMNS = ("a1", "a2", "a3", "a4", "c1", "c2", "omega")

#: pyGCC works in calories. OBIGT applies a row's energy unit to the HKF
#: coefficients as well as to G, H and S -- every one of a1-a4, c1, c2 and
#: omega carries an energy in its unit -- so a joule row converts all ten.
#: Only the charge does not. Leaving a coefficient in joules would be wrong
#: only away from 25 C, which is why tests/test_external.py checks the
#: converted coefficients against the source's own heat capacity.
_TO_CALORIES = {"cal": 1.0, "J": 1.0 / 4.184}


@dataclass(frozen=True)
class ExternalSource:
    """One imported database, and what it is."""

    name: str
    file: str
    kind: str
    energy_units: str
    citation: str
    url: str = ""
    license: str = ""
    retrieved: str = ""
    note: str = ""


@dataclass(frozen=True)
class ExternalSpecies:
    """One species imported from an external source."""

    name: str
    formula: str
    entry: tuple
    source: str
    citation: str
    ref_key: str = ""

    @property
    def provenance(self) -> str:
        primary = f" -- {self.citation}" if self.citation else ""
        return f"{self.source}{primary}"


@lru_cache(maxsize=1)
def external_sources() -> tuple:
    text = resources.files(_PACKAGE).joinpath(_MANIFEST).read_text()
    return tuple(
        ExternalSource(
            name=item["name"],
            file=item["file"],
            kind=item["kind"],
            energy_units=item["energy_units"],
            citation=" ".join(item.get("citation", "").split()),
            url=item.get("url", ""),
            license=item.get("license", ""),
            retrieved=str(item.get("retrieved", "")),
            note=" ".join(item.get("note", "").split()),
        )
        for item in yaml.safe_load(text)["sources"]
    )


@lru_cache(maxsize=1)
def external_species() -> dict:
    """Every imported species, keyed by name.

    Earlier sources in the manifest win, so the file order is the precedence
    order -- stated here because it is the kind of thing that otherwise gets
    discovered by accident.
    """
    out: dict[str, ExternalSpecies] = {}
    for source in external_sources():
        if source.kind != "hkf_aqueous":
            continue
        text = resources.files(_PACKAGE).joinpath(source.file).read_text()
        for row in csv.DictReader(text.splitlines()):
            if row["name"] in out:
                continue
            # pyGCC wants [formula, reference, G, H, S, a1..a4, c1, c2, omega, z].
            entry = [row["formula"], f" ref:{row.get('ref_key', '')} {row.get('date', '')}"]
            entry += _entry_values(row, source.energy_units)
            out[row["name"]] = ExternalSpecies(
                name=row["name"],
                formula=row["formula"],
                entry=tuple(entry),
                source=source.name,
                citation=row.get("citation", ""),
                ref_key=row.get("ref_key", ""),
            )
    return out


def _entry_values(row: dict, energy_units: str) -> list:
    """G, H, S, the seven HKF coefficients and the charge, in calories."""
    if energy_units not in _TO_CALORIES:
        raise ValueError(f"unknown energy_units {energy_units!r} in {_MANIFEST}")
    scale = _TO_CALORIES[energy_units]
    columns = _ENTRY_COLUMNS[energy_units] + _HKF_COLUMNS
    return [float(row[column]) * scale for column in columns] + [float(row["z"])]


def source_named(name: str) -> ExternalSource:
    for source in external_sources():
        if source.name == name:
            return source
    raise KeyError(name)
