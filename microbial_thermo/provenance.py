"""Where a number came from.

A free energy computed here rests on four things: the software that did the
arithmetic, the database the formation energies were read from, the
conventions applied on the way through, and the conditions asked for. A result
quoted without those is not reproducible, and in a teaching tool it is not
checkable either.

This module collects all four into one record. It is deliberately boring:
version strings, a file hash, a per-species source line, and the standard-state
choices. The point is that a student can paste it under a figure and someone
else can tell whether they would get the same number.

The database hash matters more than it looks. pyGCC ships several SUPCRT and
GWB files, they are revised between releases, and nothing in a computed value
records which one was read. A hash pins it.
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field
from datetime import UTC, datetime

from .units import KJ_PER_MOL_STR, as_magnitude

#: How much of a database file to hash. These run to tens of megabytes and
#: hashing all of it on every result is wasteful; the first megabyte is more
#: than enough to distinguish revisions, and the size is recorded alongside.
_HASH_BYTES = 1 << 20


@dataclass(frozen=True)
class SpeciesSource:
    """Where one species' formation energy came from."""

    name: str
    backend_name: str
    source: str
    verified: bool = True

    def __str__(self) -> str:
        flag = "" if self.verified else "  [UNVERIFIED]"
        return f"{self.name} <- {self.source}{flag}"


@dataclass(frozen=True)
class DatabaseFile:
    """A database file, identified well enough to tell revisions apart."""

    name: str
    path: str | None = None
    sha256_prefix: str | None = None
    size_bytes: int | None = None

    @classmethod
    def describe(cls, path: str | None, name: str) -> DatabaseFile:
        if not path or not os.path.exists(path):
            return cls(name=name)
        with open(path, "rb") as handle:
            digest = hashlib.sha256(handle.read(_HASH_BYTES)).hexdigest()
        return cls(
            name=name,
            path=path,
            sha256_prefix=digest[:16],
            size_bytes=os.path.getsize(path),
        )

    def __str__(self) -> str:
        if self.sha256_prefix is None:
            return f"{self.name} (not located on disk)"
        return f"{self.name}  sha256:{self.sha256_prefix}…  {self.size_bytes} bytes"


@dataclass(frozen=True)
class Provenance:
    """Everything needed to reproduce one computed result."""

    equation: str
    delta_G_standard_prime: float
    n_electrons: float
    conditions: str
    versions: dict
    databases: tuple
    species: tuple
    conventions: tuple
    generated: str = field(default_factory=lambda: datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC"))

    @property
    def unverified(self) -> list:
        return [s for s in self.species if not s.verified]

    def to_dict(self) -> dict:
        return {
            "equation": self.equation,
            "delta_G_standard_prime_kJ_per_mol": self.delta_G_standard_prime,
            "n_electrons": self.n_electrons,
            "conditions": self.conditions,
            "generated": self.generated,
            "versions": dict(self.versions),
            "databases": [
                {
                    "name": d.name,
                    "sha256_prefix": d.sha256_prefix,
                    "size_bytes": d.size_bytes,
                }
                for d in self.databases
            ],
            "species": [
                {
                    "name": s.name,
                    "backend_name": s.backend_name,
                    "source": s.source,
                    "verified": s.verified,
                }
                for s in self.species
            ],
            "conventions": list(self.conventions),
        }

    def to_text(self) -> str:
        lines = [
            "Provenance",
            "=" * 60,
            f"reaction   {self.equation}",
            f"dG0'       {self.delta_G_standard_prime:+.2f} kJ/mol per {self.n_electrons:g} e-",
            f"conditions {self.conditions}",
            f"generated  {self.generated}",
            "",
            "Software",
        ]
        lines += [f"  {k}: {v}" for k, v in self.versions.items()]
        lines += ["", "Databases"]
        lines += [f"  {d}" for d in self.databases]
        lines += ["", "Species"]
        lines += [f"  {s}" for s in self.species]
        lines += ["", "Conventions"]
        lines += [f"  {c}" for c in self.conventions]
        if self.unverified:
            names = ", ".join(s.name for s in self.unverified)
            lines += [
                "",
                "WARNING: rests on hand-entered, unverified values for " + names + ".",
            ]
        return "\n".join(lines)

    def to_bibtex(self, key: str | None = None) -> str:
        """A citable record of the calculation.

        ``@misc`` rather than ``@software``: what is being cited is this
        particular calculation, of which the software is one input.
        """
        key = key or _bibtex_key(self.equation)
        database = "; ".join(str(d) for d in self.databases)
        note = (
            f"Computed with microbial_thermo "
            f"{self.versions.get('microbial_thermo', 'unknown')} "
            f"on pyGCC {self.versions.get('pygcc', 'unknown')}. "
            f"Databases: {database}. Conditions: {self.conditions}. "
            f"{'; '.join(self.conventions)}."
        )
        if self.unverified:
            note += (
                " Depends on hand-entered, unverified formation energies for "
                + ", ".join(s.name for s in self.unverified)
                + "."
            )
        year = self.generated[:4]
        return "\n".join(
            [
                f"@misc{{{key},",
                f"  title  = {{{{{self.equation}}}}},",
                "  author = {{microbial\\_thermo}},",
                f"  year   = {{{year}}},",
                f"  note   = {{{_escape(note)}}}",
                "}",
            ]
        )


def _escape(text: str) -> str:
    for char in ("&", "%", "#", "_"):
        text = text.replace(char, "\\" + char)
    return text


def _bibtex_key(equation: str) -> str:
    digest = hashlib.sha256(equation.encode()).hexdigest()[:8]
    return f"mthermo{digest}"


#: Conventions that shape every number this library produces. Recorded because
#: they are exactly what a reader comparing against another source needs to
#: know, and exactly what nobody writes down.
def _conventions(reaction) -> tuple:
    conditions = reaction.conditions
    return (
        "standard state: dGf(H+) = 0 and dGf(e-) = 0, i.e. the SHE scale",
        "dG0' is the species-level (microbial bioenergetics) convention with "
        "explicit protons, not the Legendre-transformed biochemical convention",
        f"activity model: {conditions.activity_model}",
        f"water: {getattr(reaction.backend, 'water_convention', 'supcrt')} datum",
        "couples written (reduced, oxidized); half reactions stored as reductions",
    )


def provenance_for(reaction) -> Provenance:
    """Build the full provenance record for a computed reaction."""
    backend = reaction.backend
    conditions = reaction.conditions

    sources = []
    for species in sorted(reaction.coefficients, key=lambda s: s.backend):
        try:
            record = backend.resolve(species.backend)
            sources.append(
                SpeciesSource(
                    name=species.backend,
                    backend_name=record.backend_name,
                    source=record.source,
                    verified=record.verified,
                )
            )
        except Exception as exc:  # a species with no describable origin
            sources.append(
                SpeciesSource(
                    name=species.backend,
                    backend_name=species.backend,
                    source=f"unresolved ({type(exc).__name__})",
                    verified=False,
                )
            )

    return Provenance(
        equation=reaction.format(),
        delta_G_standard_prime=as_magnitude(reaction.delta_G_standard_prime, KJ_PER_MOL_STR),
        n_electrons=float(reaction.n_electrons),
        conditions=(
            f"pH {conditions.pH:g}, {conditions.temperature_c:g} °C, "
            f"activity model {conditions.activity_model}"
        ),
        versions=_versions(),
        databases=_databases(backend),
        species=tuple(sources),
        conventions=_conventions(reaction),
    )


def _versions() -> dict:
    from . import versions

    return versions()


def _databases(backend) -> tuple:
    from .backends.gwb import default_gwb_path

    hkf_path = getattr(backend, "_database", None) or _default_hkf_path()
    gwb_path = getattr(backend, "_mineral_database", None) or default_gwb_path()
    return (
        DatabaseFile.describe(hkf_path, getattr(backend, "database_name", "unknown")),
        DatabaseFile.describe(gwb_path, getattr(backend, "mineral_database_name", "unknown")),
    )


def _default_hkf_path() -> str | None:
    try:
        import pygcc
    except ImportError:
        return None
    return os.path.join(os.path.dirname(pygcc.__file__), "default_db", "speq21.dat")
