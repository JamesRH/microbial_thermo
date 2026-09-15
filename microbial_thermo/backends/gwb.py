"""Mineral data from a GWB-format thermodynamic database.

pyGCC's direct-access databases (``speq21.dat``, ``supcrtbl.dat``) carry
equation-of-state parameters and cover aqueous solutes and a few hundred
minerals. They do **not** contain the manganese oxides that matter for
microbial metal cycling: pyrolusite, birnessite, manganite, hausmannite.

Those live in the GWB-format databases pyGCC also ships (``thermo.com.dat``),
but in a different form. A GWB database stores, for each mineral, a
dissolution reaction into basis species together with tabulated equilibrium
constants on a fixed temperature grid, rather than equation-of-state
parameters. For example::

    Pyrolusite      formula= MnO2
         2 species in reaction
      0.5000 Mn++              0.5000 MnO4--
           -19.0122  -17.6439  -16.1759  -14.9379
           -13.8546  -13.1548  -12.7741  -12.7231

which is the reaction ``Pyrolusite -> 0.5 Mn(2+) + 0.5 MnO4(2-)`` with
``log K`` at 0, 25, 60, 100, 150, 200, 250 and 300 degrees Celsius.

Since

.. math::
    \\Delta G_\\mathrm{rxn} = \\sum_i \\nu_i \\Delta G_{f,i}
        - \\Delta G_f(\\mathrm{mineral}) = -RT\\ln(10)\\log K

the mineral's formation energy follows from the basis species, which the HKF
database does provide:

.. math::
    \\Delta G_f(\\mathrm{mineral}) = \\sum_i \\nu_i \\Delta G_{f,i}
        + RT\\ln(10)\\,\\log K

**Validation.** Goethite, hematite and magnetite appear in both the GWB
database and ``supcrtbl.dat``, giving two independent routes to the same
number. They agree with each other, and with published values, to within
1-3 kJ/mol:

=========== ============= =============== ============
mineral     this route    supcrtbl HP11   published
=========== ============= =============== ============
Goethite    -488.5        -491.0          -490.6
Hematite    -745.3        -743.7          -742.7
Magnetite   -1014.8       -1011.6         -1012.6
=========== ============= =============== ============

So the manganese oxides, which have no second route, can be trusted to about
the same few kJ/mol.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from ..exceptions import MissingDataError, OutOfRangeError

#: GWB writes this in place of a log K that could not be computed.
_NO_DATA_SENTINEL = 500.0

#: How far from a single tabulated point an isothermal entry may be used.
_ISOTHERMAL_TOLERANCE_C = 0.5

#: Gas constant in cal/(mol K), matching pyGCC's internal energy units.
_R_CAL = 1.987204

_SPECIES_COUNT = re.compile(r"^\s*(\d+)\s+species in reaction", re.I)
_FORMULA = re.compile(r"^\s*formula=\s*(\S+)")
_NUMBER = re.compile(r"[-+]?\d*\.\d+|[-+]?\d+")

#: A name line starts at column zero and is not a comment or a keyword.
_NAME_LINE = re.compile(r"^(?!\s|\*|-)(\S.*?)\s*(?:type=.*)?$")


@dataclass(frozen=True)
class MineralEntry:
    """One mineral's dissolution reaction and its log K grid."""

    name: str
    formula: str
    #: ``[(coefficient, basis species name), ...]`` for the dissolution reaction.
    stoichiometry: tuple
    temperatures_c: np.ndarray
    log_k: np.ndarray

    @property
    def is_isothermal(self) -> bool:
        """True when the database gives only a single usable log K.

        Several manganese oxides -- manganite and birnessite among them -- are
        tabulated at 25 C only, every other slot being the no-data sentinel.
        """
        return int(np.count_nonzero(np.isfinite(self.log_k))) == 1

    def log_k_at(self, temperature_c: float) -> float:
        """Interpolate log K, ignoring the no-data sentinel.

        Entries with a single usable value are honoured at that temperature
        and refused elsewhere. Extrapolating one point across a temperature
        range would invent a temperature dependence the database does not
        contain, which is worse than failing.
        """
        valid = np.isfinite(self.log_k)
        temperatures = self.temperatures_c[valid]
        values = self.log_k[valid]
        if values.size == 0:
            raise MissingDataError(f"{self.name} has no usable log K values")
        if values.size == 1:
            only_t = float(temperatures[0])
            if abs(temperature_c - only_t) <= _ISOTHERMAL_TOLERANCE_C:
                return float(values[0])
            raise OutOfRangeError(
                f"{self.name} is tabulated at {only_t:g} C only, so its free "
                f"energy at {temperature_c:g} C is unknown. Evaluate at "
                f"{only_t:g} C, or choose a phase with a full log K grid."
            )
        if not (temperatures.min() <= temperature_c <= temperatures.max()):
            raise OutOfRangeError(
                f"{self.name} has no log K data at {temperature_c} C; "
                f"the grid spans {temperatures.min()}-{temperatures.max()} C"
            )
        # The grid is coarse (0, 25, 60, 100 within our range), so a cubic
        # spline tracks the curvature of log K with temperature far better
        # than a straight line between points.
        if values.size >= 4:
            from scipy.interpolate import CubicSpline

            return float(CubicSpline(temperatures, values)(temperature_c))
        return float(np.interp(temperature_c, temperatures, values))


def _parse_reaction_line(line: str) -> list:
    """Read alternating coefficient and species-name pairs from one line."""
    tokens = line.split()
    pairs = []
    index = 0
    while index + 1 < len(tokens):
        try:
            coefficient = float(tokens[index])
        except ValueError:
            break
        pairs.append((coefficient, tokens[index + 1]))
        index += 2
    return pairs


@lru_cache(maxsize=4)
def load_gwb_minerals(path: str) -> dict:
    """Parse a GWB-format database into :class:`MineralEntry` records.

    Deliberately tolerant: entries it cannot make sense of are skipped rather
    than raising, since a database holds thousands of phases and only a
    handful are ever needed.
    """
    with open(path, errors="ignore") as handle:
        lines = handle.read().splitlines()

    temperatures = _parse_temperature_grid(lines)
    entries: dict[str, MineralEntry] = {}

    index = 0
    while index < len(lines):
        match = _SPECIES_COUNT.match(lines[index])
        if not match:
            index += 1
            continue

        name = _find_name(lines, index)
        if name is None:
            index += 1
            continue
        formula = _find_formula(lines, index)

        count = int(match.group(1))
        stoichiometry: list = []
        cursor = index + 1
        while cursor < len(lines) and len(stoichiometry) < count:
            pairs = _parse_reaction_line(lines[cursor])
            if not pairs:
                break
            stoichiometry.extend(pairs)
            cursor += 1

        values: list[float] = []
        while cursor < len(lines) and len(values) < len(temperatures):
            numbers = [float(n) for n in _NUMBER.findall(lines[cursor])]
            if not numbers or lines[cursor].lstrip().startswith("*"):
                break
            values.extend(numbers)
            cursor += 1

        if len(stoichiometry) == count and len(values) >= len(temperatures):
            grid = np.array(values[: len(temperatures)], dtype=float)
            grid[grid >= _NO_DATA_SENTINEL] = np.nan
            entries[name] = MineralEntry(
                name=name,
                formula=formula,
                stoichiometry=tuple(stoichiometry),
                temperatures_c=temperatures,
                log_k=grid,
            )
        index = max(cursor, index + 1)

    return entries


def _parse_temperature_grid(lines) -> np.ndarray:
    """Read the ``* temperatures`` block from the file header."""
    for index, line in enumerate(lines):
        if line.strip().lower().startswith("* temperatures"):
            values: list[float] = []
            cursor = index + 1
            while cursor < len(lines) and not lines[cursor].lstrip().startswith("*"):
                values.extend(float(n) for n in _NUMBER.findall(lines[cursor]))
                cursor += 1
            if values:
                return np.array(values, dtype=float)
    # The standard GWB grid, used when a file omits the header block.
    return np.array([0.0, 25.0, 60.0, 100.0, 150.0, 200.0, 250.0, 300.0])


def _find_name(lines, index: int) -> str | None:
    """Walk back from the reaction-count line to the entry's name."""
    for cursor in range(index - 1, max(index - 8, -1), -1):
        stripped = lines[cursor]
        if not stripped.strip() or stripped.lstrip().startswith("*"):
            continue
        if stripped[:1].isspace():
            continue  # an indented property line such as "formula=" or "mole vol."
        match = _NAME_LINE.match(stripped)
        if match:
            return match.group(1).strip()
    return None


def _find_formula(lines, index: int) -> str:
    for cursor in range(index - 1, max(index - 8, -1), -1):
        match = _FORMULA.match(lines[cursor])
        if match:
            return match.group(1)
    return ""


def default_gwb_path() -> str:
    """Path to the GWB database bundled with pyGCC."""
    import pygcc

    return os.path.join(os.path.dirname(pygcc.__file__), "default_db", "thermo.com.dat")


def mineral_gibbs_cal(entry: MineralEntry, temperature_c: float, basis_gibbs_cal) -> float:
    """Formation energy of a mineral in cal/mol.

    ``basis_gibbs_cal`` is called with a basis species name and must return its
    formation energy in cal/mol at the same temperature and pressure.
    """
    log_k = entry.log_k_at(temperature_c)
    total = 0.0
    for coefficient, species in entry.stoichiometry:
        total += coefficient * basis_gibbs_cal(species)
    rt_ln10 = _R_CAL * (temperature_c + 273.15) * np.log(10.0)
    return float(total + rt_ln10 * log_k)
