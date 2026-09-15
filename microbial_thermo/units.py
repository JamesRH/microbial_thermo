"""Unit registry and physical constants.

One shared :mod:`pint` registry lives here. Multiple registries produce cryptic
cross-registry errors, so every module must import ``ureg`` from this file rather
than constructing its own.

Convention used throughout the library:

* ``pint`` quantities at the public API boundary (arguments, return values, axis
  labels).
* Bare ``float``/``ndarray`` on internal numeric hot paths, because ``scipy``
  routines do not accept ``Quantity`` objects.

Use :func:`as_magnitude` to cross that boundary explicitly rather than calling
``.magnitude`` ad hoc, so the expected unit is always stated.
"""

from __future__ import annotations

import pint

ureg = pint.UnitRegistry()
Quantity = ureg.Quantity

# --- Physical constants -------------------------------------------------------
# Taken from the pint registry so they carry CODATA values and correct units.
# Multiplied by 1 so these are Quantities, not bare Units -- a Unit has no
# .to() and the difference surfaces far from here.
R = 1 * ureg.molar_gas_constant
FARADAY = 1 * ureg.faraday_constant

#: pyGCC reports energies in thermochemical calories per mole.
CAL_PER_MOL = ureg.calorie / ureg.mol
KJ_PER_MOL = ureg.kilojoule / ureg.mol
KJ_PER_MOL_STR = "kJ/mol"

#: Absolute zero offset, kept explicit so temperature conversions are auditable.
KELVIN_OFFSET = 273.15

# --- Domain defaults ----------------------------------------------------------
#: Standard biochemical pH for the primed (E-standard-prime) convention.
STANDARD_BIOCHEMICAL_PH = 7.0

#: Free energy of ATP hydrolysis under typical cellular conditions. Used only to
#: express a result in "ATP equivalents"; it is a modelling choice, not a
#: measured constant, so it is exposed as an overridable default everywhere.
DEFAULT_DELTA_G_ATP = Quantity(-50.0, KJ_PER_MOL_STR)

#: Minimum free energy increment thought to sustain life (the "biological energy
#: quantum"). Roughly one third of an ATP; literature values span -10 to -20.
DEFAULT_BIOLOGICAL_ENERGY_QUANTUM = Quantity(-20.0, KJ_PER_MOL_STR)


def as_magnitude(value, unit: str) -> float:
    """Return ``value`` as a bare float in ``unit``.

    Accepts a plain number (assumed to already be in ``unit``) or a pint
    Quantity (converted). Crossing from the unit-aware API into bare numeric
    code should always go through here so the assumed unit is written down.
    """
    if isinstance(value, Quantity):
        return float(value.to(unit).magnitude)
    return float(value)


def to_quantity(value, unit: str) -> Quantity:
    """Coerce ``value`` to a Quantity in ``unit``, converting if already united."""
    if isinstance(value, Quantity):
        return value.to(unit)
    return Quantity(float(value), unit)


def celsius_to_kelvin(t_c: float) -> float:
    return t_c + KELVIN_OFFSET


def kelvin_to_celsius(t_k: float) -> float:
    return t_k - KELVIN_OFFSET
