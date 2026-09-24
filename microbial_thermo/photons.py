"""Light as an energy input.

A phototroph runs reactions that do not pay. Photoferrotrophy and
photoarsenotrophy both fix CO2 using a donor far too weak to do it in the
dark -- the reaction is uphill by tens to hundreds of kJ/mol -- and the
shortfall is made up by absorbed photons. This module supplies the photon
side of that ledger so the two can be put on one axis.

The arithmetic is trivial; the honesty is not. Three things to be clear about
before quoting a number from here:

**A photon's energy is not all available as work.** ``photon_energy`` returns
:math:`N_A hc/\\lambda`, which is an *upper bound*. Radiation carries entropy,
the reaction centre captures only part of the excitation, and some is lost as
heat within picoseconds. Treat a calculation at ``efficiency=1.0`` as asking
"could light conceivably cover this?", not "does it".

**Real quantum requirements are several photons per electron.** Anoxygenic
phototrophs run cyclic electron flow and reverse electron transport, so the
photons per CO2 fixed is an empirical number in the range of four to ten, not
something derivable here. ``photons_required`` answers the narrower
thermodynamic question -- how many photons' worth of energy the reaction is
short by -- and that is a floor, not a prediction.

**Adding photon energy to a Gibbs energy is a bookkeeping convenience.** It
treats light as a work input to an otherwise closed system. That is the
standard teaching treatment and it gets the sign and scale right, but a proper
account uses the excited-state redox potential of the reaction centre.
"""

from __future__ import annotations

from .exceptions import MicrobialThermoError
from .units import KJ_PER_MOL_STR, Quantity, ureg

#: Reaction-centre absorption maxima, nm. The names are the pigments'
#: conventional labels: P870 is bacteriochlorophyll a in purple bacteria,
#: which is what photoferrotrophs and the Mono Lake photoarsenotroph use.
REACTION_CENTRES = {
    "P680": 680.0,  # chlorophyll a, photosystem II
    "P700": 700.0,  # chlorophyll a, photosystem I
    "P840": 840.0,  # bacteriochlorophyll a, green sulfur bacteria
    "P870": 870.0,  # bacteriochlorophyll a, purple bacteria
    "P960": 960.0,  # bacteriochlorophyll b
}

#: The default: purple-bacterial bacteriochlorophyll a, the pigment of both
#: the photoferrotrophs and *Ectothiorhodospira* sp. PHS-1.
DEFAULT_WAVELENGTH_NM = REACTION_CENTRES["P870"]


def photon_energy(wavelength_nm: float | str = DEFAULT_WAVELENGTH_NM) -> Quantity:
    """Energy of one mole of photons at ``wavelength_nm``, as kJ/mol.

    Accepts a reaction-centre name as well as a number, so ``"P870"`` and
    ``870`` are the same request. This is the thermodynamic ceiling on what a
    photon can contribute, not what a cell gets from it.
    """
    nanometres = _as_nanometres(wavelength_nm)
    energy = (
        1
        * ureg.avogadro_constant
        * ureg.planck_constant
        * ureg.speed_of_light
        / (nanometres * ureg.nanometer)
    )
    return energy.to(KJ_PER_MOL_STR)


def _as_nanometres(wavelength_nm) -> float:
    if isinstance(wavelength_nm, str):
        try:
            return REACTION_CENTRES[wavelength_nm.upper()]
        except KeyError:
            known = ", ".join(sorted(REACTION_CENTRES))
            raise MicrobialThermoError(
                f"unknown reaction centre {wavelength_nm!r}; known: {known}"
            ) from None
    value = float(wavelength_nm)
    if value <= 0:
        raise MicrobialThermoError("wavelength must be positive")
    return value


def light_available(
    photons: float = 1.0,
    wavelength_nm: float | str = DEFAULT_WAVELENGTH_NM,
    efficiency: float = 1.0,
) -> Quantity:
    """Energy ``photons`` moles of photons could supply, in kJ/mol.

    ``efficiency`` scales the ideal photon energy. Leave it at 1 for the
    thermodynamic ceiling; set it to something like 0.3-0.5 to ask a question
    closer to what a reaction centre actually delivers.
    """
    if photons < 0:
        raise MicrobialThermoError("photon count cannot be negative")
    if not 0 < efficiency <= 1:
        raise MicrobialThermoError("efficiency must be in (0, 1]")
    return photons * efficiency * photon_energy(wavelength_nm)


def net_with_light(
    delta_g,
    photons: float = 1.0,
    wavelength_nm: float | str = DEFAULT_WAVELENGTH_NM,
    efficiency: float = 1.0,
) -> Quantity:
    """``delta_g`` after crediting the energy of ``photons`` moles of photons.

    Light is an input, so the credit is subtracted: an uphill reaction becomes
    less positive. Returns kJ/mol.
    """
    from .units import as_magnitude

    cost = as_magnitude(delta_g, KJ_PER_MOL_STR)
    credit = light_available(photons, wavelength_nm, efficiency)
    return Quantity(cost - credit.magnitude, KJ_PER_MOL_STR)


def photons_required(
    delta_g,
    wavelength_nm: float | str = DEFAULT_WAVELENGTH_NM,
    efficiency: float = 1.0,
    energy_quantum=None,
) -> float:
    """Photons needed to bring ``delta_g`` down to the energy quantum.

    Not to zero: a reaction sitting exactly at zero does no useful work, so
    the target is the biological energy quantum, the same threshold the
    figures shade. Returns 0 for a reaction that already pays, and a
    fractional count otherwise -- fractions are meaningful here because the
    question is "how many photons' worth of energy", not "how many photons".
    """
    from .units import DEFAULT_BIOLOGICAL_ENERGY_QUANTUM, as_magnitude

    cost = as_magnitude(delta_g, KJ_PER_MOL_STR)
    target = as_magnitude(energy_quantum or DEFAULT_BIOLOGICAL_ENERGY_QUANTUM, KJ_PER_MOL_STR)
    shortfall = cost - target
    if shortfall <= 0:
        return 0.0
    per_photon = light_available(1.0, wavelength_nm, efficiency).magnitude
    return shortfall / per_photon
