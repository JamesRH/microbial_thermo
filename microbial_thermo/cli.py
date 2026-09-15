"""Command line interface.

Every command is a thin wrapper over the library: argument parsing is
separated from the work so the same functions can be imported into a notebook
without the CLI getting in the way.
"""

from __future__ import annotations

import click

from . import configure, versions
from .conditions import Conditions
from .reaction import Couple, HalfReactionResult, Reaction
from .species import default_registry

#: Dependencies whose versions are reported by --version, for traceability.
_REPORTED = (
    "pygcc",
    "numpy",
    "pandas",
    "scipy",
    "chempy",
    "sympy",
    "pint",
    "matplotlib",
    "plotly",
)


def _print_versions(context, parameter, value):
    if not value or context.resilient_parsing:
        return
    table = versions()
    click.echo(f"microbial_thermo {table['microbial_thermo']}")
    click.echo("computational dependencies:")
    for package in _REPORTED:
        click.echo(f"  {package:<12} {table.get(package, 'not installed')}")
    context.exit()


_version_option = click.option(
    "--version",
    is_flag=True,
    callback=_print_versions,
    expose_value=False,
    is_eager=True,
    help="Show the library and dependency versions, then exit.",
)


def _conditions(temperature, ph, ionic_strength, concentration, partial_pressure):
    """Build a Conditions from repeated KEY=VALUE options."""

    def parse(pairs, what):
        out = {}
        for item in pairs:
            if "=" not in item:
                raise click.BadParameter(f"{what} must be given as SPECIES=VALUE, got {item!r}")
            name, value = item.split("=", 1)
            try:
                out[name.strip()] = float(value)
            except ValueError as exc:
                raise click.BadParameter(f"could not read a number from {item!r}") from exc
        return out

    concentrations = parse(concentration, "--concentration")
    pressures = parse(partial_pressure, "--partial-pressure")
    model = "bdot" if ionic_strength > 0 else ("ideal" if concentrations or pressures else "unit")
    return Conditions(
        temperature_c=temperature,
        pH=ph,
        ionic_strength=ionic_strength,
        concentrations=concentrations,
        partial_pressures=pressures,
        activity_model=model,
    )


_common_options = [
    click.option(
        "--temperature",
        "-T",
        default=25.0,
        show_default=True,
        help="Temperature in Celsius (0.01 to 100).",
    ),
    click.option("--ph", default=7.0, show_default=True, help="pH."),
    click.option(
        "--ionic-strength",
        default=0.0,
        show_default=True,
        help="Ionic strength in molal; enables the B-dot activity model.",
    ),
    click.option(
        "--concentration",
        "-c",
        multiple=True,
        metavar="SPECIES=MOLAL",
        help="Set a species concentration; repeatable.",
    ),
    click.option(
        "--partial-pressure",
        "-p",
        multiple=True,
        metavar="SPECIES=BAR",
        help="Set a gas partial pressure; repeatable.",
    ),
]


def _with_common(function):
    for option in reversed(_common_options):
        function = option(function)
    return function


@click.group()
@_version_option
@click.option("--database", default=None, help="Path to an alternative pyGCC database file.")
def main(database):
    """Thermodynamics for microbial physiology and biogeochemistry."""
    if database:
        configure("pygcc", database=database)


@main.command()
@click.argument("reduced")
@click.argument("oxidized")
@click.option(
    "--key-element",
    default=None,
    help="Redox-active element, when it cannot be detected automatically.",
)
@_with_common
@_version_option
def couple(
    reduced, oxidized, key_element, temperature, ph, ionic_strength, concentration, partial_pressure
):
    """Balance one redox couple and report its potentials.

    REDUCED and OXIDIZED name the two forms, e.g. `mthermo couple HS- SO4-2`.
    """
    conditions = _conditions(temperature, ph, ionic_strength, concentration, partial_pressure)
    from . import get_backend

    pair = Couple.make(reduced, oxidized, key_element)
    result = HalfReactionResult(
        couple=pair,
        half=pair.half_reaction(),
        conditions=conditions,
        backend=get_backend(),
    )
    oxidation_state, reduced_state = result.oxidation_states()
    click.echo(f"reduction   {result.half.format('reduction')}")
    click.echo(f"oxidation   {result.half.format('oxidation')}")
    click.echo(f"n electrons {result.n_electrons}")
    click.echo(f"{result.half.key_element} oxidation state {oxidation_state} -> {reduced_state}")
    click.echo(f"E0          {result.E_standard.to('V').magnitude:+.3f} V")
    click.echo(f"E0'         {result.E_standard_prime.to('V').magnitude:+.3f} V")


@main.command()
@click.option(
    "--donor",
    "-d",
    nargs=2,
    required=True,
    metavar="REDUCED OXIDIZED",
    help="Electron donor couple.",
)
@click.option(
    "--acceptor",
    "-a",
    nargs=2,
    required=True,
    metavar="REDUCED OXIDIZED",
    help="Electron acceptor couple.",
)
@click.option(
    "--normalize-to",
    default=None,
    help="electron_pair (default), electron, donor, acceptor, integer, or a species name.",
)
@click.option(
    "--electrons",
    "-n",
    default=2,
    show_default=True,
    help="Electrons transferred, unless --normalize-to overrides it.",
)
@click.option("--show-work", is_flag=True, help="Print the full derivation.")
@_with_common
@_version_option
def reaction(
    donor,
    acceptor,
    normalize_to,
    electrons,
    show_work,
    temperature,
    ph,
    ionic_strength,
    concentration,
    partial_pressure,
):
    """Compute a whole reaction from a donor and an acceptor couple.

    Example: mthermo reaction -d H2(g) H+ -a HS- SO4-2 --normalize-to integer
    """
    conditions = _conditions(temperature, ph, ionic_strength, concentration, partial_pressure)
    result = Reaction.from_couples(
        donor=donor,
        acceptor=acceptor,
        conditions=conditions,
        n_electrons=electrons,
        normalize_to=normalize_to,
    )
    click.echo(result.summary())
    if show_work:
        click.echo()
        click.echo(result.show_work().to_text())


@main.command()
@click.option("--donor", "-d", nargs=2, required=True, metavar="REDUCED OXIDIZED")
@click.option("--acceptor", "-a", nargs=2, required=True, metavar="REDUCED OXIDIZED")
@click.option(
    "--kind",
    type=click.Choice(["half", "tower", "explorer"]),
    default="half",
    show_default=True,
    help="Which figure to draw.",
)
@click.option(
    "--out", "-o", required=True, metavar="PATH", help="Output path without an extension."
)
@click.option("--normalize-to", default="integer", show_default=True)
@click.option("--electrons", "-n", default=2, show_default=True)
@_with_common
@_version_option
def figure(
    donor,
    acceptor,
    kind,
    out,
    normalize_to,
    electrons,
    temperature,
    ph,
    ionic_strength,
    concentration,
    partial_pressure,
):
    """Render a figure to SVG and PNG, or to standalone HTML for the explorer."""
    conditions = _conditions(temperature, ph, ionic_strength, concentration, partial_pressure)
    result = Reaction.from_couples(
        donor=donor,
        acceptor=acceptor,
        conditions=conditions,
        n_electrons=electrons,
        normalize_to=normalize_to,
    )
    if kind == "explorer":
        from .figures.explorer import plot_energy_explorer

        plot_energy_explorer(result, save_html=out)
        click.echo(f"wrote {out}.html")
        return

    import matplotlib

    matplotlib.use("Agg")
    if kind == "half":
        from .figures.halfreaction import plot_half_reactions

        plot_half_reactions(result, save=out)
    else:
        from .figures.tower import plot_redox_tower

        plot_redox_tower(result, save=out)
    click.echo(f"wrote {out}.svg and {out}.png")


@main.command()
@_with_common
@_version_option
def tower(temperature, ph, ionic_strength, concentration, partial_pressure):
    """Print the reference redox couples, ordered by potential."""
    from .tower import tower_table

    conditions = _conditions(temperature, ph, ionic_strength, concentration, partial_pressure)
    click.echo(tower_table(conditions).to_string(index=False))


@main.command()
@click.argument("family", required=False, default="")
@_with_common
@_version_option
def speciation(family, temperature, ph, ionic_strength, concentration, partial_pressure):
    """Show acid-base speciation at a given pH and temperature.

    With no FAMILY, lists every family and its pKa ladder. With one, prints the
    distribution across its members, e.g. `mthermo speciation sulfide --ph 7`.
    """
    from .speciation import default_families, pKa_ladder, speciation_table

    if not family:
        click.echo(f"pKa values at {temperature:g} C:")
        for entry in default_families():
            ladder = ", ".join(f"{v:.2f}" for v in pKa_ladder(entry, temperature_c=temperature))
            members = " / ".join(m.backend for m in entry.members)
            click.echo(f"  {entry.name:<11} {ladder:<22} {members}")
        return

    table = speciation_table(family, ph, temperature_c=temperature)
    click.echo(f"{family} at pH {ph:g}, {temperature:g} C:")
    click.echo(table.to_string(index=False))


@main.command()
@click.argument("pattern", required=False, default="")
@_version_option
def species(pattern):
    """List known species, optionally filtered by a substring."""
    registry = default_registry()
    needle = pattern.lower()
    rows = [
        s
        for s in registry.all_species()
        if not needle or needle in s.backend.lower() or any(needle in a.lower() for a in s.aliases)
    ]
    if not rows:
        click.echo(f"no species matching {pattern!r}")
        return
    click.echo(f"{'backend name':<18} {'formula':<10} {'phase':<6} aliases")
    for entry in sorted(rows, key=lambda s: s.backend):
        click.echo(
            f"{entry.backend:<18} {entry.formula:<10} {entry.phase:<6} "
            f"{', '.join(entry.aliases[:4])}"
        )


if __name__ == "__main__":
    main()
