"""Streamlit web frontend for microbial_thermo.

Three teaching tools over the existing library -- compound formation-energy
lookup, half-reaction potentials, and a free-text reaction balancer -- with
one shared conditions panel. See webfrontend.md for the full specification
this implements.

Run with:

    streamlit run app.py
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

import microbial_thermo as mt

MAX_HISTORY = 10
ACTIVITY_MODELS = ("ideal", "bdot", "unit")


# --- backend and cached lookups ---------------------------------------------
# Conditions objects are not clean cache keys (a mutable dataclass with dict
# fields), so every cached function below takes scalar arguments and a
# conditions "key" tuple, and builds Conditions inside. See webfrontend.md §5.


@st.cache_resource
def backend():
    return mt.get_backend()


def conditions_key(conditions: mt.Conditions) -> tuple:
    return (
        conditions.temperature_c,
        conditions.pH,
        conditions.pressure_bar,
        conditions.ionic_strength,
        tuple(sorted(conditions.concentrations.items())),
        tuple(sorted(conditions.total_concentrations.items())),
        tuple(sorted(conditions.partial_pressures.items())),
        conditions.activity_model,
        conditions.water_activity,
    )


def conditions_from_key(key: tuple) -> mt.Conditions:
    (
        temperature_c,
        pH,
        pressure_bar,
        ionic_strength,
        concentrations,
        total_concentrations,
        partial_pressures,
        activity_model,
        water_activity,
    ) = key
    return mt.Conditions(
        temperature_c=temperature_c,
        pH=pH,
        pressure_bar=pressure_bar,
        ionic_strength=ionic_strength,
        concentrations=dict(concentrations),
        total_concentrations=dict(total_concentrations),
        partial_pressures=dict(partial_pressures),
        activity_model=activity_model,
        water_activity=water_activity,
    )


@st.cache_data(show_spinner=False)
def lookup_species(name: str, phase: str | None = None):
    return mt.default_registry().resolve(name, phase=phase)


@st.cache_data(show_spinner=False)
def lookup_species_record(species_backend: str):
    return backend().resolve(species_backend)


@st.cache_data(show_spinner=False)
def lookup_dgf(species_backend: str, temperature_c: float, pressure_bar: float | None):
    return backend().delta_Gf(species_backend, temperature_c, pressure_bar)


@st.cache_data(show_spinner=False)
def lookup_half_reaction(reduced: tuple[str, ...], oxidized: tuple[str, ...], key: tuple):
    conditions = conditions_from_key(key)
    r = reduced[0] if len(reduced) == 1 else list(reduced)
    o = oxidized[0] if len(oxidized) == 1 else list(oxidized)
    return mt.half_reaction(r, o, conditions=conditions)


@st.cache_data(show_spinner=False)
def lookup_reaction(equation: str, key: tuple):
    conditions = conditions_from_key(key)
    return mt.Reaction.from_equation(equation, conditions=conditions)


@st.cache_data(show_spinner=False)
def lookup_reaction_from_couples(
    donor_reduced: str,
    donor_oxidized: str,
    acceptor_reduced: str,
    acceptor_oxidized: str,
    key: tuple,
):
    conditions = conditions_from_key(key)
    return mt.Reaction.from_couples(
        donor=(donor_reduced, donor_oxidized),
        acceptor=(acceptor_reduced, acceptor_oxidized),
        conditions=conditions,
    )


# --- small helpers -----------------------------------------------------------


def _parse_multi(text: str) -> list[str]:
    return [s.strip() for s in text.split(",") if s.strip()]


def _normalise(name: str) -> str:
    return str(name).strip().lower().replace(" ", "_")


def _rows_to_dict(df: pd.DataFrame, key_col: str, val_col: str) -> dict:
    out: dict = {}
    if df is None or df.empty:
        return out
    for _, row in df.iterrows():
        name = row.get(key_col)
        if name is None or (isinstance(name, float) and pd.isna(name)):
            continue
        name = str(name).strip()
        value = row.get(val_col)
        if not name or value is None or (isinstance(value, float) and pd.isna(value)):
            continue
        try:
            out[name] = float(value)
        except (TypeError, ValueError):
            continue
    return out


def _select_suggestion(target_key: str, value: str) -> None:
    st.session_state[target_key] = value


def _suggestion_buttons(suggestions: list[str], key_prefix: str, target_key: str) -> None:
    if not suggestions:
        return
    st.write("Did you mean:")
    cols = st.columns(len(suggestions))
    for col, suggestion in zip(cols, suggestions, strict=True):
        # Set session_state via on_click (which runs before the script body
        # re-executes), not after the fact -- a widget's key cannot be
        # written once that widget has been instantiated this run.
        col.button(
            suggestion,
            key=f"{key_prefix}_{suggestion}",
            on_click=_select_suggestion,
            args=(target_key, suggestion),
        )


def _unverified_footnote(provenance) -> None:
    unverified = provenance.unverified
    if not unverified:
        return
    names = ", ".join(s.name for s in unverified)
    st.markdown(f":red[**UNVERIFIED** -- rests on hand-entered value(s) for {names}.]")


# --- shared conditions panel (sidebar, §6.4) --------------------------------


def render_conditions_sidebar() -> mt.Conditions:
    st.sidebar.header("Conditions")
    st.sidebar.caption("Feeds tools 2 and 3. Tool 1 only uses temperature and pressure.")

    if st.sidebar.button("Reset to standard state"):
        for k in list(st.session_state.keys()):
            if k.startswith("cond_"):
                del st.session_state[k]
        st.rerun()

    defaults = mt.Conditions()

    temperature_c = st.sidebar.number_input(
        "Temperature (°C)",
        min_value=0.01,
        max_value=100.0,
        value=defaults.temperature_c,
        step=1.0,
        key="cond_temperature",
        help="pyGCC's supported window is 0.01-100 °C at near-surface pressure.",
    )
    pH = st.sidebar.number_input("pH", value=defaults.pH, step=0.1, key="cond_pH")

    activity_model = st.sidebar.selectbox(
        "Activity model",
        ACTIVITY_MODELS,
        index=ACTIVITY_MODELS.index(defaults.activity_model),
        key="cond_activity_model",
    )
    ionic_strength = defaults.ionic_strength
    if activity_model == "bdot":
        ionic_strength = st.sidebar.number_input(
            "Ionic strength (mol/kg)",
            min_value=0.0001,
            value=0.25,
            step=0.05,
            key="cond_ionic_strength",
            help="Required by the 'bdot' activity model.",
        )

    set_pressure = st.sidebar.checkbox(
        "Set pressure explicitly", key="cond_set_pressure", value=False
    )
    pressure_bar = None
    if set_pressure:
        pressure_bar = st.sidebar.number_input(
            "Pressure (bar)", value=1.0, step=0.1, key="cond_pressure"
        )

    st.sidebar.caption("Concentrations (species → molality)")
    conc_rows = st.sidebar.data_editor(
        pd.DataFrame({"species": pd.Series(dtype="str"), "molality": pd.Series(dtype="float")}),
        num_rows="dynamic",
        key="cond_concentrations_editor",
        hide_index=True,
        use_container_width=True,
    )
    concentrations = _rows_to_dict(conc_rows, "species", "molality")

    st.sidebar.caption(
        "Total concentrations (family → total molality): sulfide, carbonate, "
        "ammonia, phosphate, acetate, lactate, formate, propanoate, "
        "butanoate, sulfite, nitrite, sulfate"
    )
    total_rows = st.sidebar.data_editor(
        pd.DataFrame(
            {"family": pd.Series(dtype="str"), "total_molality": pd.Series(dtype="float")}
        ),
        num_rows="dynamic",
        key="cond_total_editor",
        hide_index=True,
        use_container_width=True,
    )
    total_concentrations = _rows_to_dict(total_rows, "family", "total_molality")

    st.sidebar.caption("Partial pressures (gas species → bar)")
    pp_rows = st.sidebar.data_editor(
        pd.DataFrame({"species": pd.Series(dtype="str"), "bar": pd.Series(dtype="float")}),
        num_rows="dynamic",
        key="cond_pp_editor",
        hide_index=True,
        use_container_width=True,
    )
    partial_pressures = _rows_to_dict(pp_rows, "species", "bar")

    st.sidebar.caption(
        "Ionic strength does not reconcile this with an Alberty-convention "
        "table; it moves ΔG at stated concentrations but not ΔG°′, which is "
        "defined at unit activity."
    )

    try:
        conditions = mt.Conditions(
            temperature_c=temperature_c,
            pH=pH,
            pressure_bar=pressure_bar,
            ionic_strength=ionic_strength,
            concentrations=concentrations,
            total_concentrations=total_concentrations,
            partial_pressures=partial_pressures,
            activity_model=activity_model,
        )
    except ValueError as exc:
        st.sidebar.error(str(exc))
        conditions = mt.Conditions()

    with st.sidebar.expander("About / versions"):
        info = mt.versions()
        info["streamlit"] = st.__version__
        st.json(info)

    return conditions


# --- tool 1: compound dGf lookup (§6.1) -------------------------------------


def render_tool1() -> None:
    st.subheader("Compound formation energy (ΔG_f°) lookup")

    if "tool1_history" not in st.session_state:
        st.session_state.tool1_history = []

    col1, col2, col3 = st.columns([3, 1, 1])
    with col1:
        name = st.text_input(
            "Species name", key="tool1_name", placeholder="e.g. HS-, SO4-2, Acetate"
        )
    with col2:
        temperature_c = st.number_input(
            "Temperature (°C)",
            min_value=0.01,
            max_value=100.0,
            value=25.0,
            step=1.0,
            key="tool1_temp",
        )
    with col3:
        use_pressure = st.checkbox("Set pressure", key="tool1_use_pressure")
        pressure_bar = None
        if use_pressure:
            pressure_bar = st.number_input(
                "Pressure (bar)", value=1.0, step=0.1, key="tool1_pressure"
            )

    if not name.strip():
        return

    try:
        species = lookup_species(name)
    except mt.SpeciesNotFoundError as exc:
        st.error(str(exc))
        _suggestion_buttons(exc.suggestions, "tool1_suggest", "tool1_name")
        return

    # Ambiguous names (H2, O2, N2, CO2, CH4, ...): offer a phase toggle,
    # defaulting to the dissolved form -- that is the form a cell sees.
    ambiguities = mt.default_registry().ambiguities()
    key = _normalise(name)
    if key in ambiguities:
        _, claimants, _ = ambiguities[key]
        phases = sorted({c.phase for c in claimants})
        if len(phases) > 1:
            default_phase = species.phase if species.phase in phases else phases[0]
            phase_choice = st.radio(
                "Phase",
                phases,
                index=phases.index(default_phase),
                horizontal=True,
                key=f"tool1_phase_{key}",
                help="Defaults to the dissolved form -- that is the form a cell sees.",
            )
            try:
                species = lookup_species(name, phase=phase_choice)
            except mt.SpeciesNotFoundError as exc:
                st.error(str(exc))
                return

    try:
        dgf = lookup_dgf(species.backend, temperature_c, pressure_bar)
    except mt.OutOfRangeError as exc:
        st.error(str(exc))
        return
    except mt.MissingDataError as exc:
        st.error(str(exc))
        return

    dgf_kj = dgf.to("kJ/mol").magnitude
    st.metric(f"ΔG_f° [{species.label}]", f"{dgf_kj:.2f} kJ/mol")
    st.caption(
        f"Backend name `{species.backend}` · formula {species.formula} · phase {species.phase}"
    )
    with st.expander("Raw value (backend's native unit)"):
        st.write(f"{dgf.to('cal/mol').magnitude:.1f} cal/mol")

    record = lookup_species_record(species.backend)
    if not record.verified:
        st.markdown(f":red[**UNVERIFIED** -- hand-entered value. Source: {record.source}]")

    st.session_state.tool1_history.insert(
        0,
        {
            "species": species.label,
            "backend name": species.backend,
            "phase": species.phase,
            "T (°C)": temperature_c,
            "ΔG_f° (kJ/mol)": round(dgf_kj, 2),
        },
    )
    st.session_state.tool1_history = st.session_state.tool1_history[:MAX_HISTORY]

    if st.session_state.tool1_history:
        st.divider()
        st.caption("Recent lookups this session")
        st.dataframe(
            pd.DataFrame(st.session_state.tool1_history), hide_index=True, use_container_width=True
        )


# --- tool 2: half-reaction E°/E°′ lookup (§6.2) -----------------------------


def _render_half_reaction(result: mt.HalfReactionResult) -> None:
    st.code(str(result), language=None)
    oxidation_state, reduced_state = result.oxidation_states()
    st.caption(
        f"n = {result.n_electrons} electrons · {result.half.key_element} oxidation "
        f"state {oxidation_state} → {reduced_state}"
    )
    c1, c2, c3 = st.columns(3)
    c1.metric("E°′ (headline)", f"{result.E_standard_prime.to('V').magnitude:+.3f} V")
    c2.metric("E°", f"{result.E_standard.to('V').magnitude:+.3f} V")
    c3.metric("E (conditions)", f"{result.E.to('V').magnitude:+.3f} V")


def render_tool2(conditions: mt.Conditions) -> None:
    st.subheader("Half-reaction E°/E°′/E lookup")

    col1, col2 = st.columns(2)
    with col1:
        reduced_input = st.text_input(
            "Reduced species (comma-separated for multiple)",
            key="tool2_reduced",
            placeholder="e.g. HS-",
        )
    with col2:
        oxidized_input = st.text_input(
            "Oxidized species (comma-separated for multiple)",
            key="tool2_oxidized",
            placeholder="e.g. SO4-2",
        )

    if not (reduced_input.strip() and oxidized_input.strip()):
        return

    reduced_ok, oxidized_ok = True, True
    for column, text, field_ok_key in (
        (col1, reduced_input, "reduced"),
        (col2, oxidized_input, "oxidized"),
    ):
        for one_name in _parse_multi(text):
            try:
                lookup_species(one_name)
            except mt.SpeciesNotFoundError as exc:
                if field_ok_key == "reduced":
                    reduced_ok = False
                else:
                    oxidized_ok = False
                with column:
                    st.error(str(exc))
                    _suggestion_buttons(
                        exc.suggestions, f"tool2_suggest_{field_ok_key}", f"tool2_{field_ok_key}"
                    )

    if not (reduced_ok and oxidized_ok):
        return

    try:
        result = lookup_half_reaction(
            tuple(_parse_multi(reduced_input)),
            tuple(_parse_multi(oxidized_input)),
            conditions_key(conditions),
        )
    except mt.MissingDataError as exc:
        st.error(str(exc))
    except mt.OutOfRangeError as exc:
        st.error(str(exc))
    except mt.ThermodynamicConsistencyError as exc:
        st.error(f"Internal consistency check failed -- this is a bug, please report it. {exc}")
    else:
        _render_half_reaction(result)


# --- tool 3: reaction balancer, reveal dG°′ (§6.3) --------------------------


def _render_ambiguity_basis(basis) -> None:
    if not (isinstance(basis, tuple) and len(basis) == 2):
        return
    donors, acceptors = basis
    rows = []
    for reduced, oxidized, element in donors:
        rows.append(
            {
                "role": "possible donor",
                "reduced": reduced.backend,
                "oxidized": oxidized.backend,
                "element": element,
            }
        )
    for reduced, oxidized, element in acceptors:
        rows.append(
            {
                "role": "possible acceptor",
                "reduced": reduced.backend,
                "oxidized": oxidized.backend,
                "element": element,
            }
        )
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)


def _render_reaction(reaction: mt.Reaction) -> None:
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Donor (oxidation)**")
        st.code(reaction.donor_half.format("oxidation"), language=None)
        st.metric("E°′", f"{reaction.donor_half.E_standard_prime.to('V').magnitude:+.3f} V")
    with c2:
        st.markdown("**Acceptor (reduction)**")
        st.code(reaction.acceptor_half.format("reduction"), language=None)
        st.metric("E°′", f"{reaction.acceptor_half.E_standard_prime.to('V').magnitude:+.3f} V")

    st.caption(f"n = {reaction.n_electrons} electrons")
    st.code(str(reaction), language=None)

    _unverified_footnote(reaction.provenance())

    show_dg = st.toggle("Reveal ΔG°′", key="tool3_reveal")
    if show_dg:
        d1, d2 = st.columns(2)
        d1.metric("ΔG°′", f"{reaction.delta_G_standard_prime.to('kJ/mol').magnitude:+.1f} kJ/mol")
        d2.metric("ΔG (in-situ)", f"{reaction.delta_G.to('kJ/mol').magnitude:+.1f} kJ/mol")


def _render_couple_picker_fallback(conditions: mt.Conditions) -> None:
    with st.expander("Or specify the couples directly"):
        cc1, cc2 = st.columns(2)
        with cc1:
            st.caption("Donor couple (reduced, oxidized)")
            donor_reduced = st.text_input("Donor reduced", key="tool3_donor_reduced")
            donor_oxidized = st.text_input("Donor oxidized", key="tool3_donor_oxidized")
        with cc2:
            st.caption("Acceptor couple (reduced, oxidized)")
            acceptor_reduced = st.text_input("Acceptor reduced", key="tool3_acceptor_reduced")
            acceptor_oxidized = st.text_input("Acceptor oxidized", key="tool3_acceptor_oxidized")

        if not (donor_reduced and donor_oxidized and acceptor_reduced and acceptor_oxidized):
            return
        if not st.button("Balance from couples", key="tool3_balance_couples"):
            return
        try:
            reaction = lookup_reaction_from_couples(
                donor_reduced,
                donor_oxidized,
                acceptor_reduced,
                acceptor_oxidized,
                conditions_key(conditions),
            )
        except mt.SpeciesNotFoundError as exc:
            st.error(str(exc))
            _suggestion_buttons(exc.suggestions, "tool3_couple_suggest", "tool3_donor_reduced")
        except mt.BalancingError as exc:
            st.error(str(exc))
        except mt.ThermodynamicConsistencyError as exc:
            st.error(f"Internal consistency check failed -- this is a bug, please report it. {exc}")
        else:
            _render_reaction(reaction)


def render_tool3(conditions: mt.Conditions) -> None:
    st.subheader("Reaction balancer (to an electron pair)")

    equation = st.text_input(
        "Equation", key="tool3_equation", placeholder="e.g. NO3- + H2 -> NH2OH"
    )
    if not equation.strip():
        return

    try:
        reaction = lookup_reaction(equation, conditions_key(conditions))
    except mt.AmbiguousReactionError as exc:
        st.error(str(exc))
        _render_ambiguity_basis(exc.basis)
        st.caption(
            "This is the one cost of free-text entry: the equation needs a "
            "constraint the library can't guess. Name the couples directly instead."
        )
        _render_couple_picker_fallback(conditions)
    except mt.SpeciesNotFoundError as exc:
        st.error(str(exc))
        _suggestion_buttons(exc.suggestions, "tool3_suggest", "tool3_equation")
    except mt.BalancingError as exc:
        st.error(str(exc))
    except mt.OutOfRangeError as exc:
        st.error(str(exc))
    except mt.ThermodynamicConsistencyError as exc:
        st.error(f"Internal consistency check failed -- this is a bug, please report it. {exc}")
    else:
        _render_reaction(reaction)


# --- page --------------------------------------------------------------------


def main() -> None:
    st.set_page_config(page_title="microbial_thermo", layout="wide")
    st.title("microbial_thermo")
    st.caption(
        "Thermodynamics for microbial physiology and biogeochemistry -- "
        "a teaching tool for CSU's Topics in Microbial Physiology and Biogeochemistry course."
    )

    conditions = render_conditions_sidebar()

    tab1, tab2, tab3 = st.tabs(
        ["1. Compound ΔG_f° lookup", "2. Half-reaction potentials", "3. Reaction balancer"]
    )
    with tab1:
        render_tool1()
    with tab2:
        render_tool2(conditions)
    with tab3:
        render_tool3(conditions)


if __name__ == "__main__":
    main()
