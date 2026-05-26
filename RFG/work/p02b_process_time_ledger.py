# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.
# Active coefficient scheme: Y-scheme; c_Y means c_Y^(Y), not X-scheme c_X.

"""
p02b: process-time ledger.

ეს ფაილი არის ცალკე ledger process-time იდეებისთვის. ის არ არის primary
FLRW/CMB metric branch და არ ამტკიცებს დამოუკიდებელ metric mode-ს.
"""

import sympy as sp


DEFAULT_PRESENT_AGE_GYR = 13.8


def get_process_time_relation():
    """
    ისტორიული process-time ხაზობრივი ესკიზის დაკონსოლიდებული ფორმა.

    აქტიური normalization ემთხვევა p10-ის ბი-კონფორმულ სკალირებას:
    d tau / dt = exp(phi/2). ეს ledger არ შედის primary FLRW/CMB metric
    branch-ში.
    """
    t = sp.Symbol("t", real=True)
    tau = sp.Function("tau")(t)
    phi = sp.Symbol("phi", real=True)

    relation = sp.Eq(sp.diff(tau, t), sp.exp(phi / 2))
    weak_field = sp.series(sp.exp(phi / 2), phi, 0, 2).removeO()

    return relation, weak_field, tau, t, phi


def process_time_scaling_ledger():
    """
    წნევის-ისტორიის process-time სკალირება.

    C(z) ნორმირებულია C(0)=1-ით და გამოიყენება მხოლოდ წინასწარ მონიშნულ
    intrinsic process-rate სიდიდეებზე.
    """
    z = sp.Symbol("z", real=True, nonnegative=True)
    phi_bg = sp.Function("phi_bg")
    C = sp.exp((phi_bg(z) - phi_bg(0)) / 2)

    return {
        "C_of_z": sp.Eq(sp.Function("C")(z), C),
        "process_clock": "d tau_proc / dt = C(z)",
        "rod_scale": "d ell(z) / d ell(0) = C(z)^(-1)",
        "mass_scale": "m_eff(z) / m_eff(0) = C(z)",
        "tail_power_perturbative": "P_tail(z) / P_tail(0) = C(z)^4",
        "normalization": "C(0)=1 by present-epoch calibration",
        "status": "bookkeeping only; not an independent metric mode",
    }


def process_time_source_requirements():
    """
    რა უნდა დაიხუროს, სანამ C(z) ფიზიკურ მოდელად ჩაითვლება.
    """
    z = sp.Symbol("z", real=True, nonnegative=True)
    phi_bg = sp.Function("phi_bg")

    return {
        "C_of_z_source": "POSTULATE_CANDIDATE_NOT_DERIVED",
        "active_definition": sp.Eq(
            sp.Function("C")(z),
            sp.exp((phi_bg(z) - phi_bg(0)) / 2),
        ),
        "candidate_closure": (
            "self-similar process-time postulate: every epoch measures the "
            "infinite internal past as the same finite age T0 in its own units"
        ),
        "postulate_to_field_bridge": (
            "self_similar_Cz_phi_bridge algebraically gives "
            "C(z)=T0/t_age(z) and phi_bg(z)-phi_bg(0)=2*log(C(z))"
        ),
        "required_before_physics_claim": [
            "derive or observationally fit phi_bg(z), not only define it algebraically",
            "state which physical process is parameterized by tau_proc",
            "show C(0)=1 calibration",
            "provide BBN/CMB/structure-growth bounds on C(z)",
            "write an energy/entropy budget if C(z) changes tail power or rates",
        ],
        "allowed_interim_status": (
            "postulate-level diagnostic calculation only; no solved JWST, H0, "
            "dark-energy, or BBN/CMB compatibility claim"
        ),
    }


def process_time_self_similar_postulate():
    """
    Minimal self-similar process-time postulate.

    This is not derived from p02 FLRW algebra. It is a narrow calculational
    ansatz: the CMB-comoving present calibration T0 is finite even if the
    internal past process coordinate is unbounded.
    """
    tau = sp.Symbol("tau", real=True)
    s = sp.Symbol("s", real=True)
    T0 = sp.Symbol("T0", positive=True, real=True)
    Delta_today = sp.Symbol("Delta_today", nonnegative=True, real=True)
    C = sp.Function("C")
    t_today = sp.Function("t_today")
    lambda_future = sp.Symbol("lambda_future", positive=True, real=True)

    C_tau = sp.exp(-tau / T0)
    today_age_at_tau = T0 * sp.exp(tau / T0)
    finite_present_age_integral = sp.Integral(sp.exp(s / T0), (s, -sp.oo, 0))
    finite_present_age_value = sp.integrate(sp.exp(s / T0), (s, -sp.oo, 0))
    future_elapsed_today_units = T0 * (sp.exp(tau / T0) - 1)
    future_unit_ratio_from_delta = 1 + Delta_today / T0
    future_age_in_future_units = (T0 + Delta_today) / future_unit_ratio_from_delta

    return {
        "status": "POSTULATE_CANDIDATE_NOT_DERIVED",
        "postulate": (
            "The internal process past may be infinite, while each epoch's own "
            "clock unit reads the full past as the same finite age T0."
        ),
        "T0_role": (
            "symbolic CMB-comoving present-epoch calibrated age; 13.8 Gyr is "
            "only today's numerical fit example, not a new fundamental constant"
        ),
        "C_tau": sp.Eq(C(tau), C_tau),
        "process_clock": "d tau_proc / dt_today = C(tau)",
        "today_unit_age_at_internal_tau": sp.Eq(t_today(tau), today_age_at_tau),
        "present_age_from_infinite_internal_past": sp.Eq(
            finite_present_age_integral,
            T0,
            evaluate=False,
        ),
        "present_age_integral_value": finite_present_age_value,
        "future_elapsed_today_units": sp.Eq(
            Delta_today,
            future_elapsed_today_units,
        ),
        "future_unit_ratio": sp.Eq(lambda_future, future_unit_ratio_from_delta),
        "future_age_measured_in_future_units": sp.Eq(
            future_age_in_future_units,
            T0,
            evaluate=False,
        ),
        "future_age_simplified": sp.simplify(future_age_in_future_units),
        "age_invariance_status": "DERIVED_WITHIN_THIS_POSTULATE_NOT_SEPARATE_POSTULATE",
        "equivalent_age_rule": "C(t_age) = T0 / t_age in present-calibrated age units",
        "JWST_use": (
            "can enlarge intrinsic formation/process time without replacing "
            "metric redshift or observed time dilation"
        ),
        "non_uniqueness_warning": (
            "infinite internal past plus present age T0 does not uniquely fix "
            "this exponential curve; self-similarity is the extra postulate"
        ),
        "forbidden_reading": "not tired-light; not an extra redshift law",
    }


def cosmic_age_invariance_result():
    """
    The age-invariance result derived from the self-similar postulate.
    """
    T0 = sp.Symbol("T0", positive=True, real=True)
    Delta_today = sp.Symbol("Delta_today", nonnegative=True, real=True)
    lambda_epoch_symbol = sp.Symbol("lambda_epoch", positive=True, real=True)
    lambda_epoch_expr = sp.Mul(T0 + Delta_today, 1 / T0, evaluate=False)
    age_before_substitution = sp.Mul(
        T0 + Delta_today,
        sp.Pow(lambda_epoch_symbol, -1, evaluate=False),
        evaluate=False,
    )
    age_after_substitution = sp.Mul(
        T0 + Delta_today,
        sp.Pow(lambda_epoch_expr, -1, evaluate=False),
        evaluate=False,
    )

    return {
        "status": "DERIVED_RESULT_NOT_NEW_POSTULATE",
        "input_postulate": "process_time_self_similar_postulate",
        "epoch_unit_ratio": sp.Eq(
            lambda_epoch_symbol,
            lambda_epoch_expr,
            evaluate=False,
        ),
        "age_in_epoch_units_before_substitution": sp.Eq(
            sp.Symbol("T_epoch_units", positive=True, real=True),
            age_before_substitution,
            evaluate=False,
        ),
        "age_in_epoch_units": sp.Eq(
            age_after_substitution,
            T0,
            evaluate=False,
        ),
        "age_in_epoch_units_simplified": sp.simplify(age_after_substitution),
        "meaning": (
            "future or past epochs read the full cosmological past as T0 in "
            "their own clock calibration"
        ),
        "T0_not_numeric_constant": (
            "leave T0 symbolic in theory; insert 13.8 Gyr only as the current "
            "CMB-comoving observational calibration example"
        ),
    }


def self_similar_Cz_phi_bridge():
    """
    Bridge the self-similar postulate to the active C(z)/phi_bg(z) notation.

    After the postulate is accepted, this part is algebraic: if the standard
    CMB-comoving age at redshift z is t_age(z), then C(z)=T0/t_age(z). The
    phi_bg(z) expression is the field-space bookkeeping needed to match the
    existing C(z)=exp((phi_bg(z)-phi_bg(0))/2) convention.
    """
    z = sp.Symbol("z", real=True, nonnegative=True)
    T0 = sp.Symbol("T0", positive=True, real=True)
    t_age = sp.Function("t_age")
    C = sp.Function("C")
    phi_bg = sp.Function("phi_bg")
    H = sp.Function("H")
    a = sp.Symbol("a", positive=True)

    cmb_age_integral = sp.Integral(1 / (a * H(a)), (a, 0, 1 / (1 + z)))
    C_z = T0 / t_age(z)
    phi_difference = 2 * sp.log(C_z)
    phi_z_with_zero_today = phi_difference
    active_C_definition = sp.exp((phi_bg(z) - phi_bg(0)) / 2)
    bridge_check = sp.simplify(sp.exp(phi_difference / 2))

    return {
        "status": "ALGEBRAIC_BRIDGE_AFTER_SELF_SIMILAR_POSTULATE",
        "standard_age_input": sp.Eq(t_age(z), cmb_age_integral),
        "normalization": sp.Eq(t_age(0), T0),
        "C_z_self_similar": sp.Eq(C(z), C_z),
        "phi_bg_difference": sp.Eq(phi_bg(z) - phi_bg(0), phi_difference),
        "gauge_choice_phi_bg_today_zero": sp.Eq(phi_bg(z), phi_z_with_zero_today),
        "active_definition_check": sp.Eq(active_C_definition, bridge_check),
        "meaning": (
            "C(z) is now fixed by the CMB-comoving age history once the "
            "self-similar postulate is accepted"
        ),
        "remaining_physics_gap": (
            "phi_bg(z) is still a diagnostic background history until derived "
            "from an RFG field equation or fitted with observational bounds"
        ),
    }


def process_rate_vs_lapse_separation():
    """
    Keep the intrinsic process-rate factor separate from the metric lapse.
    """
    C_proc, alpha = sp.symbols("C_proc alpha", positive=True, real=True)
    dt_cmb = sp.Symbol("dt_CMB", positive=True, real=True)
    d_tau_proc = sp.Symbol("d_tau_proc", positive=True, real=True)
    d_tau_local = sp.Symbol("d_tau_local", positive=True, real=True)

    return {
        "status": "SEPARATION_REQUIRED",
        "process_rate_factor": sp.Eq(d_tau_proc, C_proc * dt_cmb),
        "local_metric_lapse": sp.Eq(d_tau_local, alpha * dt_cmb),
        "not_identical_by_default": sp.Ne(C_proc, alpha),
        "identification_allowed_only_if": (
            "a separate stress/lapse/substrate bridge proves that the same RFG "
            "field controls both the intrinsic process clock and the local lapse"
        ),
        "double_counting_rule": (
            "never apply both C_proc and alpha to the same observed time "
            "interval unless the variables and clock channel are explicitly distinct"
        ),
    }


def process_time_channel_law():
    """
    Channel law for applying C(z).

    The self-similar C(z) is dangerous if treated as a universal multiplier for
    all physics. This table restricts it to explicitly tagged intrinsic
    substrate/formation channels and blocks direct use in metric, photon,
    atomic, and early-universe standard-clock sectors until bounded.
    """
    return {
        "status": "CHANNEL_RESTRICTION_REQUIRED",
        "core_rule": (
            "C(z) may multiply exactly one explicitly tagged intrinsic "
            "process-rate channel; it is not a universal replacement for dt"
        ),
        "allowed_active_channels": [
            {
                "channel": "oscillon_tail_internal_relaxation",
                "rate_tag": "intrinsic_process_time_rate",
                "rule": "dX/dt_CMB = C(z) * dX/dtau_proc",
                "claim_level": "allowed diagnostic channel",
                "needed_closure": "derive the X evolution law in tau_proc",
            },
            {
                "channel": "early_galaxy_internal_formation_clock",
                "rate_tag": "formation_history_intrinsic_rate",
                "rule": "Delta_tau_proc controls internal maturation budget",
                "claim_level": "JWST budget support only",
                "needed_closure": "fit stellar mass, SFR efficiency, metallicity, dust, halo abundance",
            },
            {
                "channel": "resonance_internal_clock",
                "rate_tag": "resonance_internal_clock_rate",
                "rule": "local resonance phase can be parameterized by tau_proc",
                "claim_level": "conditional",
                "needed_closure": "show no conflict with observed spectral redshift/time dilation",
            },
        ],
        "blocked_standard_clock_channels": [
            {
                "channel": "photon_propagation_redshift",
                "reason": "would become tired-light or double-count metric redshift",
            },
            {
                "channel": "observed_supernova_or_quasar_time_dilation",
                "reason": "observed light-curve stretch remains metric (1+z)",
            },
            {
                "channel": "Einstein_Boltzmann_CMB_background",
                "reason": "would move p02b into the primary CMB metric branch",
            },
            {
                "channel": "BBN_nuclear_reaction_clock",
                "reason": "would change light-element yields unless separately bounded",
            },
            {
                "channel": "atomic_transition_frequency_clock",
                "reason": "would imply varying constants/laboratory-clock bounds",
            },
            {
                "channel": "stellar_evolution_microphysics",
                "reason": "would alter nuclear burning times unless a microphysical bridge is derived",
            },
        ],
        "conditional_channels_need_bound": [
            "tail_power_scaling",
            "entropy_production",
            "structure_growth_rate",
            "black-hole-near-lapse analogies",
        ],
        "single_factor_rule": (
            "If a rate is written per tau_proc, convert once with C(z). If a "
            "rate is already per FLRW/proper/observed time, do not add C(z)."
        ),
    }


def process_time_channel_gate(channel, rate_tag=None):
    """Gate a proposed channel use before applying C(z)."""
    law = process_time_channel_law()
    allowed = {row["channel"]: row for row in law["allowed_active_channels"]}
    blocked = {row["channel"]: row for row in law["blocked_standard_clock_channels"]}

    if channel in allowed:
        row = allowed[channel]
        tag_ok = rate_tag is None or rate_tag == row["rate_tag"]
        return {
            "channel": channel,
            "status": "PASS" if tag_ok else "BLOCKED",
            "reason": (
                "allowed tagged intrinsic process channel"
                if tag_ok
                else f"wrong rate_tag; expected {row['rate_tag']}"
            ),
            "rule": row["rule"],
            "needed_closure": row["needed_closure"],
        }

    if channel in blocked:
        row = blocked[channel]
        return {
            "channel": channel,
            "status": "BLOCKED",
            "reason": row["reason"],
            "rule": "do not apply C(z)",
        }

    return {
        "channel": channel,
        "status": "BLOCKED",
        "reason": "unknown channel; add explicit channel-law entry before applying C(z)",
        "rule": "no implicit C(z) factor",
    }


def clock_pressure_index_definition():
    """
    Operational index built from local cosmic-age reading.

    This is a measurement definition. Identifying it with RFG physical
    pressure requires a separate bridge to stress/lapse/substrate invariants.
    """
    T0 = sp.Symbol("T0", positive=True, real=True)
    T_local = sp.Symbol("T_local", positive=True, real=True)
    Pi_clock = sp.Symbol("Pi_clock", positive=True, real=True)

    return {
        "status": "OPERATIONAL_DEFINITION_NOT_NEW_TIME_POSTULATE",
        "definition": sp.Eq(Pi_clock, T0 / T_local),
        "earth_calibration": sp.Eq(
            sp.Symbol("Pi_clock_earth", positive=True, real=True),
            1,
            evaluate=False,
        ),
        "earth_calibration_substitution": sp.Eq(
            T_local,
            T0,
            evaluate=False,
        ),
        "earth_calibration_result": sp.Eq(
            Pi_clock,
            1,
            evaluate=False,
        ),
        "reading": (
            "smaller local age reading means larger clock-compression index"
        ),
        "derived_input": "uses T0 from cosmic_age_invariance_result",
        "not_yet_claimed": (
            "Pi_clock equals physical stress pressure only after an explicit "
            "bridge to RFG stress-energy, lapse, or substrate invariants"
        ),
    }


def clock_pressure_lapse_bridge():
    """
    Geometric bridge from local age reading to lapse.

    This strengthens Pi_clock as a lapse/compression index. It still does not
    prove equality with thermodynamic or stress-energy pressure.
    """
    T0 = sp.Symbol("T0", positive=True, real=True)
    T_local = sp.Symbol("T_local", positive=True, real=True)
    alpha = sp.Symbol("alpha", positive=True, real=True)
    Pi_clock = sp.Symbol("Pi_clock", positive=True, real=True)
    Phi = sp.Symbol("Phi", real=True)

    return {
        "status": "GEOMETRIC_LAPSE_BRIDGE_CONDITIONAL",
        "local_age_lapse_rule": sp.Eq(T_local, alpha * T0),
        "clock_pressure_index": sp.Eq(Pi_clock, 1 / alpha),
        "RFG_biconformal_lapse_if_gtt_exp_Phi": sp.Eq(alpha, sp.exp(Phi / 2)),
        "RFG_field_form_if_bridge_accepted": sp.Eq(Pi_clock, sp.exp(-Phi / 2)),
        "weak_field": sp.series(sp.exp(-Phi / 2), Phi, 0, 2).removeO(),
        "strengthened_claim": (
            "Pi_clock is no longer just a ratio; it is the inverse local lapse "
            "when local cosmic-age readings are controlled by a metric clock"
        ),
        "remaining_physics_gap": (
            "Pi_clock becomes physical RFG pressure only after a stress-energy "
            "or substrate-invariant derivation"
        ),
        "sign_warning": (
            "the sign of Phi must follow the active metric convention g_tt=exp(Phi); "
            "do not mix this with X/Y coefficient sign conventions"
        ),
    }


def schwarzschild_clock_pressure_reference():
    """
    GR reference example for local age compression outside a black hole.
    """
    r, r_s, T0 = sp.symbols("r r_s T0", positive=True, real=True)
    alpha = sp.sqrt(1 - r_s / r)
    T_local = T0 * alpha
    pi_clock_expr = sp.Pow(alpha, -1, evaluate=False)

    return {
        "status": "GR_ANALOGY_STATIC_EXTERIOR_ONLY",
        "validity": "r > r_s for a static exterior observer; not on the horizon",
        "lapse": sp.Eq(sp.Symbol("alpha", positive=True, real=True), alpha),
        "local_age_reading": sp.Eq(
            sp.Symbol("T_local", positive=True, real=True),
            T_local,
        ),
        "clock_pressure_index": sp.Eq(
            sp.Symbol("Pi_clock", positive=True, real=True),
            pi_clock_expr,
            evaluate=False,
        ),
        "interpretation": (
            "near the horizon alpha -> 0, so the local cosmic-age reading "
            "shrinks and Pi_clock grows"
        ),
        "guardrail": "do not treat this as a new RFG pressure law without the bridge",
    }


def process_time_age_conversion_example(
    present_age_gyr=DEFAULT_PRESENT_AGE_GYR,
    future_elapsed_today_gyr=2.0,
):
    """
    Convert a future measurement between present and future clock units.
    """
    T0 = sp.Float(present_age_gyr)
    Delta = sp.Float(future_elapsed_today_gyr)
    future_unit_ratio = sp.simplify((T0 + Delta) / T0)
    C_future = sp.simplify(1 / future_unit_ratio)
    internal_elapsed = sp.simplify(T0 * sp.log(future_unit_ratio))
    future_age_future_units = sp.simplify((T0 + Delta) / future_unit_ratio)

    return {
        "status": "NUMERICAL_CALIBRATION_EXAMPLE_ONLY",
        "present_age_gyr": float(T0),
        "future_elapsed_in_today_units_gyr": float(Delta),
        "future_age_in_today_units_gyr": float(T0 + Delta),
        "future_unit_ratio_to_today": float(future_unit_ratio),
        "C_future": float(C_future),
        "internal_process_elapsed_gyr": float(internal_elapsed),
        "future_age_measured_in_future_units_gyr": float(future_age_future_units),
        "interpretation": (
            "with the self-similar postulate, every epoch reads the full past "
            "as T0 in its own clock units"
        ),
    }


def _lcdm_high_z_age_approx_gyr(z, H0_km_s_mpc=67.4, omega_m=0.315):
    """
    Matter-era LCDM age approximation in Gyr.

    Valid only as a high-z diagnostic. It ignores radiation corrections and
    late dark-energy corrections, so it must not replace a Boltzmann/CLASS fit.
    """
    H0_inverse_gyr = 9.778 / (H0_km_s_mpc / 100.0)
    prefactor = 2.0 * H0_inverse_gyr / (3.0 * float(sp.sqrt(omega_m)))
    return prefactor * (1.0 + float(z)) ** (-1.5)


def jwst_process_time_diagnostic_table(
    redshifts=(7.3, 10.0, 12.0, 14.44, 20.0, 30.0),
    present_age_gyr=DEFAULT_PRESENT_AGE_GYR,
):
    """
    Diagnostic high-z process-time budget for early-galaxy discussions.
    """
    T0 = float(present_age_gyr)
    rows = []
    for z in redshifts:
        t_lcdm = _lcdm_high_z_age_approx_gyr(z)
        tau_internal = T0 * float(sp.log(t_lcdm / T0))
        local_process_factor = T0 / t_lcdm
        rows.append(
            {
                "z": float(z),
                "lcdm_age_gyr_approx": round(t_lcdm, 4),
                "internal_tau_from_present_gyr": round(tau_internal, 3),
                "local_process_factor_C": round(local_process_factor, 2),
                "status": "diagnostic_only_not_observational_fit",
            }
        )

    return {
        "postulate_used": "self-similar process-time",
        "T0_role": "symbolic; numeric value here is CMB-comoving present calibration",
        "present_age_calibration_gyr": T0,
        "age_approximation": (
            "high-z matter-era LCDM approximation with H0=67.4, Omega_m=0.315"
        ),
        "rows": rows,
        "claim_limit": (
            "supports a possible intrinsic formation-time budget; does not by "
            "itself solve JWST galaxy abundance, stellar mass, or metallicity"
        ),
    }


def jwst_process_time_interval_table(
    intervals=((30.0, 20.0), (20.0, 14.44), (15.0, 10.0), (12.0, 7.3)),
    present_age_gyr=DEFAULT_PRESENT_AGE_GYR,
):
    """
    Compare standard coordinate-time windows with process-time windows.
    """
    T0 = float(present_age_gyr)
    rows = []
    for z_start, z_end in intervals:
        t_start = _lcdm_high_z_age_approx_gyr(z_start)
        t_end = _lcdm_high_z_age_approx_gyr(z_end)
        standard_window = t_end - t_start
        process_window = T0 * float(sp.log(t_end / t_start))
        average_factor = process_window / standard_window
        rows.append(
            {
                "z_window": f"{z_start:g}->{z_end:g}",
                "lcdm_coordinate_window_gyr": round(standard_window, 4),
                "process_window_gyr": round(process_window, 3),
                "average_process_factor": round(average_factor, 1),
                "status": "formation_clock_diagnostic_only",
            }
        )

    return {
        "T0_role": "symbolic; numeric value here is CMB-comoving present calibration",
        "present_age_calibration_gyr": T0,
        "rows": rows,
        "interpretation": (
            "short early LCDM windows can correspond to much larger intrinsic "
            "process-time budgets under the postulate"
        ),
        "guardrail": "do not apply these factors to observed redshift time dilation",
    }


def jwst_formation_budget_test(
    required_formation_times_gyr=(0.3, 1.0, 3.0),
    intervals=((30.0, 20.0), (20.0, 14.44), (15.0, 10.0), (12.0, 7.3)),
    present_age_gyr=DEFAULT_PRESENT_AGE_GYR,
):
    """
    Turn the JWST diagnostic into a falsifiable formation-time budget test.
    """
    interval_rows = jwst_process_time_interval_table(
        intervals=intervals,
        present_age_gyr=present_age_gyr,
    )["rows"]

    tests = []
    for row in interval_rows:
        process_window = row["process_window_gyr"]
        coordinate_window = row["lcdm_coordinate_window_gyr"]
        for required_time in required_formation_times_gyr:
            tests.append(
                {
                    "z_window": row["z_window"],
                    "required_formation_time_gyr": float(required_time),
                    "lcdm_coordinate_window_gyr": coordinate_window,
                    "process_window_gyr": process_window,
                    "coordinate_budget_pass": coordinate_window >= required_time,
                    "process_budget_pass": process_window >= required_time,
                    "process_margin_gyr": round(process_window - required_time, 3),
                    "status": "budget_test_only_not_population_fit",
                }
            )

    return {
        "status": "FITTING_READY_BUDGET_TEST",
        "criterion": "required formation time <= Delta tau_proc(z_start,z_end)",
        "tests": tests,
        "still_required_for_JWST_claim": [
            "stellar mass and star-formation efficiency model",
            "metallicity and dust history",
            "halo abundance / selection function",
            "spectroscopic-redshift confirmed sample",
        ],
    }


def flrw_vs_process_time_separation():
    """
    წესი, რომელიც process-time ledger-ს primary FLRW/CMB branch-ისგან ყოფს.
    """
    return {
        "metric_FLRW_branch": "phi_0(t)=0 matter-clock cosmic-time gauge-ში",
        "process_diagnostic": "phi_bg(z)=coarse-grained pressure-history diagnostic",
        "not_a_second_metric_mode": True,
        "CMB_rule": (
            "phi_bg(z) არ შეიტანო Einstein-Boltzmann/CMB სექტორში, როგორც "
            "დამატებითი ჰომოგენური metric perturbation."
        ),
        "where_it_acts": (
            "C(z) გამოიყენე მხოლოდ მონიშნულ intrinsic matter/resonance "
            "process-rate კანონებში, formation-history ledger-ში და tail-emission "
            "bookkeeping-ში."
        ),
    }


def process_time_allowed_rate_tags():
    """Allowed/blocked use map for C(z)."""
    return {
        "allowed": {
            "intrinsic_process_time_rate",
            "formation_history_intrinsic_rate",
            "resonance_internal_clock_rate",
            "tail_emission_bookkeeping_rate",
        },
        "blocked": {
            "metric_FLRW_background",
            "einstein_boltzmann_CMB_background",
            "coordinate_cosmic_time_rate",
            "observed_redshift_time_dilation",
            "photon_propagation_redshift",
            "BBN_nuclear_reaction_clock",
            "atomic_transition_frequency_clock",
            "stellar_evolution_microphysics",
            "cluster_crossing_time",
            "merger_crossing_time",
        },
    }


def process_time_use_gate(rate_tag, *, enters_metric_branch=False, adds_extra_redshift=False):
    """
    Decide whether a proposed use of C(z) is allowed.

    This is a guardrail against double-counting and against silently turning
    process time into a second metric mode.
    """
    tags = process_time_allowed_rate_tags()

    reasons = []
    if enters_metric_branch:
        reasons.append("blocked: would enter the primary FLRW/CMB metric branch")
    if adds_extra_redshift:
        reasons.append("blocked: would double-count observed metric redshift/time dilation")
    if rate_tag in tags["blocked"]:
        reasons.append(f"blocked: {rate_tag} is not an intrinsic process-time channel")
    if rate_tag not in tags["allowed"] and rate_tag not in tags["blocked"]:
        reasons.append(f"blocked: unknown rate tag {rate_tag}; add an explicit ledger entry first")

    allowed = not reasons
    return {
        "rate_tag": rate_tag,
        "allowed": allowed,
        "status": "PASS" if allowed else "BLOCKED",
        "reasons": reasons or ["allowed only as a single C(z) factor in the tagged intrinsic channel"],
    }


def process_time_channel_audit_examples():
    """Representative channel-gate examples for the ledger output."""
    return {
        "allowed_galaxy_internal": process_time_channel_gate(
            "early_galaxy_internal_formation_clock",
            "formation_history_intrinsic_rate",
        ),
        "blocked_BBN": process_time_channel_gate("BBN_nuclear_reaction_clock"),
        "blocked_photon_redshift": process_time_channel_gate("photon_propagation_redshift"),
        "blocked_unknown": process_time_channel_gate("unregistered_process_channel"),
    }


def process_time_integrals():
    """Coordinate lookback და pressure-weighted process-time ინტეგრალები."""
    z, zp = sp.symbols("z z_prime", real=True, nonnegative=True)
    H = sp.Function("H")
    C = sp.Function("C")

    lookback = sp.Integral(1 / ((1 + zp) * H(zp)), (zp, 0, z))
    process = sp.Integral(C(zp) / ((1 + zp) * H(zp)), (zp, 0, z))
    enhancement = sp.Symbol("A_proc")

    return {
        "coordinate_lookback_time": sp.Eq(sp.Function("t_lb")(z), lookback),
        "cumulative_process_time": sp.Eq(sp.Function("tau_proc")(z), process),
        "enhancement_factor": sp.Eq(enhancement, process / lookback),
        "past_reading": "თუ C(z)>1 წარსულში, intrinsic process history იზრდება",
        "future_reading": "თუ C(t)->0 საკმარისად სწრაფად, process-time ბიუჯეტი შეიძლება სასრული გახდეს",
    }


def process_time_observational_anchor_requirements():
    """
    Observational anchors needed before C(z) can be treated as a fitted model.
    """
    return {
        "BBN": (
            "BLOCKED_AS_STANDARD_CLOCK: do not apply C(z) to nuclear reaction "
            "rates unless a separate bounded microphysical bridge is derived"
        ),
        "CMB": "BLOCKED_METRIC_BRANCH: C(z) must not enter as metric background mode",
        "structure_growth": "OPEN_BOUND: formation-rate use needs galaxy/halo benchmark",
        "local_epoch": "CALIBRATION: C(0)=1 by definition",
        "JWST_high_z": (
            "SUPPORT_CALC_ONLY: self-similar curve can enlarge intrinsic "
            "formation time, but needs galaxy abundance/mass/metallicity fit"
        ),
        "H0_tension": "NO_CLAIM: process-time ledger does not alter FLRW distance ladder by itself",
    }


def process_time_energy_budget_guardrail():
    """
    Energy/entropy accounting needed for tail-power or formation-rate claims.
    """
    return {
        "tail_power_scaling": "P_tail(z)/P_tail(0)=C(z)^4 is bookkeeping only",
        "required_budget": [
            "identify which reservoir supplies/removes changed tail power",
            "show no conflict with FLRW stress-energy bookkeeping in p02",
            "separate internal rate acceleration from observed luminosity/redshift effects",
            "state whether entropy production is changed or only reparameterized",
        ],
        "status": "OPEN_ENERGY_BUDGET",
    }


def rate_conversion_no_double_counting():
    """რომელ სიჩქარეში შედის C(z) და სად არ უნდა შევიდეს."""
    return [
        {
            "rate_tag": "intrinsic_process_time_rate",
            "symbolic_rule": "dX/dt = C(z) * dX/dtau_proc",
            "use_for": "კანონები, რომლებიც პირდაპირ tau_proc-ზეა განსაზღვრული",
            "double_counting_check": "გამოიყენე ზუსტად ერთი C(z) ფაქტორი",
        },
        {
            "rate_tag": "local_proper_time_rate",
            "symbolic_rule": "გადაყვანა მხოლოდ local clock calibration-ის მითითების შემდეგ",
            "use_for": "ადგილობრივი decay/oscillon/tail emission ერთ ეპოქაში",
            "double_counting_check": "არ დაუმატო process-time ფაქტორი, თუ კანონი tau_proc-ზე არ გადაიწერა",
        },
        {
            "rate_tag": "coordinate_cosmic_time_rate",
            "symbolic_rule": "დამატებითი C(z) არ შედის",
            "use_for": "FLRW background, Boltzmann evolution, merger crossing times",
            "double_counting_check": "dot უკვე დროის ცვლადს აფიქსირებს",
        },
        {
            "rate_tag": "observed_time_rate",
            "symbolic_rule": "რჩება სტანდარტული metric redshift/time dilation",
            "use_for": "დაკვირვებული სპექტრები, light curves, inferred rates",
            "double_counting_check": "observed redshift stretch არ გადაითარგმნოს process-time-ად",
        },
    ]


def process_time_application_map():
    """process-time ledger-ის გამოყენების რუკა სხვა კოსმოლოგიურ ბლოკებთან."""
    return {
        "MOND_vortex_maturation": (
            "coordinate-time H_eff რჩება ძირითად კალიბრაციად; process weighting "
            "მხოლოდ მონიშნულ intrinsic formation-rate კანონებში შედის"
        ),
        "JWST_high_z_galaxies": (
            "self-similar postulate-ის ქვეშ მოკლე LCDM ფანჯარა შეიძლება "
            "გადაითარგმნოს დიდ intrinsic process-time ბიუჯეტად; ეს არის "
            "formation-clock support calculation, არა solved abundance claim"
        ),
        "resonant_tail_dark_energy": (
            "tail-power C^4 არის bookkeeping input; P_tail-ის დროის ცვლადი "
            "ინტეგრაციამდე უნდა გამოცხადდეს"
        ),
        "cluster_mergers": (
            "merger crossing და redistribution times ადგილობრივი ეპოქის სიდიდეებია; "
            "არ გაამრავლო Bullet/cluster crossing times C(z)-ზე"
        ),
        "CMB": (
            "primary CMB იყენებს locked FLRW branch-ს; C(z) არ არის დამატებითი "
            "background metric mode"
        ),
        "BBN_atomic_clocks": (
            "standard nuclear/atomic clocks are blocked channels unless a "
            "separate varying-constants/microphysics bridge is derived and bounded"
        ),
    }


def stage_a2_process_time_status():
    """p02b-ში დარჩენილი process-time ledger-ის სრული სტატუსი."""
    return {
        "source_provenance": "OLD/18-derived ledger; not active authority",
        "source_requirements": process_time_source_requirements(),
        "self_similar_postulate": process_time_self_similar_postulate(),
        "cosmic_age_invariance": cosmic_age_invariance_result(),
        "Cz_phi_bridge": self_similar_Cz_phi_bridge(),
        "process_rate_vs_lapse": process_rate_vs_lapse_separation(),
        "channel_law": process_time_channel_law(),
        "channel_audit_examples": process_time_channel_audit_examples(),
        "clock_pressure_index": clock_pressure_index_definition(),
        "clock_pressure_lapse_bridge": clock_pressure_lapse_bridge(),
        "schwarzschild_reference": schwarzschild_clock_pressure_reference(),
        "age_conversion_example": process_time_age_conversion_example(),
        "JWST_diagnostic_table": jwst_process_time_diagnostic_table(),
        "JWST_interval_table": jwst_process_time_interval_table(),
        "JWST_budget_test": jwst_formation_budget_test(),
        "scaling": process_time_scaling_ledger(),
        "flrw_separation": flrw_vs_process_time_separation(),
        "allowed_rate_tags": process_time_allowed_rate_tags(),
        "observational_anchors": process_time_observational_anchor_requirements(),
        "energy_budget": process_time_energy_budget_guardrail(),
        "integrals": process_time_integrals(),
        "rate_table": rate_conversion_no_double_counting(),
        "applications": process_time_application_map(),
        "integration_status": "separate ledger in p02b; not part of p02 FLRW proof",
    }


def process_time_claim_gate():
    """ცხადი claim gate: რას აკეთებს p02b და რას კრძალავს."""
    separation = flrw_vs_process_time_separation()
    source = process_time_source_requirements()
    energy = process_time_energy_budget_guardrail()
    postulate = process_time_self_similar_postulate()
    age_result = cosmic_age_invariance_result()
    cz_bridge = self_similar_Cz_phi_bridge()
    channel_law = process_time_channel_law()
    pressure_index = clock_pressure_index_definition()
    lapse_bridge = clock_pressure_lapse_bridge()
    budget_test = jwst_formation_budget_test()

    return {
        "metric_independence": "PASS" if separation["not_a_second_metric_mode"] else "FAIL",
        "C_z_source": source["C_of_z_source"],
        "self_similar_curve": postulate["status"],
        "age_invariance": age_result["status"],
        "Cz_phi_bridge": cz_bridge["status"],
        "channel_law": channel_law["status"],
        "clock_pressure_index": pressure_index["status"],
        "clock_pressure_lapse_bridge": lapse_bridge["status"],
        "energy_budget": energy["status"],
        "intrinsic_rate_use": "CONDITIONAL_ONLY_WITH_EXPLICIT_RATE_TAG",
        "JWST_formation_time_support": budget_test["status"],
        "observational_status": "OPEN_BOUNDS_AND_FIT",
        "allowed_example": process_time_use_gate("intrinsic_process_time_rate"),
        "blocked_metric_example": process_time_use_gate(
            "metric_FLRW_background",
            enters_metric_branch=True,
        ),
        "do_not_claim": [
            "do not claim JWST high-z galaxy solution without abundance/mass/metallicity fit",
            "do not claim H0-tension solution",
            "do not claim dark-energy mechanism",
            "do not claim BBN/CMB compatibility",
            "do not apply C(z) to BBN, atomic, photon, or stellar microphysics clocks",
            "do not read process time as tired-light or redshift replacement",
            "do not call Pi_clock physical pressure before the stress/lapse bridge",
            "do not promote 13.8 Gyr to a new fundamental constant; keep T0 symbolic",
            "do not insert C(z) into Einstein-Boltzmann background",
            "do not multiply observed redshift time dilation by C(z)",
        ],
    }


if __name__ == "__main__":
    print("=" * 72)
    print("p02b: process-time ledger")
    print("=" * 72)

    relation, weak_field, *_ = get_process_time_relation()
    print("process-time relation:", relation)
    print("weak-field limit:", weak_field)
    print("FLRW/CMB separation:", flrw_vs_process_time_separation()["not_a_second_metric_mode"])
    print("C(z) source status:", process_time_source_requirements()["C_of_z_source"])
    print("self-similar status:", process_time_self_similar_postulate()["status"])
    print("age invariance:", cosmic_age_invariance_result()["status"])
    print("C(z)-phi bridge:", self_similar_Cz_phi_bridge()["status"])
    print("channel law:", process_time_channel_law()["status"])
    print("clock-pressure index:", clock_pressure_index_definition()["status"])
    print("clock-pressure lapse bridge:", clock_pressure_lapse_bridge()["status"])
    print("Schwarzschild reference:", schwarzschild_clock_pressure_reference()["status"])
    print(
        "2 Gyr future age in future units:",
        process_time_age_conversion_example()["future_age_measured_in_future_units_gyr"],
    )
    print("JWST interval diagnostics:")
    for row in jwst_process_time_interval_table()["rows"]:
        print(" ", row)
    print("JWST budget test:", jwst_formation_budget_test()["status"])
    print("energy budget:", process_time_energy_budget_guardrail()["status"])
    print("BBN channel gate:", process_time_channel_gate("BBN_nuclear_reaction_clock")["status"])
    print("intrinsic rate gate:", process_time_use_gate("intrinsic_process_time_rate")["status"])
    print("metric rate gate:", process_time_use_gate("metric_FLRW_background", enters_metric_branch=True)["status"])
    print("claim gate:", process_time_claim_gate()["observational_status"])
    print("integration status:", stage_a2_process_time_status()["integration_status"])
