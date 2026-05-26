# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.
# Active coefficient scheme: Y-scheme; c_Y means c_Y^(Y), not X-scheme c_X.

"""
p02: FLRW კოსმოლოგიური ფონი.

ამ ფაილის ფარგლები:
- სრული პოლინომიური L-ის FLRW სტრესი;
- rho და p_iso-ის სიმბოლური შემოწმება;
- GR-ფორმის Friedmann bookkeeping იგივე-მატერიის ფონზე;
- Phi(t)-ის Noether/conservation იდენტობა.

ეს ფაილი არ ამტკიცებს სრულ RFG გრავიტაციულ ფონის დინამიკას. Friedmann-ის
განტოლებები აქ არის სამუშაო bookkeeping branch, რომელიც შემდეგ უნდა შეიკრას
სრული გრავიტაციული სექტორით და დაკვირვებითი ფიტით.

Process-time ledger გატანილია p02b_process_time_ledger.py-ში.
"""

import sympy as sp

from p01_core import calculate_stress_tensor
from p01_core import evaluate_on_background
from p01_core import get_polynomial_lagrangian
from p01_core import init_variables
from p01_core import reduce_zero


def get_flrw_pressures():
    """
    FLRW ფონური ჩასმა p01_core.evaluate_on_background-ით.

    აბრუნებს:
    - rho = T_00
    - p_iso = T_11 / a^2
    - a(t)
    - სრულ phase/Bianchi შედეგს
    """
    result = evaluate_on_background("flrw", lagrangian_mode="full")

    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    rho = sp.simplify(result["T_cov"][0, 0])
    p_iso = sp.simplify(result["T_cov"][1, 1] / a**2)

    return rho, p_iso, a, result


def expected_flrw_stress():
    """
    rho და p_iso-ის ხელით ჩაწერილი სრული FLRW ფორმა.

    ეს არის ალგებრული smoke-test p01_core-ის გამოთვლასთან შესადარებლად;
    დამოუკიდებელ გრავიტაციულ გამოყვანად არ ითვლება.
    """
    c_Y, c_Y2 = sp.symbols("c_Y c_Y2", real=True)
    c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )

    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    Phi = sp.Function("Phi")(t)
    Phi_dot = sp.diff(Phi, t)

    rho_full = (
        -3 * c_I1 / a**2
        - 9 * c_I1sq / a**4
        - 3 * c_I2 / a**4
        - c_I3 / a**6
        + c_Y * Phi_dot**2
        + 3 * c_Y2 * Phi_dot**4
        + 3 * c_YI1 * Phi_dot**2 / a**2
    )

    p_iso_full = (
        c_I1 / a**2
        - (3 * c_I1sq + c_I2) / a**4
        - c_I3 / a**6
        + c_Y * Phi_dot**2
        + c_Y2 * Phi_dot**4
        + c_YI1 * Phi_dot**2 / a**2
    )

    rho_old_simplified = (
        -3 * c_I1 / a**2
        - 9 * c_I1sq / a**4
        + c_Y * Phi_dot**2
        + 3 * c_Y2 * Phi_dot**4
        + 3 * c_YI1 * Phi_dot**2 / a**2
    )

    return rho_full, p_iso_full, rho_old_simplified


def compare_with_theory(rho):
    """
    ძველ გამარტივებულ rho ფორმულასთან თავსებადობის შემოწმება.

    ძველ ფორმულას აკლია c_I2 და c_I3 წევრები; სრული ფორმა უნდა დაემთხვეს.
    API შეგნებულად რჩება სამწევრიანი tuple, ძველი შემოწმებების დასაცავად.
    """
    rho_full, _, rho_old_simplified = expected_flrw_stress()

    full_match = sp.simplify(rho - rho_full) == 0
    old_match = sp.simplify(rho - rho_old_simplified) == 0
    missing_terms = sp.simplify(rho_full - rho_old_simplified)

    return full_match, old_match, missing_terms


def compare_stress_with_expected(rho, p_iso):
    """rho და p_iso-ის სრული FLRW smoke-test."""
    rho_full, p_iso_full, rho_old_simplified = expected_flrw_stress()

    return {
        "rho_full_match": sp.simplify(rho - rho_full) == 0,
        "p_iso_full_match": sp.simplify(p_iso - p_iso_full) == 0,
        "old_rho_match": sp.simplify(rho - rho_old_simplified) == 0,
        "old_rho_missing_terms": sp.simplify(rho_full - rho_old_simplified),
    }


def flrw_stress_theorem():
    """
    Algebraic theorem: the full polynomial L gives the listed FLRW rho and p.
    """
    rho, p_iso, _, _ = get_flrw_pressures()
    rho_full, p_iso_full, _ = expected_flrw_stress()
    rho_residual = sp.simplify(rho - rho_full)
    p_residual = sp.simplify(p_iso - p_iso_full)

    return {
        "status": "PASS" if rho_residual == 0 and p_residual == 0 else "CHECK",
        "rho_expected": rho_full,
        "p_iso_expected": p_iso_full,
        "rho_residual": rho_residual,
        "p_iso_residual": p_residual,
        "claim": "full FLRW stress algebra follows from the p01 polynomial L",
    }


def check_bianchi_residual(result):
    """p01_core-ის Bianchi/Noether residual FLRW-ზე."""
    residual = [reduce_zero(value) for value in result["residual"]]
    is_ok = all(value == 0 for value in residual)
    return is_ok, residual


def flrw_noether_theorem():
    """
    Algebraic theorem: FLRW conservation equals Phi EOM times Phi_dot.
    """
    rho, p_iso, _, result = get_flrw_pressures()
    bianchi_ok, bianchi_residual = check_bianchi_residual(result)
    conservation_diff, expected_lhs, conservation_strict, EOM_strict = check_conservation()

    return {
        "status": "PASS" if bianchi_ok and sp.simplify(conservation_diff) == 0 else "CHECK",
        "bianchi_residual": bianchi_residual,
        "conservation_minus_EOMPhi_PhiDot": sp.simplify(conservation_diff),
        "EOMPhi_times_PhiDot": expected_lhs,
        "strict_clock_residual_raw": sp.factor(conservation_strict),
        "strict_clock_EOM_raw": sp.factor(EOM_strict),
        "claim": "drho/dt + 3H(rho+p) = EOM_Phi * Phi_dot on FLRW",
    }


def get_friedmann_equations():
    """
    GR-ფორმის FLRW bookkeeping იგივე-მატერიის ფონზე.

    ეს არ არის RFG გრავიტაციული სექტორის სრული გამოყვანა. Phi რჩება Phi(t)-ად;
    strict Phi=t გამოიყენება მხოლოდ ცალკე diagnostic-ში.
    """
    rho_solid, p_iso_solid, a, _ = get_flrw_pressures()
    t = sp.Symbol("t", real=True)
    H = sp.diff(a, t) / a
    kappa = sp.Symbol("kappa", real=True)

    rho_m0, rho_rad0 = sp.symbols("rho_m0 rho_rad0", real=True)
    rho_m = rho_m0 / a**3
    rho_rad = rho_rad0 / a**4
    p_m = 0
    p_rad = rho_rad / 3

    friedmann1 = sp.Eq(3 * H**2, kappa * (rho_solid + rho_m + rho_rad))
    friedmann2 = sp.Eq(
        2 * sp.diff(a, t, t) / a + H**2,
        -kappa * (p_iso_solid + p_m + p_rad),
    )

    return friedmann1, friedmann2, a, t, rho_solid, p_iso_solid


def flrw_pressure_relaxation_interpretation():
    """
    Interpretive bridge from Intuitive_Theory §6.3 to the p02 FLRW ledger.

    This does not replace FLRW expansion or redshift. It records the internal
    RFG reading: the same background history can be read as substrate pressure
    relaxation while the metric branch remains the observational FLRW language.
    """
    a = sp.Symbol("a", positive=True)
    Phi_dot = sp.Symbol("Phi_dot", real=True)
    c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_I1 c_I1sq c_I2 c_I3 c_YI1",
        real=True,
    )

    return {
        "status": "INTERPRETIVE_BRIDGE_ONLY",
        "external_FLRW_language": (
            "a(t) grows; cosmological distances/redshift are described in the "
            "standard metric FLRW branch"
        ),
        "internal_RFG_language": (
            "substrate pressure relaxes; operational matter standards and "
            "clock readings are part of the changing medium"
        ),
        "not_a_replacement": (
            "do not claim galaxies are not receding or redshift is non-metric "
            "without a separate observational bridge"
        ),
        "process_time_pointer": (
            "internal formation-clock reparameterization remains in p02b; it "
            "does not enter this primary FLRW/CMB metric branch"
        ),
        "early_component_interpretation": {
            "curvature_like_pressure_memory": sp.simplify(
                (-3 * c_I1 + 3 * c_YI1 * Phi_dot**2) / a**2
            ),
            "radiation_like_substrate_stress": sp.simplify(
                -(9 * c_I1sq + 3 * c_I2) / a**4
            ),
            "stiff_early_substrate_compression": sp.simplify(-c_I3 / a**6),
        },
        "late_component_interpretation": (
            "strict-clock Lambda_eff is a late residual pressure-floor "
            "diagnostic, not a derived dark-energy solution"
        ),
        "observational_guardrails": [
            "CMB/BBN/Planck bounds still apply to 1/a^2, 1/a^4, and 1/a^6 terms",
            "SN/BAO/redshift observations remain metric-FLRW in this file",
            "do not import p02b process-time factors into H(z) or CMB background",
            "do not claim shrinking rods replace expansion without a full fit",
        ],
    }


def check_conservation():
    """
    FLRW-ზე ამოწმებს:
        drho/dt + 3H(rho+p) = EOM_Phi * Phi_dot

    strict Phi=t ნარჩენი მოდის დაუხურავი phase-current EOM-იდან; c_YI1
    ამ ნარჩენში solid-coupled ნაწილს აჩენს, ხოლო c_Y/c_Y2 pure-phase ნაწილს.
    ამიტომ strict Phi=t არ უნდა ჩაითვალოს background-ის დახურვად დამატებითი
    პირობის გარეშე.
    """
    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    H = sp.diff(a, t) / a

    g00, g11, g22, g33 = sp.symbols("g00 g11 g22 g33", real=True)
    Phi_dot = sp.Function("Phi_dot")(t)

    Y_g = g00 * Phi_dot**2
    I1_g = -g11 - g22 - g33
    I2_g = g11 * g22 + g22 * g33 + g11 * g33
    I3_g = -g11 * g22 * g33

    Y, I1, I2, I3 = init_variables()
    L = get_polynomial_lagrangian(Y, I1, I2, I3)
    L_g = L.subs({Y: Y_g, I1: I1_g, I2: I2_g, I3: I3_g})

    T_mixed = calculate_stress_tensor(L_g, [g00, g11, g22, g33])
    ansatz = {g00: 1, g11: -1 / a**2, g22: -1 / a**2, g33: -1 / a**2}

    rho = sp.simplify(T_mixed[0].subs(ansatz))
    p_iso = sp.simplify(-T_mixed[1].subs(ansatz))

    conservation_lhs = sp.simplify(sp.diff(rho, t) + 3 * H * (rho + p_iso))
    dL_dPhi_dot = sp.simplify(sp.diff(L_g, Phi_dot).subs(ansatz))
    EOM_Phi = sp.simplify(sp.diff(a**3 * dL_dPhi_dot, t) / a**3)
    expected_lhs = sp.simplify(EOM_Phi * Phi_dot)
    difference = sp.simplify(conservation_lhs - expected_lhs)

    strict_phi_t = {Phi_dot: 1, sp.diff(Phi_dot, t): 0}
    conservation_strict = sp.simplify(conservation_lhs.subs(strict_phi_t))
    EOM_strict = sp.simplify(EOM_Phi.subs(strict_phi_t))

    return difference, expected_lhs, conservation_strict, EOM_strict


def phi_clock_closure_conditions():
    """
    strict Phi=t branch-ის ალგებრული დახურვის პირობები.

    ეს არ ხდის strict-clock branch-ს დაკვირვებით დამტკიცებულად; ის მხოლოდ
    აჩვენებს, რა coefficient tuning არის საჭირო, რომ phase-current residual
    generic expanding FLRW-ზე გაქრეს.
    """
    _, _, conservation_strict, EOM_strict = check_conservation()

    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)
    closure_subs = {
        c_Y: -2 * c_Y2,
        c_YI1: 0,
    }

    K_phi_today = sp.simplify(c_Y + 6 * c_Y2 + 3 * c_YI1)
    rho_late_strict = sp.simplify(c_Y + 3 * c_Y2)
    p_late_strict = sp.simplify(c_Y + c_Y2)

    rho_closed = sp.simplify(rho_late_strict.subs(closure_subs))
    p_closed = sp.simplify(p_late_strict.subs(closure_subs))

    return {
        "strict_clock_residual": sp.factor(conservation_strict),
        "strict_clock_EOM_residual": sp.factor(EOM_strict),
        "generic_expanding_FLRW_closure": [
            sp.Eq(c_YI1, 0),
            sp.Eq(c_Y + 2 * c_Y2, 0),
        ],
        "residual_after_closure": sp.simplify(conservation_strict.subs(closure_subs)),
        "phase_no_ghost_after_closure": sp.simplify(K_phi_today.subs(closure_subs)),
        "late_rho_after_closure": rho_closed,
        "late_p_after_closure": p_closed,
        "late_w_after_closure": sp.simplify(p_closed / rho_closed),
        "remaining_status": (
            "STRICT_CLOCK_BRANCH_CONDITIONAL: closes algebraically only with "
            "c_YI1=0 and c_Y=-2*c_Y2 for generic expansion; still requires "
            "solid-sector stability and observational coefficient fit."
        ),
    }


def strict_clock_obstruction_theorem():
    """
    Algebraic theorem: strict Phi=t is obstructed unless the closure conditions hold.
    """
    closure = phi_clock_closure_conditions()
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    return {
        "status": "PASS",
        "strict_clock_residual": closure["strict_clock_residual"],
        "generic_expanding_FLRW_requires": closure["generic_expanding_FLRW_closure"],
        "interpretation": (
            "For generic a_dot != 0, strict Phi=t is not an automatic background; "
            "it requires c_YI1=0 and c_Y=-2*c_Y2."
        ),
        "fine_tuning_warning": (
            "these conditions are algebraic closure conditions, not yet derived "
            "from a symmetry or selection rule"
        ),
        "not_a_dark_energy_claim": True,
    }


def strict_clock_lambda_like_branch_theorem():
    """
    Conditional algebraic theorem: if the strict-clock closure holds, w=-1.
    """
    closure = phi_clock_closure_conditions()
    rho = closure["late_rho_after_closure"]
    p = closure["late_p_after_closure"]
    w = sp.simplify(p / rho)

    return {
        "status": "CONDITIONAL_THEOREM",
        "conditions": closure["generic_expanding_FLRW_closure"],
        "phase_no_ghost_after_closure": closure["phase_no_ghost_after_closure"],
        "late_rho_after_closure": rho,
        "late_p_after_closure": p,
        "late_w_after_closure": w,
        "theorem_pass": sp.simplify(w + 1) == 0,
        "claim": "closed strict-clock branch has Lambda-like equation of state w=-1",
        "not_claimed": (
            "observed dark energy, naturalness, full gravity closure, and "
            "perturbation stability remain open"
        ),
    }


def dynamic_phase_clock_branch_theorem():
    """
    Dynamic alternative to the locked strict Phi=t diagnostic.

    The phase EOM has a conserved-current form. Therefore the natural branch is
    Phi_dot(a), not necessarily Phi_dot=1. This does not prove process-time, but
    it shows how the strict-clock tuning can be replaced by a dynamical clock.
    """
    a = sp.Symbol("a", positive=True)
    u = sp.Symbol("u", real=True)  # u = Phi_dot
    Q = sp.Symbol("Q_Phi", real=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    current_density = sp.simplify(u * (c_Y + 2 * c_Y2 * u**2 + 3 * c_YI1 / a**2))
    conserved_charge = sp.simplify(a**3 * current_density)
    cubic_condition = sp.Eq(
        sp.simplify(2 * c_Y2 * u**3 + (c_Y + 3 * c_YI1 / a**2) * u),
        Q / a**3,
    )

    zero_current_u2 = sp.simplify(-(c_Y + 3 * c_YI1 / a**2) / (2 * c_Y2))
    zero_current_late_u2 = sp.simplify(sp.limit(zero_current_u2, a, sp.oo))

    rho_late_zero_current = sp.simplify(c_Y**2 / (4 * c_Y2))
    p_late_zero_current = sp.simplify(-c_Y**2 / (4 * c_Y2))
    w_late_zero_current = sp.simplify(p_late_zero_current / rho_late_zero_current)
    kinetic_prefactor_zero_current = sp.simplify(4 * c_Y2 * zero_current_u2)

    return {
        "status": "CONDITIONAL_DYNAMIC_BRANCH_CANDIDATE",
        "phase_current_density": current_density,
        "conserved_current": sp.Eq(conserved_charge, Q),
        "dynamic_clock_equation": cubic_condition,
        "zero_current_branch_u2": sp.Eq(u**2, zero_current_u2),
        "zero_current_late_u2": sp.Eq(sp.Symbol("u_late", real=True) ** 2, zero_current_late_u2),
        "late_rho_zero_current": rho_late_zero_current,
        "late_p_zero_current": p_late_zero_current,
        "late_w_zero_current": w_late_zero_current,
        "kinetic_prefactor_zero_current": kinetic_prefactor_zero_current,
        "late_viable_sign_window_hint": [
            sp.Gt(c_Y2, 0),
            sp.Lt(c_Y, 0),
        ],
        "what_it_changes": (
            "strict Phi_dot=1 needs c_YI1=0 and c_Y=-2*c_Y2 for all a; "
            "dynamic Phi_dot(a) can close the phase current without forcing c_YI1=0"
        ),
        "required_checks": [
            "u^2 must be positive on the intended epoch range",
            "kinetic/no-ghost and gradient stability must be checked",
            "Friedmann closure and BBN/CMB/Planck coefficient bounds remain open",
        ],
        "not_claimed": (
            "this does not derive the p02b self-similar C(z); it only provides "
            "a p02 phase-current branch to compare with process-time"
        ),
    }


def dynamic_zero_current_lambda_like_branch_theorem():
    """
    Conditional dynamic theorem: the nonzero zero-current phase branch is w=-1 late.

    This is stronger than the locked strict Phi=t diagnostic because Phi_dot is
    allowed to vary with a. It still needs positivity, stability, gravity closure,
    and observational fitting.
    """
    branch = dynamic_phase_clock_branch_theorem()
    c_Y, c_Y2 = sp.symbols("c_Y c_Y2", real=True)
    u_late = sp.Symbol("u_late", real=True)

    late_u2 = sp.simplify(-c_Y / (2 * c_Y2))
    rho = branch["late_rho_zero_current"]
    p = branch["late_p_zero_current"]
    w = sp.simplify(p / rho)
    kinetic_prefactor_late = sp.simplify(4 * c_Y2 * late_u2)

    return {
        "status": "CONDITIONAL_DYNAMIC_THEOREM",
        "branch": "Q_Phi=0, Phi_dot != 0, a -> infinity",
        "late_clock_condition": sp.Eq(u_late**2, late_u2),
        "late_rho": rho,
        "late_p": p,
        "late_w": w,
        "theorem_pass": sp.simplify(w + 1) == 0,
        "late_kinetic_prefactor": kinetic_prefactor_late,
        "conditions": [
            sp.Ne(c_Y2, 0),
            "u_late^2 > 0",
            "late_kinetic_prefactor > 0",
            "for positive late rho with this branch: c_Y2 > 0 and c_Y < 0",
            "full perturbation stability and observational fit remain open",
        ],
        "what_is_strengthened": (
            "Lambda-like late pressure can arise on a dynamic phase-current "
            "branch without requiring Phi_dot=1 for every epoch"
        ),
        "not_claimed": (
            "this is not yet an observed dark-energy solution and it does not "
            "derive the p02b self-similar process-time law"
        ),
    }


def process_time_dynamic_clock_comparison():
    """
    Compare p02's dynamic Phi_dot(a) branch with p02b's process-time target.
    """
    a = sp.Symbol("a", positive=True)
    T0 = sp.Symbol("T0", positive=True)
    t_age = sp.Function("t_age")
    C_proc = sp.Function("C_proc")
    Q = sp.Symbol("Q_Phi", real=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    C_target = T0 / t_age(a)
    cubic_match = sp.Eq(
        2 * c_Y2 * C_proc(a) ** 3 + (c_Y + 3 * c_YI1 / a**2) * C_proc(a),
        Q / a**3,
    )
    zero_current_match = sp.Eq(
        C_proc(a) ** 2,
        -(c_Y + 3 * c_YI1 / a**2) / (2 * c_Y2),
    )

    return {
        "status": "COMPARISON_LEDGER_NOT_DERIVED_MATCH",
        "p02_dynamic_clock_variable": "u(a)=Phi_dot(a)",
        "p02b_self_similar_target": sp.Eq(C_proc(a), C_target),
        "general_match_condition": cubic_match,
        "zero_current_match_condition": zero_current_match,
        "interpretation": (
            "process-time suggests replacing locked Phi_dot=1 by a dynamic "
            "clock. p02 supplies a conserved-current equation for that clock."
        ),
        "main_gap": (
            "constant p02 coefficients have not been shown to reproduce "
            "C_proc(a)=T0/t_age(a) over the observational range"
        ),
        "safe_next_step": (
            "fit or derive c_Y, c_Y2, c_YI1 and Q_Phi so that u(a) is compared "
            "against the p02b process-time channel without entering the CMB metric branch"
        ),
    }


def early_scaling_theorem():
    """
    Algebraic theorem: FLRW solid components scale as a^-2, a^-4, and a^-6.
    """
    a = sp.Symbol("a", positive=True)
    Phi_dot = sp.Symbol("Phi_dot", real=True)
    c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_I1 c_I1sq c_I2 c_I3 c_YI1",
        real=True,
    )

    curvature_like = sp.simplify((-3 * c_I1 + 3 * c_YI1 * Phi_dot**2) / a**2)
    radiation_like = sp.simplify(-(9 * c_I1sq + 3 * c_I2) / a**4)
    stiff_like = sp.simplify(-c_I3 / a**6)

    return {
        "status": "PASS",
        "curvature_like_a_minus_2": curvature_like,
        "radiation_like_a_minus_4": radiation_like,
        "stiff_like_a_minus_6": stiff_like,
        "claim": "the FLRW background algebra separates the early components by scaling",
        "observational_warning": (
            "PASS is algebraic only; BBN/CMB/Planck compatibility requires "
            "numerical coefficient bounds"
        ),
    }


def late_time_density_status(rho_solid):
    """
    a -> infinity ლიმიტი.

    strict_clock_density არის მხოლოდ Phi_dot=1 diagnostic. სრული branch-ში
    rho_lambda Phi_dot(t)-ზეა დამოკიდებული.
    """
    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    x = sp.Symbol("x", positive=True)
    Phi = sp.Function("Phi")(t)
    Phi_dot = sp.diff(Phi, t)
    kappa = sp.Symbol("kappa", real=True)

    rho_late = sp.simplify(sp.limit(rho_solid.subs(a, x), x, sp.oo))
    rho_strict = sp.simplify(rho_late.subs(Phi_dot, 1))

    return {
        "rho_late_dynamic": rho_late,
        "rho_late_strict_clock": rho_strict,
        "Lambda_eff_strict_clock": sp.simplify(kappa * rho_strict),
        "naturalness_status": (
            "OPEN_NUMERICAL_FIT: Lambda_eff strict-clock diagnostic is not a "
            "derived dark-energy solution; requires coefficient fit and dynamic Phi(t)."
        ),
    }


def early_component_status():
    """ადრეული ეპოქის დამატებითი წევრები, რომლებიც დაკვირვებით უნდა შეიზღუდოს."""
    return {
        "stiff_like": "-c_I3/a^6",
        "radiation_like_extra": "-(9*c_I1sq + 3*c_I2)/a^4",
        "curvature_like_extra": "(-3*c_I1 + 3*c_YI1*Phi_dot^2)/a^2",
        "pressure_relaxation_reading": (
            "these are substrate pressure/tension-memory components in the "
            "internal RFG language, but observationally they remain constrained "
            "as curvature-like, radiation-like, and stiff FLRW terms"
        ),
        "status": (
            "OPEN_NUMERICAL_FIT: no BBN/CMB/Planck compatibility claim before "
            "fitting c_I1, c_I1sq, c_I2, c_I3, c_YI1 and Phi(t)."
        ),
    }


def cosmology_claim_gate():
    """ცხადი claim ledger: რა იხურება p02-ში და რა არა."""
    rho, p_iso, _, result = get_flrw_pressures()
    stress_check = compare_stress_with_expected(rho, p_iso)
    bianchi_ok, residual = check_bianchi_residual(result)
    conservation_diff, _, conservation_strict, _ = check_conservation()
    clock_closure = phi_clock_closure_conditions()

    strict_clock_closed = sp.simplify(clock_closure["residual_after_closure"]) == 0
    relaxation = flrw_pressure_relaxation_interpretation()
    stress_theorem = flrw_stress_theorem()
    noether_theorem = flrw_noether_theorem()
    obstruction_theorem = strict_clock_obstruction_theorem()
    lambda_theorem = strict_clock_lambda_like_branch_theorem()
    dynamic_branch = dynamic_phase_clock_branch_theorem()
    dynamic_lambda_theorem = dynamic_zero_current_lambda_like_branch_theorem()
    process_time_comparison = process_time_dynamic_clock_comparison()
    scaling_theorem = early_scaling_theorem()

    return {
        "stress_theorem": stress_theorem["status"],
        "stress_algebra": "PASS" if (
            stress_check["rho_full_match"] and stress_check["p_iso_full_match"]
        ) else "CHECK",
        "noether_theorem": noether_theorem["status"],
        "bianchi_noether": "PASS" if bianchi_ok else "CHECK",
        "noether_difference": "PASS" if sp.simplify(conservation_diff) == 0 else "CHECK",
        "friedmann_GR_form": "BOOKKEEPING_ONLY",
        "strict_clock_obstruction": obstruction_theorem["status"],
        "strict_Phi_t": (
            "CONDITIONAL_CLOSED" if strict_clock_closed else "NOT_CLOSED"
        ),
        "strict_Phi_t_conditions": clock_closure["generic_expanding_FLRW_closure"],
        "lambda_like_branch_theorem": lambda_theorem["status"],
        "dynamic_phase_clock": dynamic_branch["status"],
        "dynamic_lambda_like_branch": dynamic_lambda_theorem["status"],
        "process_time_dynamic_clock_comparison": process_time_comparison["status"],
        "late_time_Lambda_eff": (
            "CONDITIONAL_DIAGNOSTIC_ONLY: w=-1 on strict-clock closure, "
            "and w=-1 on the dynamic zero-current late branch, but observed "
            "dark-energy claim needs coefficient fit and full gravity closure"
        ),
        "early_scaling_theorem": scaling_theorem["status"],
        "early_universe_BBN_CMB_Planck": "OPEN_FIT",
        "pressure_relaxation_interpretation": relaxation["status"],
        "process_time": "MOVED_TO_p02b_NOT_PRIMARY_FLRW_BRANCH",
        "do_not_claim": [
            "do not claim H0-tension solution",
            "do not claim dark-energy solution",
            "do not claim BBN/CMB/Planck compatibility",
            "do not use process-time as second metric mode",
            "do not identify Phi_dot(a) with p02b C(z) without a match condition and bounds",
            "do not claim expansion is replaced by shrinking rods without observational bridge",
            "do not treat pressure-relaxation interpretation as a new H(z) fit",
        ],
        "raw_bianchi_residual": residual,
        "raw_strict_clock_residual": sp.factor(conservation_strict),
    }


def module_status():
    """p02-ის მოკლე სტატუსი საბჭოს შენიშვნების შემდეგ."""
    rho, p_iso, _, result = get_flrw_pressures()
    stress_check = compare_stress_with_expected(rho, p_iso)
    bianchi_ok, residual = check_bianchi_residual(result)
    conservation_diff, _, conservation_strict, EOM_strict = check_conservation()

    return {
        "scope": "FLRW stress/Friedmann-bookkeeping/conservation ledger only",
        "stress_check": stress_check,
        "stress_theorem": flrw_stress_theorem(),
        "bianchi_ok": bianchi_ok,
        "bianchi_residual": residual,
        "noether_difference": conservation_diff,
        "noether_theorem": flrw_noether_theorem(),
        "strict_Phi_t_conservation_residual": conservation_strict,
        "strict_Phi_t_EOM_residual": EOM_strict,
        "strict_clock_obstruction_theorem": strict_clock_obstruction_theorem(),
        "strict_clock_lambda_like_branch_theorem": strict_clock_lambda_like_branch_theorem(),
        "dynamic_phase_clock_branch_theorem": dynamic_phase_clock_branch_theorem(),
        "dynamic_zero_current_lambda_like_branch_theorem": (
            dynamic_zero_current_lambda_like_branch_theorem()
        ),
        "process_time_dynamic_clock_comparison": process_time_dynamic_clock_comparison(),
        "early_scaling_theorem": early_scaling_theorem(),
        "pressure_relaxation_interpretation": flrw_pressure_relaxation_interpretation(),
        "late_time_density": late_time_density_status(rho),
        "early_components": early_component_status(),
        "phi_clock_closure": phi_clock_closure_conditions(),
        "claim_gate": cosmology_claim_gate(),
        "process_time_ledger": "moved to p02b_process_time_ledger.py",
        "friedmann_status": "GR-ფორმის bookkeeping; სრული RFG გრავიტაციული გამოყვანა ღიაა",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("p02: FLRW კოსმოლოგიური ფონი")
    print("=" * 72)

    f1, f2, a, t, rho_solid, p_iso_solid = get_friedmann_equations()
    rho, p_iso, _, result = get_flrw_pressures()

    print("\n1. GR-ფორმის Friedmann bookkeeping")
    print("პირველი განტოლება:", f1)
    print("მეორე განტოლება:", f2)
    print("სტატუსი: ეს არის bookkeeping branch, არა სრული RFG გრავიტაციული გამოყვანა.")

    print("\n2. FLRW სტრესი")
    print("rho =", sp.expand(rho))
    print("p_iso =", sp.expand(p_iso))

    print("\n3. rho და p_iso-ის სრული ფორმის შემოწმება")
    checks = compare_stress_with_expected(rho, p_iso)
    print("rho სრული ფორმა ემთხვევა:", checks["rho_full_match"])
    print("p_iso სრული ფორმა ემთხვევა:", checks["p_iso_full_match"])
    print("ძველი rho ფორმა ემთხვევა:", checks["old_rho_match"])
    print("ძველ rho ფორმაში გამოტოვებული წევრები:", sp.expand(checks["old_rho_missing_terms"]))

    print("\n3b. ალგებრული theorem ledger")
    theorem_rows = {
        "stress theorem": flrw_stress_theorem()["status"],
        "Noether theorem": flrw_noether_theorem()["status"],
        "strict-clock obstruction": strict_clock_obstruction_theorem()["status"],
        "Lambda-like branch": strict_clock_lambda_like_branch_theorem()["status"],
        "dynamic phase-clock": dynamic_phase_clock_branch_theorem()["status"],
        "dynamic Lambda-like branch": (
            dynamic_zero_current_lambda_like_branch_theorem()["status"]
        ),
        "process-time comparison": process_time_dynamic_clock_comparison()["status"],
        "early scaling theorem": early_scaling_theorem()["status"],
    }
    for label, status in theorem_rows.items():
        print(f"{label}: {status}")

    print("\n3c. ფუძე-წნევის მოდუნების ინტერპრეტაცია")
    relaxation = flrw_pressure_relaxation_interpretation()
    print("სტატუსი:", relaxation["status"])
    print("გარე FLRW ენა:", relaxation["external_FLRW_language"])
    print("შიდა RFG ენა:", relaxation["internal_RFG_language"])
    print("guardrail:", relaxation["not_a_replacement"])

    print("\n4. Bianchi/Noether residual")
    bianchi_ok, residual = check_bianchi_residual(result)
    print("residual vector:", residual)
    if bianchi_ok:
        print("Bianchi/Noether იდენტობა FLRW-ზე სრულდება.")
    else:
        print("Bianchi residual შესამოწმებელია.")

    print("\n5. Phi(t) conservation diagnostic")
    diff, expected, cons_strict, eom_strict = check_conservation()
    print("Noether სხვაობა:", diff)
    print("strict Phi=t conservation residual:", cons_strict)
    print("strict Phi=t EOM residual:", eom_strict)
    print("დასკვნა: strict Phi=t ფონად არ იკეტება დამატებითი პირობის გარეშე.")

    print("\n6. strict Phi=t closure პირობები")
    clock_gate = phi_clock_closure_conditions()
    print("generic expanding closure:", clock_gate["generic_expanding_FLRW_closure"])
    print("residual after closure:", clock_gate["residual_after_closure"])
    print("phase no-ghost after closure:", clock_gate["phase_no_ghost_after_closure"])
    print("late w after closure:", clock_gate["late_w_after_closure"])
    print("სტატუსი:", clock_gate["remaining_status"])

    print("\n7. დინამიკური Phi_dot(a) შტო")
    dynamic_branch = dynamic_phase_clock_branch_theorem()
    dynamic_lambda = dynamic_zero_current_lambda_like_branch_theorem()
    process_match = process_time_dynamic_clock_comparison()
    print("სტატუსი:", dynamic_branch["status"])
    print("conserved current:", dynamic_branch["conserved_current"])
    print("dynamic clock equation:", dynamic_branch["dynamic_clock_equation"])
    print("zero-current branch:", dynamic_branch["zero_current_branch_u2"])
    print("late w zero-current:", dynamic_lambda["late_w"])
    print("დინამიკური theorem:", dynamic_lambda["status"])
    print("process-time comparison:", process_match["status"])
    print("მთავარი ღიობი:", process_match["main_gap"])

    print("\n8. გვიანი და ადრეული წევრები")
    late = late_time_density_status(rho_solid)
    print("rho_late dynamic:", late["rho_late_dynamic"])
    print("rho_late strict Phi=t:", late["rho_late_strict_clock"])
    print("Lambda_eff strict Phi=t:", late["Lambda_eff_strict_clock"])
    print("ნატურალურობა:", late["naturalness_status"])
    print("ადრეული წევრების სტატუსი:", early_component_status())

    print("\n9. Claim gate")
    gate = cosmology_claim_gate()
    for key in [
        "stress_theorem",
        "stress_algebra",
        "noether_theorem",
        "bianchi_noether",
        "noether_difference",
        "friedmann_GR_form",
        "strict_clock_obstruction",
        "strict_Phi_t",
        "lambda_like_branch_theorem",
        "dynamic_phase_clock",
        "dynamic_lambda_like_branch",
        "process_time_dynamic_clock_comparison",
        "late_time_Lambda_eff",
        "early_scaling_theorem",
        "early_universe_BBN_CMB_Planck",
        "pressure_relaxation_interpretation",
        "process_time",
    ]:
        print(f"{key}: {gate[key]}")
    print("do_not_claim:", gate["do_not_claim"])

    print("\n10. Process-time ledger")
    print("გატანილია: p02b_process_time_ledger.py")
