# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
p02: FLRW კოსმოლოგიური ფონი.

ამ ფაილის ფარგლები:
- სრული პოლინომიური L-ის FLRW სტრესი;
- rho და p_iso-ის სიმბოლური შემოწმება;
- GR-ფორმის Friedmann bookkeeping იგივე-მატერიის ფონზე;
- Phi(t)-ის Noether/conservation იდენტობა.
- process-time ledger კოსმოლოგიურ გამოყენებებთან ერთად.

ეს ფაილი არ ამტკიცებს სრულ RFG გრავიტაციულ ფონის დინამიკას. Friedmann-ის
განტოლებები აქ არის სამუშაო bookkeeping branch, რომელიც შემდეგ უნდა შეიკრას
სრული გრავიტაციული სექტორით და დაკვირვებითი ფიტით.

შერეული თემები ამ ფაილში შეგნებულად რჩება: p02 არის კოსმოლოგიის სამუშაო
საბჭო/ledger, არა მხოლოდ ერთი იზოლირებული ფორმულის ფაილი.
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


def check_bianchi_residual(result):
    """p01_core-ის Bianchi/Noether residual FLRW-ზე."""
    residual = [reduce_zero(value) for value in result["residual"]]
    is_ok = all(value == 0 for value in residual)
    return is_ok, residual


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


def check_conservation():
    """
    FLRW-ზე ამოწმებს:
        drho/dt + 3H(rho+p) = EOM_Phi * Phi_dot

    strict Phi=t ზოგადად ნარჩენს ტოვებს; ამიტომ იგი background-ის დახურვად
    არ უნდა ჩაითვალოს დამატებითი პირობის გარეშე.
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
        "naturalness_status": "არ არის გამოყვანილი",
    }


def early_component_status():
    """ადრეული ეპოქის დამატებითი წევრები, რომლებიც დაკვირვებით უნდა შეიზღუდოს."""
    return {
        "stiff_like": "-c_I3/a^6",
        "radiation_like_extra": "-(9*c_I1sq + 3*c_I2)/a^4",
        "curvature_like_extra": "(-3*c_I1 + 3*c_YI1*Phi_dot^2)/a^2",
        "status": "საჭიროა დაკვირვებითი ფიტი",
    }


def get_process_time_relation():
    """
    ისტორიული process-time ხაზობრივი ესკიზის დაკონსოლიდებული ფორმა.

    აქტიური normalization ემთხვევა p10-ის ბი-კონფორმულ სკალირებას:
    d tau / dt = exp(phi/2). ეს ბლოკი p02-ში რჩება კოსმოლოგიური ledger-ისთვის,
    მაგრამ არ შედის primary FLRW/CMB metric branch-ში.
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
        "status": "bookkeeping; დამოუკიდებელი metric mode არ არის",
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
            "შეიძლება შეამციროს საჭირო intrinsic development time მხოლოდ მაშინ, "
            "თუ formation rate process-controlled-ად არის მონიშნული"
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
    }


def stage_a2_process_time_status():
    """p02-ში დარჩენილი process-time ledger-ის სრული სტატუსი."""
    return {
        "source_provenance": "OLD/18-დან გადმოტანილი provenance; აქტიური ავტორიტეტი არ არის",
        "scaling": process_time_scaling_ledger(),
        "flrw_separation": flrw_vs_process_time_separation(),
        "integrals": process_time_integrals(),
        "rate_table": rate_conversion_no_double_counting(),
        "applications": process_time_application_map(),
        "integration_status": "დიზაინით რჩება p02-ში",
    }


def module_status():
    """p02-ის მოკლე სტატუსი საბჭოს შენიშვნების შემდეგ."""
    rho, p_iso, _, result = get_flrw_pressures()
    stress_check = compare_stress_with_expected(rho, p_iso)
    bianchi_ok, residual = check_bianchi_residual(result)
    conservation_diff, _, conservation_strict, EOM_strict = check_conservation()

    return {
        "scope": "FLRW კოსმოლოგია და process-time ledger",
        "stress_check": stress_check,
        "bianchi_ok": bianchi_ok,
        "bianchi_residual": residual,
        "noether_difference": conservation_diff,
        "strict_Phi_t_conservation_residual": conservation_strict,
        "strict_Phi_t_EOM_residual": EOM_strict,
        "late_time_density": late_time_density_status(rho),
        "early_components": early_component_status(),
        "process_time": stage_a2_process_time_status(),
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

    print("\n6. გვიანი და ადრეული წევრები")
    late = late_time_density_status(rho_solid)
    print("rho_late dynamic:", late["rho_late_dynamic"])
    print("rho_late strict Phi=t:", late["rho_late_strict_clock"])
    print("Lambda_eff strict Phi=t:", late["Lambda_eff_strict_clock"])
    print("ნატურალურობა:", late["naturalness_status"])
    print("ადრეული წევრების სტატუსი:", early_component_status())

    print("\n7. Process-time ledger")
    relation, weak_field, *_ = get_process_time_relation()
    print("process-time relation:", relation)
    print("weak-field limit:", weak_field)
    print("FLRW/CMB separation:", flrw_vs_process_time_separation()["not_a_second_metric_mode"])
    print("integration status:", stage_a2_process_time_status()["integration_status"])
