# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.
# Intuitive_Theory.md section 0: measurability axiom.
# Free substrate is not directly measurable; measurable physics begins when an
# oscillon-tail trace makes the substrate locally distinguishable.  p01 uses
# Y=1, B=delta as the effective normalized background representation of that
# axiom and tests consistency around it.

import sys
from enum import Enum
from typing import NamedTuple

import sympy as sp


class NoGhostResult(NamedTuple):
    K_PhiPhi_FLRW: sp.Expr
    K_pipi_FLRW: sp.Expr
    K_PhiPhi_Minkowski: sp.Expr
    K_pipi_Minkowski: sp.Expr


class SphericalInvariants(NamedTuple):
    r: sp.Expr
    theta: sp.Expr
    A: sp.Expr
    B: sp.Expr
    C: sp.Expr
    f: sp.Expr
    Y: sp.Expr
    I1: sp.Expr
    I2: sp.Expr
    I3: sp.Expr


class BackgroundName(str, Enum):
    MINKOWSKI = "minkowski"
    FLRW = "flrw"
    BIANCHI_I = "bianchi_i"
    SCHWARZSCHILD = "schwarzschild"


class BackgroundData(NamedTuple):
    coords: list[sp.Expr]
    g_cov_diag: list[sp.Expr]
    g_inv_diag: list[sp.Expr]
    sqrt_minus_g: sp.Expr
    fields: list[sp.Expr]


P01_MAIN_SECTIONS = {
    "base",
    "spherical",
    "moduli",
    "stress",
    "horndeski",
    "hyperbolicity",
    "eft",
    "lorentz",
    "old",
    "audit",
}


def _requested_main_sections():
    """Keep import and default script execution light; run demos by section name."""
    args = {arg.lower() for arg in sys.argv[1:]}
    if "all" in args:
        return P01_MAIN_SECTIONS
    return args & P01_MAIN_SECTIONS


def _should_run_main_section(section_name):
    return section_name in _requested_main_sections()

def init_variables():
    # ფაზური ინვარიანტი
    Y = sp.Symbol('Y', real=True)
    
    # ელასტიური ინვარიანტები
    I1 = sp.Symbol('I1', real=True)
    I2 = sp.Symbol('I2', real=True)
    I3 = sp.Symbol('I3', real=True)
    
    return Y, I1, I2, I3

def get_lagrangian(Y, I1, I2, I3):
    # განვსაზღვროთ ზოგადი ლაგრანჟიანი როგორც ფუნქცია
    P = sp.Function('P')(Y, I1, I2, I3)
    return P

def get_polynomial_lagrangian(Y, I1, I2, I3):
    # მინიმალური პოლინომიალური ფორმა
    c_Y = sp.Symbol('c_Y', real=True)
    c_Y2 = sp.Symbol('c_Y2', real=True)
    c_I1 = sp.Symbol('c_I1', real=True)
    c_I1sq = sp.Symbol('c_I1sq', real=True) # I1-ის კვადრატული წევრი
    c_I2 = sp.Symbol('c_I2', real=True)
    c_I3 = sp.Symbol('c_I3', real=True)
    
    c_YI1 = sp.Symbol('c_YI1', real=True) # შერეული წევრი (ფაზა-ელასტიურობა)
    L_poly = c_Y * Y + c_Y2 * Y**2 + c_I1 * I1 + c_I1sq * I1**2 + c_I2 * I2 + c_I3 * I3 + c_YI1 * Y * I1
    return L_poly

def get_energy_density(L, Y):
    # ენერგიის სიმკვრივის ზოგადი ფორმულა ფაზური ცვლადის მიმართ
    rho = 2 * Y * sp.diff(L, Y) - L
    return sp.simplify(rho)

def calculate_stress_tensor(L, metric_inverse_diagonal):
    """
    Diagonal mixed stress tensor helper used by consolidated FLRW modules.

    For independent diagonal inverse-metric variables q_i = g^{ii}, the local
    convention is T^i_i = 2*q_i*dL/dq_i - L.
    """
    return [
        sp.simplify(2 * q_i * sp.diff(L, q_i) - L)
        for q_i in metric_inverse_diagonal
    ]

def analyze_no_ghost() -> NoGhostResult:
    a = sp.Symbol('a', real=True, positive=True) # FLRW scale factor
    dPhi_dot = sp.Symbol('dPhi_dot', real=True)
    pi1_dot, pi2_dot, pi3_dot = sp.symbols('pi1_dot pi2_dot pi3_dot', real=True)
    pi_dot_sq = pi1_dot**2 + pi2_dot**2 + pi3_dot**2
    
    # ფლუქტუაციები
    # შენიშვნა: Y=1, B=δ არის Intuitive_Theory.md §0-ის გაზომვადობის
    # აქსიომის ეფექტური ნორმალიზებული ფონი. p01 ამ აქსიომას არ ამტკიცებს;
    # ის ამოწმებს, თავსებადია თუ არა ეფექტური პოლინომი ამ ფონის გარშემო.
    Y_pert = 1 + 2*dPhi_dot + dPhi_dot**2
    I1_pert = 3/a**2 - pi_dot_sq
    I2_pert = 3/a**4 - 2/a**2 * pi_dot_sq
    I3_pert = 1/a**6 - 1/a**4 * pi_dot_sq
    
    Y_s, I1_s, I2_s, I3_s = sp.symbols('Y I1 I2 I3', real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    L_pert = L_poly.subs({Y_s: Y_pert, I1_s: I1_pert, I2_s: I2_pert, I3_s: I3_pert})
    
    # კინეტიკური მატრიცის დიაგონალური წევრები (მეორე რიგის წარმოებულები)
    K_PhiPhi = sp.simplify(sp.diff(L_pert, dPhi_dot, 2) / 2)
    K_pipi = sp.simplify(sp.diff(L_pert, pi1_dot, 2) / 2)
    
    # შეფასება ფონის წერტილში (ფლუქტუაციები ნულზე)
    bg_subs = {dPhi_dot: 0, pi1_dot: 0, pi2_dot: 0, pi3_dot: 0}
    K_PhiPhi = sp.simplify(K_PhiPhi.subs(bg_subs))
    K_pipi = sp.simplify(K_pipi.subs(bg_subs))
    
    K_PhiPhi_Mink = sp.simplify(K_PhiPhi.subs(a, 1))
    K_pipi_Mink = sp.simplify(K_pipi.subs(a, 1))
    
    return NoGhostResult(K_PhiPhi, K_pipi, K_PhiPhi_Mink, K_pipi_Mink)

def analyze_lorentz_constrained_stability():
    """
    ამოწმებს No-Ghost პირობებს Lorentz (T_01=0) და PPN (gamma=1) 
    კონსტრეინტების ჩასმის შემდეგ.
    """
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols('c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1', real=True)
    K_PhiPhi, K_pipi, K_PhiPhi_Mink, K_pipi_Mink = analyze_no_ghost()
    
    # Lorentz + PPN გვაძლევს ორ პირობას. ამ ეტაპზე ეს არის exact tuning
    # condition; დამცავი სიმეტრია/RG-stability ჯერ ცალკე დასამტკიცებელია.
    # 1. c_Y2 = c_I1sq
    # 2. c_I1 = c_Y - 4*c_Y2 + 2*c_YI1 - 2*c_I2 - c_I3
    subs_dict = {
        c_I1sq: c_Y2,
        c_I1: c_Y - 4*c_Y2 + 2*c_YI1 - 2*c_I2 - c_I3
    }
    
    K_Phi_constr = sp.simplify(K_PhiPhi_Mink.subs(subs_dict))
    K_pi_constr = sp.simplify(K_pipi_Mink.subs(subs_dict))
    
    return K_Phi_constr, K_pi_constr


def foundational_axiom_bridge():
    """
    Record the foundational axiom bridge used by p01.

    Intuitive_Theory.md section 0 states the measurability axiom: free
    substrate is not directly measurable; measurable physics begins only when
    an oscillon-tail trace creates local distinguishability.  p01 does not
    derive this primitive.  It represents its unexcited effective background
    by Y=1, B=delta and tests consistency around that normalized branch.
    """
    return {
        "source": "Intuitive_Theory.md section 0",
        "foundational_axiom": (
            "Free substrate is not directly measurable; it becomes measurable "
            "only through oscillon-tail-induced distinguishable nodal, phase, "
            "or pressure traces."
        ),
        "formal_role_in_p01": (
            "Y=1, B^{AB}=delta^{AB} is the normalized local homogeneous "
            "effective background used for stress, no-ghost and hyperbolicity tests."
        ),
        "not_claimed_here": (
            "p01 does not prove or derive the foundational axiom from the "
            "polynomial core."
        ),
        "p01_obligation": (
            "declare the axiom boundary, keep the background normalized, and "
            "check algebraic zero-stress plus perturbative consistency conditions."
        ),
    }


def normalized_substrate_background_audit():
    """
    Make the Y=1, B=delta background status explicit.

    This is the effective representation of the foundational measurability
    axiom, not a dynamically derived minimum.  The function computes the
    algebraic zero-stress constraints that must hold when p01 is expanded
    around the normalized background.
    """
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        'c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1', real=True
    )

    rho0 = sp.simplify(
        c_Y + 3 * c_Y2 + 3 * c_YI1
        - 3 * c_I1 - 9 * c_I1sq - 3 * c_I2 - c_I3
    )
    p0 = sp.simplify(
        c_Y + c_Y2 + c_YI1
        + c_I1 - 3 * c_I1sq - c_I2 - c_I3
    )
    solutions = sp.solve([rho0, p0], [c_Y, c_I1], dict=True)
    solution = solutions[0] if solutions else {}

    K_phi = sp.simplify(c_Y + 6 * c_Y2 + 3 * c_YI1)
    K_pi = sp.simplify(-c_I1 - 6 * c_I1sq - 2 * c_I2 - c_I3 - c_YI1)
    K_phi_on_vacuum = sp.simplify(K_phi.subs(solution)) if solution else None
    K_pi_on_vacuum = sp.simplify(K_pi.subs(solution)) if solution else None

    return {
        "background": "Y=1, B^{AB}=delta^{AB}",
        "status": "FOUNDATIONAL_AXIOM_EFFECTIVE_BACKGROUND",
        "axiom_bridge": foundational_axiom_bridge(),
        "axiom_boundary": (
            "This is a declared primitive of the theory stack, not a p01 proof "
            "target and not a polynomial minimum claim."
        ),
        "p01_consistency_obligation": (
            "p01 must check zero-stress, no-ghost, gradient/mixed-mode and "
            "stress/Lorentz-background consistency around this normalized input."
        ),
        "rho0": rho0,
        "p0": p0,
        "zero_stress_solution_for_cY_cI1": solution,
        "K_phi_on_zero_stress_branch": K_phi_on_vacuum,
        "K_pi_on_zero_stress_branch": K_pi_on_vacuum,
        "next_test_target": (
            "keep strengthening the effective consistency checks around the "
            "axiomatic normalized background"
        ),
    }


def homogeneous_vacuum_derivation_attempt():
    """
    Diagnostic guardrail against deriving Y=1, B=b*delta from the polynomial core.

    Result:
    - Constant-gradient field equations do not select y=b=1 by themselves;
      they are automatically solved when the currents are constant.
    - Zero stress gives an algebraic vacuum branch.
    - If we also demand phase stationarity L_y=0 at y=b=1, the solid kinetic
      coefficient K_pi becomes zero on that branch.

    This is a guardrail against over-claiming, not a required closure route.
    In the current theory stack, y=b=1 is the effective normalized background
    of the foundational measurability axiom, not a minimum derived inside p01.
    """
    y, b = sp.symbols("y b", positive=True)
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )

    L_iso = sp.simplify(
        c_Y * y
        + c_Y2 * y**2
        + 3 * c_I1 * b
        + (9 * c_I1sq + 3 * c_I2) * b**2
        + c_I3 * b**3
        + 3 * c_YI1 * y * b
    )
    rho_iso = sp.simplify(2 * y * sp.diff(L_iso, y) - L_iso)
    p_iso = sp.simplify(L_iso - sp.Rational(2, 3) * b * sp.diff(L_iso, b))

    at_vac = {y: 1, b: 1}
    rho0 = sp.simplify(rho_iso.subs(at_vac))
    p0 = sp.simplify(p_iso.subs(at_vac))
    L_y0 = sp.simplify(sp.diff(L_iso, y).subs(at_vac))
    L_b0 = sp.simplify(sp.diff(L_iso, b).subs(at_vac))

    zero_stress_solution = sp.solve([rho0, p0], [c_Y, c_I1], dict=True)[0]
    phase_stationary_solution = sp.solve(
        [rho0, p0, L_y0],
        [c_Y, c_I1, c_I2],
        dict=True,
    )[0]
    full_modulus_stationary_solution = sp.solve(
        [rho0, p0, L_y0, L_b0],
        [c_Y, c_I1, c_I2, c_I3],
        dict=True,
    )
    full_modulus_stationary_solution = (
        full_modulus_stationary_solution[0]
        if full_modulus_stationary_solution
        else {}
    )

    K_phi = sp.simplify(c_Y + 6 * c_Y2 + 3 * c_YI1)
    K_pi = sp.simplify(-c_I1 - 6 * c_I1sq - 2 * c_I2 - c_I3 - c_YI1)

    return {
        "isotropic_L_y_b": L_iso,
        "rho_y_b": rho_iso,
        "p_y_b": p_iso,
        "rho0": rho0,
        "p0": p0,
        "phase_stationarity_L_y0": L_y0,
        "solid_scale_stationarity_L_b0": L_b0,
        "zero_stress_solution": zero_stress_solution,
        "K_phi_zero_stress": sp.simplify(K_phi.subs(zero_stress_solution)),
        "K_pi_zero_stress": sp.simplify(K_pi.subs(zero_stress_solution)),
        "phase_stationary_solution": phase_stationary_solution,
        "K_phi_phase_stationary": sp.simplify(K_phi.subs(phase_stationary_solution)),
        "K_pi_phase_stationary": sp.simplify(K_pi.subs(phase_stationary_solution)),
        "full_modulus_stationary_solution": full_modulus_stationary_solution,
        "K_phi_full_modulus_stationary": sp.simplify(
            K_phi.subs(full_modulus_stationary_solution)
        ) if full_modulus_stationary_solution else None,
        "K_pi_full_modulus_stationary": sp.simplify(
            K_pi.subs(full_modulus_stationary_solution)
        ) if full_modulus_stationary_solution else None,
        "obstruction": (
            "zero-stress plus L_y=0 gives K_pi=0; therefore the current "
            "polynomial core must not claim to derive the healthy Y=1, "
            "B=delta background as a simple homogeneous modulus extremum."
        ),
        "guardrail": (
            "do not add homogeneous phase-stationarity as a p01 axiom-closure "
            "requirement; keep Y=1, B=delta as axiomatic normalized input and "
            "test its effective consequences"
        ),
    }

def analyze_sound_speeds(solve_roots=False):
    """
    Minkowski ფონზე ვითვლით ტრანსვერსულ და შერეულ (ფაზა+გრძივი) ხმის სიჩქარეებს.
    """
    # Minkowski ფონზე
    dPhi_dot, dPhi_z = sp.symbols('dPhi_dot dPhi_z', real=True)
    pi1_dot, pi2_dot, pi3_dot = sp.symbols('pi1_dot pi2_dot pi3_dot', real=True)
    pi1_z, pi2_z, pi3_z = sp.symbols('pi1_z pi2_z pi3_z', real=True)
    
    Y_pert = 1 + 2*dPhi_dot + dPhi_dot**2 - dPhi_z**2
    
    B11 = 1 - pi1_dot**2 + pi1_z**2
    B22 = 1 - pi2_dot**2 + pi2_z**2
    B33 = 1 - pi3_dot**2 + 2*pi3_z + pi3_z**2
    # შენიშვნა: 2*pi3_z მოდის ფონის phi^3=z რუკიდან და არის გრძივი პერტურბაციის ხაზოვანი წევრი.
    B12 = -pi1_dot*pi2_dot + pi1_z*pi2_z
    B13 = -pi1_dot*pi3_dot + pi1_z*(1 + pi3_z)
    B23 = -pi2_dot*pi3_dot + pi2_z*(1 + pi3_z)
    
    B = sp.Matrix([[B11, B12, B13],
                   [B12, B22, B23],
                   [B13, B23, B33]])
    
    I1_pert = sp.simplify(B.trace())
    I2_pert = sp.simplify(sp.Rational(1,2) * (I1_pert**2 - (B*B).trace()))
    I3_pert = sp.simplify(B.det())
    
    Y_s, I1_s, I2_s, I3_s = sp.symbols('Y I1 I2 I3', real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    
    eps = sp.Symbol('eps', real=True)
    subs_dict = {
        dPhi_dot: eps*dPhi_dot, dPhi_z: eps*dPhi_z,
        pi1_dot: eps*pi1_dot, pi2_dot: eps*pi2_dot, pi3_dot: eps*pi3_dot,
        pi1_z: eps*pi1_z, pi2_z: eps*pi2_z, pi3_z: eps*pi3_z
    }
    
    L_eval = L_poly.subs({Y_s: Y_pert, I1_s: I1_pert, I2_s: I2_pert, I3_s: I3_pert}).subs(subs_dict)
    L_O2 = sp.simplify(sp.series(L_eval, eps, 0, 3).coeff(eps, 2))
    
    # ტრანსვერსული მოდი (pi1) - რჩება ცალკე
    K_T = sp.simplify(L_O2.coeff(pi1_dot**2))
    G_T = sp.simplify(-L_O2.coeff(pi1_z**2))
    cs2_T = sp.simplify(G_T / K_T)
    
    # შერეული 2x2 მატრიცა ფაზისა და გრძივი მოდისთვის: {dPhi, pi3}
    A = sp.simplify(L_O2.coeff(dPhi_dot**2))
    B_pi3 = sp.simplify(L_O2.coeff(pi3_dot**2))
    C = sp.simplify(L_O2.coeff(dPhi_z**2))
    D = sp.simplify(L_O2.coeff(pi3_z**2))
    M_mix = sp.simplify(L_O2.coeff(dPhi_dot * pi3_z) + L_O2.coeff(pi3_dot * dPhi_z))
    
    # დეტერმინანტი det(K*cs² + G) = 0 მოგვცემს კვადრატულ განტოლებას cs2-სთვის.
    # M_mix=0 ლიმიტში ფესვები უნდა იყოს -C/A და -D/B.
    cs2 = sp.Symbol('cs2', real=True)
    eq_cs2 = sp.factor(sp.simplify((A*cs2 + C) * (B_pi3*cs2 + D) - M_mix**2 * cs2))
    
    # sp.solve აბრუნებს უზარმაზარ ფესვებს, ამიტომ ვაბრუნებთ მატრიცის კოეფიციენტებს და დამახასიათებელ განტოლებას
    coeffs = {'A': A, 'B_pi3': B_pi3, 'C': C, 'D': D, 'M_mix': M_mix}
    
    # სრული სიმბოლური ფესვები მძიმეა; საჭიროებისას გაიშვება ცალკე.
    cs2_roots = sp.solve(eq_cs2, cs2) if solve_roots else None
    
    return cs2_T, eq_cs2, coeffs, cs2_roots


def sound_speed_decoupled_limit_check():
    """Check M_mix=0 limit: roots are -C/A and -D/B."""
    _cs2_T, eq_cs2, coeffs, _roots = analyze_sound_speeds(solve_roots=False)
    cs2 = sp.Symbol("cs2", real=True)
    c_YI1 = sp.Symbol("c_YI1", real=True)
    eq_decoupled = sp.simplify(eq_cs2.subs(c_YI1, 0))
    A = sp.simplify(coeffs["A"].subs(c_YI1, 0))
    B_long = sp.simplify(coeffs["B_pi3"].subs(c_YI1, 0))
    C = sp.simplify(coeffs["C"].subs(c_YI1, 0))
    D = sp.simplify(coeffs["D"].subs(c_YI1, 0))
    return {
        "phase_root_residual": sp.simplify(eq_decoupled.subs(cs2, -C / A)),
        "longitudinal_root_residual": sp.simplify(eq_decoupled.subs(cs2, -D / B_long)),
    }

if __name__ == "__main__" and _should_run_main_section("base"):
    Y, I1, I2, I3 = init_variables()
    L_poly = get_polynomial_lagrangian(Y, I1, I2, I3)
    rho_poly = get_energy_density(L_poly, Y)
    
    print("პოლინომიალური ლაგრანჟიანი:", L_poly)
    print("ენერგიის სიმკვრივე:", rho_poly)

    K_PhiPhi, K_pipi, K_PhiPhi_Mink, K_pipi_Mink = analyze_no_ghost()
    print("\n--- No-ghost პირობები (კინეტიკური მატრიცის დიაგონალი ფონზე) ---")
    print("Minkowski ფონი:")
    print("K_PhiPhi > 0 =>", K_PhiPhi_Mink, "> 0")
    print("K_pipi > 0 =>", K_pipi_Mink, "> 0")
    print("\nFLRW ფონი:")
    print("K_PhiPhi > 0 =>", K_PhiPhi, "> 0")
    print("K_pipi > 0 =>", K_pipi, "> 0")

    print("\n--- Foundational axiom effective-background audit ---")
    vacuum_audit = normalized_substrate_background_audit()
    print("ფონი:", vacuum_audit["background"])
    print("სტატუსი:", vacuum_audit["status"])
    print("აქსიომის ხიდი:", vacuum_audit["axiom_bridge"]["formal_role_in_p01"])
    print("p01 ვალდებულება:", vacuum_audit["p01_consistency_obligation"])
    print("rho0 =", vacuum_audit["rho0"])
    print("p0   =", vacuum_audit["p0"])
    print("zero-stress solution:", vacuum_audit["zero_stress_solution_for_cY_cI1"])
    print("K_Phi zero-stress branch:", vacuum_audit["K_phi_on_zero_stress_branch"])
    print("K_pi zero-stress branch:", vacuum_audit["K_pi_on_zero_stress_branch"])
    print("შემდეგი ტესტ-სამიზნე:", vacuum_audit["next_test_target"])

    print("\n--- Polynomial-minimum no-go diagnostic ---")
    derivation_attempt = homogeneous_vacuum_derivation_attempt()
    print("phase stationarity L_y0 =", derivation_attempt["phase_stationarity_L_y0"])
    print("solid scale stationarity L_b0 =", derivation_attempt["solid_scale_stationarity_L_b0"])
    print("phase-stationary solution:", derivation_attempt["phase_stationary_solution"])
    print("K_Phi on phase-stationary branch:", derivation_attempt["K_phi_phase_stationary"])
    print("K_pi on phase-stationary branch:", derivation_attempt["K_pi_phase_stationary"])
    print("obstruction:", derivation_attempt["obstruction"])
    print("guardrail:", derivation_attempt["guardrail"])

    cs2_T, eq_cs2, coeffs, cs2_roots = analyze_sound_speeds(solve_roots=False)
    print("\n--- Sound Speeds (c_s^2) ---")
    print("Transverse Elastic Mode (pi_T):", cs2_T)
    print("\nMixed Phase + Longitudinal Mode 2x2 System (dPhi, pi_3):")
    print("Characteristic Equation for cs^2:", eq_cs2, "= 0")
    print("Matrix Coefficients (K and G parts):")
    print(f"  K_PhiPhi (A) = {coeffs['A']}")
    print(f"  K_L (B_pi3) = {coeffs['B_pi3']}")
    print(f"  G_PhiPhi (C) = {-coeffs['C']}")
    print(f"  G_L (D) = {-coeffs['D']}")
    print(f"  Mixing term (M_mix) = {coeffs['M_mix']}")
    print("Decoupled-limit residuals:", sound_speed_decoupled_limit_check())
    print("\nსაკუთრივი მნიშვნელობები (Eigenmode Speeds c_s^2):")
    if cs2_roots:
        print("Root 1:", cs2_roots[0])
        print("Root 2:", cs2_roots[1])
    else:
        print("სიმბოლური ფესვები default რეჟიმში არ იხსნება; გამოიყენება პოლინომი და decoupled ლიმიტი.")

    K_Phi_c, K_pi_c = analyze_lorentz_constrained_stability()
    print("\n--- აქსიომატურად ნორმალიზებული ფონის სტაბილურობა ---")
    print("კონსტრეინტების (c_Y2 = c_I1sq და PPN) ჩასმის შემდეგ No-Ghost პირობები:")
    print(f"K_PhiPhi > 0 => {K_Phi_c} > 0")
    print(f"K_pipi > 0   => {K_pi_c} > 0")
    print("დასკვნა: შერჩეულ ფონზე ორი დიაგონალური no-ghost პირობა ერთდროულად სრულდება, თუ c_Y2 > 0 და")
    print("-6*c_Y2 < (c_Y + 3*c_YI1) < -2*c_Y2.")
    print("ეს არის აუცილებელი ფანჯარა; gradient/cross-mode/eigenmode სტაბილურობა ცალკე ღიაა.")

    print("\n--- აგენტთა საბჭოს შენიშვნები ---")
    print("- ფონი (Y=1, B=δ) არის გაზომვადობის აქსიომის ნორმალიზებული ეფექტური ფონი;")
    print("  p01 აქსიომას არ ამტკიცებს, მხოლოდ ამოწმებს მის გარშემო ეფექტურ პირობებს.")
    print("- 2*pi3_z მოდის ფონის phi^3=z რუკიდან და არის გრძივი პერტურბაციის ხაზოვანი წევრი.")
    print("- ჯვარედინი შერევა დეტალურად აისახა 2x2 მატრიცის დეტერმინანტით.")
    print("- cs^2 განტოლება ჩაწერილია det(K*cs^2 + G)=0 ფორმით;")
    print("  M_mix=0 ლიმიტში ფესვები არის -C/A და -D/B.")

# ===================== CONSOLIDATED PHASE SECTIONS =====================


# ===================== merged from p01_core.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 1 (tensor): სუპერსოლიდის სტრესი სფერული ანზაცით
================================================================================

რეფერენცია: NOTATION.md, p01_core.py

ეს ფაილი იყენებს NOTATION.md-ის აქტიურ კონვენციას:
- სიგნატურა (+---)
- B^{AB} = -g^{mu nu} * d_mu phi^A * d_nu phi^B
- T_{mu nu} = 2*dL/dg^{mu nu} - g_{mu nu}*L (unique symmetric variables)
- off-diagonal: ფაქტორი 1 (არა 2), რადგან g^{mn}=g^{nm}

phase22 ფარდდება Bianchi/Noether იდენტობას სამ ფონზე (Minkowski, FLRW,
Schwarzschild) stress-tensor coefficient sanity-check-ით. ეს ფაილი იყენებს იმავე კონვენციებს
სფერული ანზაცის Δp გენერაციისთვის (MOND-ის ფორმულა §4-ში).

სფერული სტატიკური ანზაცი:
    ds^2 = B(r)*dt^2 - A(r)*dr^2 - C(r)*dOmega^2
ფაზური ველი: Phi = t (სტატიკური)
ელასტიური ველი: phi^1 = f(r) (რადიალური დეფორმაცია), phi^2 = theta, phi^3 = phi

შედეგი: rho, p_rad, p_tan, Δp = p_tan - p_rad
"""

import sympy as sp
# merged import removed: from p01_core import get_polynomial_lagrangian


def get_spherical_invariants() -> SphericalInvariants:
    """
    სფერული სტატიკური ანზაცის ინვარიანტები NOTATION.md-ის კონვენციით.

    g_inv დიაგონალური: (1/B, -1/A, -1/C, -1/(C*sin^2(theta)))
    f(r) — SO(3)-კოვარიანტული რადიალური დეფორმაცია phi^A = f(r) n^A(theta, phi)
    """
    r = sp.Symbol("r", real=True, positive=True)
    theta = sp.Symbol("theta", real=True, positive=True)
    A = sp.Function("A")(r)
    B = sp.Function("B")(r)
    C = sp.Function("C")(r)
    f = sp.Function("f")(r)
    f_prime = sp.diff(f, r)

    # ფაზური ველი — სტატიკური Phi = t
    # ე.ი. d_0 Phi = 1, d_i Phi = 0
    # Y = g^{00} * (d_0 Phi)^2 = 1/B
    Y = 1 / B

    # ელასტიური ველი — phi^A=f(r)n^A. მისი eigenvalue-ებია:
    # lambda_r=f'^2/A და lambda_t=f^2/C (ორჯერ). ეს არჩევანი ხსნის
    # theta-დამოკიდებულ კოორდინატულ არტეფაქტს phi^A=(f,theta,phi) რუკიდან.
    lambda_r = f_prime**2 / A
    lambda_t = f**2 / C

    I1 = sp.simplify(lambda_r + 2 * lambda_t)
    I2 = sp.simplify(2 * lambda_r * lambda_t + lambda_t**2)
    I3 = sp.simplify(lambda_r * lambda_t**2)

    return SphericalInvariants(r, theta, A, B, C, f, Y, I1, I2, I3)


def get_lagrangian_spherical():
    """L = poly(Y, I1, I2, I3) სფერული ანზაცის ცვლადებში."""
    r, theta, A, B, C, f, Y, I1, I2, I3 = get_spherical_invariants()
    Y_s, I1_s, I2_s, I3_s = sp.symbols("Y I1 I2 I3", real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    L = sp.simplify(L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3}))
    return r, theta, A, B, C, f, L, (Y, I1, I2, I3)


def get_stress_tensor():
    """
    T_{mu nu} გენერალური q-ცვლადებით (NOTATION.md-ის კონვენცია), შემდეგ
    სფერული ჩასმა. T-ის სქელეტი — დიაგონალური მხოლოდ (off-diagonal სფერულ
    სტატიკურ ანზაცზე ნულია).
    """
    r, theta, A, B, C, f, _, _ = get_lagrangian_spherical()
    f_prime = sp.diff(f, r)

    # q = (g^00, g^11, g^22, g^33) სიმბოლური
    q = sp.symbols("q0 q1 q2 q3", real=True, nonzero=True)
    # SO(3)-კოვარიანტული phi^A=f(r)n^A რუკა:
    # lambda_r=-q1*f'^2, lambda_theta=-q2*f^2,
    # lambda_phi=-q3*f^2*sin(theta)^2. სფერული ჩასმის შემდეგ
    # lambda_theta=lambda_phi=f^2/C.
    Y_sym = q[0]
    lambda_r = -q[1] * f_prime**2
    lambda_theta = -q[2] * f**2
    lambda_phi = -q[3] * f**2 * sp.sin(theta)**2
    I1_sym = lambda_r + lambda_theta + lambda_phi
    I2_sym = (
        lambda_r * lambda_theta
        + lambda_r * lambda_phi
        + lambda_theta * lambda_phi
    )
    I3_sym = lambda_r * lambda_theta * lambda_phi

    Y_s, I1_s, I2_s, I3_s = sp.symbols("Ys I1s I2s I3s", real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    L_sym = L_poly.subs({Y_s: Y_sym, I1_s: I1_sym, I2_s: I2_sym, I3_s: I3_sym})

    # T_{mu mu} = 2 * dL/dq[mu] - L/q[mu]    (q[mu] = g^{mu mu}, g_{mu mu} = 1/q[mu])
    T_cov_general = [
        2 * sp.diff(L_sym, q[mu]) - L_sym / q[mu]
        for mu in range(4)
    ]

    # სფერული ჩასმა
    subs_sph = {
        q[0]: 1 / B,
        q[1]: -1 / A,
        q[2]: -1 / C,
        q[3]: -1 / (C * sp.sin(theta)**2),
    }

    T_cov = [sp.simplify(expr.subs(subs_sph)) for expr in T_cov_general]

    return r, theta, A, B, C, f, T_cov


def get_pressures():
    """
    rho, p_rad, p_tan, Δp სფერული ანზაციდან NOTATION.md-ის კონვენციით.

    rho = T^{0}_{0} = g^{00} * T_{00} = (1/B) * T_{00}
    p_rad = -T^{1}_{1} = -g^{11} * T_{11} = -(-1/A) * T_{11} = T_{11}/A
    p_tan = -T^{2}_{2} = -g^{22} * T_{22} = T_{22}/C
    Δp = p_tan - p_rad
    """
    r, theta, A, B, C, f, T_cov = get_stress_tensor()

    rho = sp.simplify(T_cov[0] / B)
    p_rad = sp.simplify(T_cov[1] / A)
    p_tan = sp.simplify(T_cov[2] / C)
    delta_p = sp.simplify(p_tan - p_rad)

    return r, theta, A, B, C, f, rho, p_rad, p_tan, delta_p


if __name__ == "__main__" and _should_run_main_section("spherical"):
    print("=" * 72)
    print("PHASE 1 (tensor): სუპერსოლიდის სტრესი სფერული ანზაცით")
    print("რეფერენცია: NOTATION.md, phase22")
    print("=" * 72)

    r, theta, A, B, C, f, rho, p_rad, p_tan, delta_p = get_pressures()

    c_I1, c_I1sq, c_YI1, c_I2, c_I3 = sp.symbols(
        "c_I1 c_I1sq c_YI1 c_I2 c_I3", real=True
    )

    print("\n1. ენერგიის სიმკვრივე (rho):")
    print(sp.simplify(rho))

    print("\n2. გრძივი წნევა (p_rad):")
    print(sp.simplify(p_rad))

    print("\n3. განივი წნევა (p_tan):")
    print(sp.simplify(p_tan))

    print("\n4. ანიზოტროპია (Δp = p_tan - p_rad):")
    delta_p_expanded = sp.expand(delta_p)
    print(sp.collect(delta_p_expanded, [c_I1, c_I1sq, c_YI1, c_I2, c_I3]))

    print("\n5. სტატუსი:")
    print("  - კონვენცია: NOTATION.md-ის აქტიური ფორმა (T_mn = 2*dL/dg^mn - g_mn*L)")
    print("  - სიგნატურა: (+---)")
    print("  - სოლიდის რუკა: SO(3)-კოვარიანტული phi^A=f(r)n^A, theta-არტეფაქტის გარეშე")
    print("  - phase22-ის Bianchi/Noether იდენტობა იყენებს იმავე კონვენციას")
    print("  - Δp გენერირდება f'(r), A, B, C-ის ფუნქციად — ეს არის სტატიისთვის გამოსატანი სამუშაო შედეგი")


# ===================== merged from p01_core.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
ლოკალური მოდულუსების შემოწმების ესკიზი.

შენიშვნა:
- ეს ფაილი არ არის ტალღების სრული მტკიცება და არ შეიცავს დიაგონალიზებულ eigenmode-ებს.
- სრული ტალღებისთვის (2x2 კინეტიკური/გრადიენტული მატრიცა) გამოიყენება p01_core.py.
- გრავიტაციული ტალღებისთვის გამოიყენება p04_gw.py (TT მეტრიკის ექსპანსია).
- H_YY არ უნდა ჩაითვალოს გრავიტაციული ტალღების სიხისტედ; c_T მოდის აინშტაინ-ჰილბერტის სექტორიდან.
"""
import sympy as sp
# merged import removed: from p01_core import init_variables, get_polynomial_lagrangian

def get_local_moduli():
    Y, I1, I2, I3 = init_variables()
    L = get_polynomial_lagrangian(Y, I1, I2, I3)
    
    # ლოკალური მოდულუსებისთვის ვითვლით ლაგრანჟიანის მეორე რიგის წარმოებულებს (ჰესიანს)
    hessian_YY = sp.diff(L, Y, 2)
    hessian_I1I1 = sp.diff(L, I1, 2)
    hessian_YI1 = sp.diff(L, Y, I1)
    
    # ფონზე შეფასება (Minkowski: Y=1, I1=3, I2=3, I3=1)
    bg_subs = {Y: 1, I1: 3, I2: 3, I3: 1}
    
    H_YY_bg = sp.simplify(hessian_YY.subs(bg_subs))
    H_I1I1_bg = sp.simplify(hessian_I1I1.subs(bg_subs))
    H_YI1_bg = sp.simplify(hessian_YI1.subs(bg_subs))
    
    # სრული 2x2 ჰესიანის დეტერმინანტი
    H_det_bg = sp.simplify(H_YY_bg * H_I1I1_bg - H_YI1_bg**2)
    
    return H_YY_bg, H_I1I1_bg, H_YI1_bg, H_det_bg

if __name__ == "__main__" and _should_run_main_section("moduli"):
    h_YY, h_I1I1, h_YI1, h_det = get_local_moduli()
    print("ფაზური ლოკალური მოდულუსი ფონზე (H_YY):", h_YY)
    print("ელასტიური ლოკალური მოდულუსი ფონზე (H_I1I1):", h_I1I1)
    print("შერეული კავშირის მოდულუსი ფონზე (H_YI1):", h_YI1)
    
    print("\n--- ლოკალური სტაბილურობის პირობები (მოდულუსების დადებითობა) ---")
    print(f"დიაგონალური ფაზური სტაბილურობა: H_YY > 0  =>  {h_YY} > 0")
    print(f"დიაგონალური ელასტიური სტაბილურობა: H_I1I1 > 0  =>  {h_I1I1} > 0")
    print(f"სრული 2x2 ჰესიანის დადებითობა (Det > 0): {h_det} > 0")
    print("\nეს აჩვენებს მხოლოდ ლოკალური დიაგონალური და შერეული მოდულუსების დადებითობის პირობებს;")
    print("სრული ტალღური სტაბილურობა და eigenmode-ები მოწმდება p01_core.py-ში.")


# ===================== merged from p01_core.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 22 (v3.1): სრული ენერგია-იმპულსის ტენზორი — Carter-Karlovini ფორმალიზმი
================================================================================

სტატუსი:
ეტაპი I-ის ფესვის გასწორების ფაილი. აქ ფიქსირდება:
- T_{mu nu} = 2*dL/dg^{mu nu} - g_{mu nu}*L კონვენცია
- I_1, I_2, I_3 Carter-Karlovini-ის ჯაჭვური წესით B^{AB}-დან
- nabla_mu T^{mu nu} = sum_A E_A * partial^nu psi_A Noether/Bianchi იდენტობა
  diagonal Minkowski/FLRW/Bianchi I ფონებზე სრული L-ით
- Schwarzschild diagonal ფონზე იგივე იდენტობის reduced Y+I1 smoke-test
- off-diagonal ვარიაციის ფაქტორის ცალკე smoke-test
- stress-tensor coefficient sanity check, რომ Noether იდენტობა მცდარი
  T-კოეფიციენტის შემთხვევაში არანულოვან residual-ს იძლევა

შენიშვნა მკითხველისთვის:
matter stress tensor-ის off-shell კოვარიანტული დივერგენცია generic-ად
ნული არ არის — ის ფიზიკურ ველთა Euler-Lagrange წყაროებს უდრის. ეს ფაილი
ამოწმებს ამ იდენტობას, არა "off-shell zero"-ს.
stress_tensor_coefficient_sanity_check() აჩვენებს, რომ T-ის არასწორი
ნორმალიზაცია ამ ფონზე არანულოვან residual-ს ტოვებს; ეს არის
კონვენციის sanity-check, არა მთელი თეორიის ფალსიფიკაცია.
"""

import sympy as sp
# merged import removed: from p01_core import get_polynomial_lagrangian


DIM = 4
NSOLID = 3


# ============================================================================
# გეომეტრიული ფუნქციები
# ============================================================================


def get_christoffel(g_cov, g_inv, coords):
    """Christoffel სიმბოლოები Gamma^lambda_{mu nu}."""
    Gamma = [sp.zeros(DIM, DIM) for _ in range(DIM)]
    half = sp.Rational(1, 2)
    for lam in range(DIM):
        for mu in range(DIM):
            for nu in range(DIM):
                term = 0
                for rho in range(DIM):
                    term += half * g_inv[lam, rho] * (
                        sp.diff(g_cov[rho, mu], coords[nu])
                        + sp.diff(g_cov[rho, nu], coords[mu])
                        - sp.diff(g_cov[mu, nu], coords[rho])
                    )
                Gamma[lam][mu, nu] = sp.simplify(term)
    return Gamma


def covariant_divergence_contra(T_contra, Gamma, coords, nu):
    """nabla_mu T^{mu nu}."""
    total = 0
    for mu in range(DIM):
        total += sp.diff(T_contra[mu, nu], coords[mu])
        for lam in range(DIM):
            total += Gamma[mu][mu, lam] * T_contra[lam, nu]
            total += Gamma[nu][mu, lam] * T_contra[mu, lam]
    return sp.simplify(total)


# ============================================================================
# ინვარიანტები და ლაგრანჟიანი
# ============================================================================


def build_invariants_from_metric(g_inv, field_grads):
    """
    field_grads[0] = d_mu Phi (ფაზური სკალარი)
    field_grads[1..3] = d_mu phi^A (ელასტიური ველები)
    """
    Y = sp.simplify(sum(
        g_inv[mu, nu] * field_grads[0][mu] * field_grads[0][nu]
        for mu in range(DIM) for nu in range(DIM)
    ))

    B = sp.zeros(NSOLID, NSOLID)
    for A in range(NSOLID):
        for Bidx in range(NSOLID):
            B[A, Bidx] = sp.simplify(sum(
                -g_inv[mu, nu] * field_grads[A + 1][mu] * field_grads[Bidx + 1][nu]
                for mu in range(DIM) for nu in range(DIM)
            ))

    I1 = sp.simplify(B.trace())
    I2 = sp.simplify(sp.Rational(1, 2) * (I1**2 - (B * B).trace()))
    I3 = sp.simplify(B.det())
    return Y, I1, I2, I3, B


def get_full_lagrangian(Y, I1, I2, I3):
    """p01_core-ის სრული პოლინომიური L."""
    Y_s, I1_s, I2_s, I3_s = sp.symbols("Y I1 I2 I3", real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    L = sp.simplify(L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3}))
    return L, L_poly, (Y_s, I1_s, I2_s, I3_s)


def get_test_lagrangian(Y, I1, I2, I3, mode="full"):
    """შესამოწმებელი L: სრული ან შემცირებული."""
    if mode == "full":
        return get_full_lagrangian(Y, I1, I2, I3)
    if mode == "reduced_y_i1":
        c_y, c_i1 = sp.symbols("c_y c_i1", real=True)
        L = c_y * Y + c_i1 * I1
        L_poly = c_y * sp.Symbol("Y", real=True) + c_i1 * sp.Symbol("I1", real=True)
        return L, L_poly, sp.symbols("Y I1 I2 I3", real=True)
    raise ValueError(f"unknown lagrangian mode: {mode}")


# ============================================================================
# Off-diagonal ვარიაციის ცდა (smoke-test)
# ============================================================================


def offdiag_variation_smoke_test():
    """
    შემოწმდება: symmetric g^{01} ცვლადის ვარიაცია ფაქტორი-2-ით
    არ უნდა გადიდდეს. შედარება — დამოუკიდებელი g^{01}, g^{10}.
    """
    q00, q11, q22, q33, q01 = sp.symbols("q00 q11 q22 q33 q01", real=True)
    r01, r10 = sp.symbols("r01 r10", real=True)
    gcov01 = sp.Symbol("gcov01", real=True)
    c_y, c_i1, c_yi1 = sp.symbols("c_y c_i1 c_yi1", real=True)
    v0, v1 = sp.symbols("v0 v1", real=True)
    e10, e11 = sp.symbols("e10 e11", real=True)

    field_grads = [
        [v0, v1, 0, 0],
        [e10, e11, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]

    g_inv_sym = sp.Matrix([
        [q00, q01, 0, 0],
        [q01, q11, 0, 0],
        [0, 0, q22, 0],
        [0, 0, 0, q33],
    ])
    Y_sym, I1_sym, _, _, _ = build_invariants_from_metric(g_inv_sym, field_grads)
    L_sym = c_y * Y_sym + c_i1 * I1_sym + c_yi1 * Y_sym * I1_sym
    T01_unique = sp.expand(sp.diff(L_sym, q01) - gcov01 * L_sym)
    T01_wrong = sp.expand(2 * sp.diff(L_sym, q01) - gcov01 * L_sym)

    g_inv_independent = sp.Matrix([
        [q00, r01, 0, 0],
        [r10, q11, 0, 0],
        [0, 0, q22, 0],
        [0, 0, 0, q33],
    ])
    Y_ind, I1_ind, _, _, _ = build_invariants_from_metric(g_inv_independent, field_grads)
    L_ind = c_y * Y_ind + c_i1 * I1_ind + c_yi1 * Y_ind * I1_ind
    T01_independent = sp.expand(2 * sp.diff(L_ind, r01) - gcov01 * L_ind)
    T01_independent = sp.expand(T01_independent.subs({r01: q01, r10: q01}))

    residual_unique = sp.expand(T01_unique - T01_independent)
    residual_wrong = sp.expand(T01_wrong - T01_independent)
    return residual_unique, residual_wrong


# ============================================================================
# ფონური მონაცემები
# ============================================================================


def normalize_background_name(name: str | BackgroundName) -> BackgroundName:
    """Accept enum values and case-insensitive strings at the API boundary."""
    if isinstance(name, BackgroundName):
        return name

    try:
        return BackgroundName(str(name).lower())
    except ValueError as exc:
        known = ", ".join(bg.value for bg in BackgroundName)
        raise ValueError(f"unknown background: {name}; expected one of: {known}") from exc


def background(name: str | BackgroundName) -> BackgroundData:
    """coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields."""
    name = normalize_background_name(name)
    if name == BackgroundName.MINKOWSKI:
        t, x, y, z = sp.symbols("t x y z", real=True)
        coords = [t, x, y, z]
        g_cov_diag = [1, -1, -1, -1]
        g_inv_diag = [1, -1, -1, -1]
        sqrt_minus_g = sp.Integer(1)
        fields = [t, x, y, z]
        return BackgroundData(coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields)

    if name == BackgroundName.FLRW:
        t, x, y, z = sp.symbols("t x y z", real=True)
        coords = [t, x, y, z]
        a = sp.Function("a")(t)
        Phi = sp.Function("Phi")(t)
        g_cov_diag = [1, -a**2, -a**2, -a**2]
        g_inv_diag = [1, -1 / a**2, -1 / a**2, -1 / a**2]
        sqrt_minus_g = a**3
        fields = [Phi, x, y, z]
        return BackgroundData(coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields)

    if name == BackgroundName.BIANCHI_I:
        t, x, y, z = sp.symbols("t x y z", real=True)
        coords = [t, x, y, z]
        a = sp.Function("a")(t)
        b = sp.Function("b")(t)
        c = sp.Function("c")(t)
        Phi = sp.Function("Phi")(t)
        g_cov_diag = [1, -a**2, -b**2, -c**2]
        g_inv_diag = [1, -1 / a**2, -1 / b**2, -1 / c**2]
        sqrt_minus_g = a * b * c
        fields = [Phi, x, y, z]
        return BackgroundData(coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields)

    if name == BackgroundName.SCHWARZSCHILD:
        t, r, th, ph = sp.symbols("t r theta phi", real=True, positive=True)
        coords = [t, r, th, ph]
        r_s = sp.Symbol("r_s", real=True, positive=True)
        f = 1 - r_s / r
        g_cov_diag = [f, -1 / f, -r**2, -r**2 * sp.sin(th)**2]
        g_inv_diag = [1 / f, -f, -1 / r**2, -1 / (r**2 * sp.sin(th)**2)]
        sqrt_minus_g = r**2 * sp.sin(th)
        fields = [t, r, th, ph]
        return BackgroundData(coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields)

    raise ValueError(f"unknown background: {name}")


# ============================================================================
# T-ის ცდის ფუნქცია (ფონური ჩასმის წინ მთლიანი L ფონის ფუნქციად)
# ============================================================================


def evaluate_on_background(name, lagrangian_mode="full", t_factor=1):
    """
    ფონური Bianchi/Noether იდენტობის ცდა.

    t_factor: T-ის ცდის გადანამრავლება. სწორი მნიშვნელობა 1; sanity
    check-ისთვის t_factor != 1 — residual მაშინ ცხადად არანულოვანი.

    ცდა:
    1. T_cov[μ,μ] = 2 * dL/dq_μ - L/q_μ  (q_μ = g^{μμ}, q-ის ფუნქცია)
    2. ფონური ჩასმის შემდეგ T_cov-ი ფონის ფუნქციად
    3. T_contra = g_inv ფონური * g_inv ფონური * T_cov (diagonal)
    4. Christoffel ფონური მეტრიკიდან
    5. div = nabla_mu T^{mu nu}
    6. eom_A = (1/sqrt(-g)) d_mu (sqrt(-g) dL/dD_{A,μ})  ფონური ჩასმის შემდეგ
    7. source[nu] = sum_A eom_A * g^{νν} * partial_nu psi_A
    8. residual = div - source — სწორი T-სთვის 0, არასწორი T-სთვის არ-0.
    """
    coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields = background(name)
    q = [sp.Symbol(f"q_{mu}", real=True, nonzero=True) for mu in range(DIM)]

    # ფონური ფუნქცია ცხადად — დარჩება ფონის ფუნქცია, არ-ფიქსირდება როგორც სიმბოლო
    D = [[0 for _ in range(DIM)] for _ in range(NSOLID + 1)]
    D_symbol = [[None for _ in range(DIM)] for _ in range(NSOLID + 1)]
    grad_subs = {}

    for A, field in enumerate(fields):
        for mu, coord in enumerate(coords):
            grad = sp.diff(field, coord)
            if grad != 0:
                sym = sp.Symbol(f"D_{name}_{A}_{mu}", real=True)
                D[A][mu] = sym
                D_symbol[A][mu] = sym
                grad_subs[sym] = grad

    g_inv_diag_symbolic = sp.diag(*q)
    Y, I1, I2, I3, _ = build_invariants_from_metric(g_inv_diag_symbolic, D)
    L, _, _ = get_test_lagrangian(Y, I1, I2, I3, lagrangian_mode)

    # T-ის ცდა q-ის ფუნქცია — q-ის მიმართ ვარიაცია სიმბოლურია; ჩასმა მერე
    T_cov_general = sp.zeros(DIM, DIM)
    for mu in range(DIM):
        T_cov_general[mu, mu] = t_factor * 2 * sp.diff(L, q[mu]) - L / q[mu]

    # current_A_mu სიმბოლური — სიმბოლური ვარიაცია D[A][mu]-ის მიმართ ჯერ
    current_symbolic = [[None for _ in range(DIM)] for _ in range(NSOLID + 1)]
    for A in range(NSOLID + 1):
        for mu in range(DIM):
            if D[A][mu] == 0:
                current_symbolic[A][mu] = sp.Integer(0)
            else:
                current_symbolic[A][mu] = sp.diff(L, D[A][mu])

    # ფონური ჩასმა — q -> g_inv ფონური, D -> ფონური წარმოებული
    subs = {q[mu]: g_inv_diag[mu] for mu in range(DIM)}
    subs.update(grad_subs)

    invariants = tuple(sp.simplify(expr.subs(subs)) for expr in (Y, I1, I2, I3))

    T_cov = sp.zeros(DIM, DIM)
    for mu in range(DIM):
        T_cov[mu, mu] = sp.simplify(T_cov_general[mu, mu].subs(subs))

    T_contra = sp.zeros(DIM, DIM)
    for mu in range(DIM):
        for nu in range(DIM):
            T_contra[mu, nu] = sp.simplify(
                g_inv_diag[mu] * g_inv_diag[nu] * T_cov[mu, nu]
            )

    g_cov = sp.diag(*g_cov_diag)
    g_inv = sp.diag(*g_inv_diag)
    Gamma = get_christoffel(g_cov, g_inv, coords)

    div = [
        covariant_divergence_contra(T_contra, Gamma, coords, nu)
        for nu in range(DIM)
    ]

    eom = []
    for A in range(NSOLID + 1):
        total = 0
        for mu, coord in enumerate(coords):
            if current_symbolic[A][mu] == 0:
                continue
            current_functional = current_symbolic[A][mu].subs(subs)
            integrand = sqrt_minus_g * current_functional
            total += sp.diff(integrand, coord)
        if sqrt_minus_g != 0:
            eom.append(sp.simplify(total / sqrt_minus_g))
        else:
            eom.append(sp.simplify(total))

    source = []
    for nu, coord in enumerate(coords):
        rhs = 0
        for A, field in enumerate(fields):
            rhs += eom[A] * g_inv_diag[nu] * sp.diff(field, coord)
        source.append(sp.simplify(rhs))

    residual = [sp.simplify(div[nu] - source[nu]) for nu in range(DIM)]

    return {
        "coords": coords,
        "fields": fields,
        "g_cov_diag": g_cov_diag,
        "g_inv_diag": g_inv_diag,
        "sqrt_minus_g": sqrt_minus_g,
        "invariants": invariants,
        "T_cov": T_cov,
        "T_contra": T_contra,
        "divergence": div,
        "eom": eom,
        "source": source,
        "residual": residual,
    }


# ============================================================================
# Stress-tensor coefficient sanity check — Noether identity catches a wrong T coefficient
# ============================================================================


def stress_tensor_coefficient_sanity_check():
    """
    სცადე T-ის ცდა გადანამრავლებული t_factor = 3 ფაქტორით.
    residual მაშინ უნდა იყოს არანულოვანი. ეს ამოწმებს სტრეს-ტენზორის
    კოეფიციენტის კონვენციას კონკრეტულ ფონზე; ეს არ არის მთელი თეორიის
    ფალსიფიკაციის ტესტი.
    """
    correct = evaluate_on_background("flrw", lagrangian_mode="reduced_y_i1", t_factor=1)
    wrong = evaluate_on_background("flrw", lagrangian_mode="reduced_y_i1", t_factor=3)
    correct_ok = all(reduce_zero(r) == 0 for r in correct["residual"])
    wrong_fails = any(reduce_zero(r) != 0 for r in wrong["residual"])
    return correct_ok, wrong_fails, correct["residual"], wrong["residual"]


# ============================================================================
# დამხმარე ფუნქციები
# ============================================================================


def convention_summary():
    return {
        "signature": "(+---)",
        "B_AB": "B^{AB} = -g^{mu nu} d_mu phi^A d_nu phi^B",
        "stress_tensor": "T_{mu nu} = 2*dL/dg^{mu nu} - g_{mu nu}*L",
        "phase_relation": (
            "p02_cosmo.py და p01_core.py stress sections იყენებენ ამავე "
            "კონვენციას; mixed convention ცვლილება უნდა დაიბლოკოს NOTATION header-ით"
        ),
        "no_ghost": (
            "იხ. NOTATION.md Full No-Ghost Window. Y-სქემაში K_Y>0 იწერება "
            "Y-ის კოეფიციენტებით; X-სქემაში Y=-2X, ამიტომ იგივე ფიზიკური "
            "კინეტიკური ნიშანი X-ის კოეფიციენტში საპირისპიროდ ჩანს."
        ),
    }


def reduce_zero(expr):
    """ნულის შემოწმება ფონური სიმეტრიების გათვალისწინებით."""
    return sp.factor(sp.together(sp.trigsimp(sp.simplify(expr))))


def is_zero_vector(values):
    return all(reduce_zero(value) == 0 for value in values)


# ============================================================================
# Main
# ============================================================================


if __name__ == "__main__" and _should_run_main_section("stress"):
    print("=" * 72)
    print("PHASE 22 (v3.1): სრული ენერგია-იმპულსის ტენზორი")
    print("=" * 72)

    print("\n1. კონვენციები")
    for key, value in convention_summary().items():
        print(f"  {key:14s}: {value}")

    print("\n2. Off-diagonal ვარიაციის smoke-test")
    offdiag_ok, offdiag_wrong = offdiag_variation_smoke_test()
    print(f"  symmetric ვარიაცია residual: {reduce_zero(offdiag_ok)}")
    print(f"  ფაქტორ-2 ცდომილების residual ნულია? {reduce_zero(offdiag_wrong) == 0}")
    print(f"  ცდის სტატუსი: {'PASS' if reduce_zero(offdiag_ok) == 0 else 'CHECK'}")

    print("\n3. Bianchi/Noether იდენტობა ფონებზე")
    backgrounds = [
        ("minkowski", "full"),
        ("flrw", "full"),
        ("bianchi_i", "full"),
        ("schwarzschild", "reduced_y_i1"),
    ]
    for name, mode in backgrounds:
        result = evaluate_on_background(name, lagrangian_mode=mode)
        residual = [reduce_zero(value) for value in result["residual"]]
        residual_ok = is_zero_vector(residual)
        print(f"\n--- {name} ({mode}) ---")
        print(f"  invariants (Y, I1, I2, I3): {result['invariants']}")
        print(f"  residual vector: {residual}")
        print(f"  status: {'PASS' if residual_ok else 'CHECK'}")

    print("\n4. Stress-tensor coefficient sanity check — Noether იდენტობა ფარდდება T-ის კოეფიციენტს?")
    correct_ok, wrong_fails, correct_res, wrong_res = stress_tensor_coefficient_sanity_check()
    print(f"  სწორი T (t_factor=1): residual ნული? {correct_ok}")
    print(f"  მცდარი T (t_factor=3): residual არანული? {wrong_fails}")
    if correct_ok and wrong_fails:
        print("  ცდის სტატუსი: PASS — იდენტობა კონკრეტულ ფონზე კოეფიციენტის შეცდომას იჭერს")
    else:
        print("  ცდის სტატუსი: CHECK — ცდა არ ფარდდება სწორ/მცდარ T-ს")

    print("\n5. FLRW ენერგია და წნევა")
    flrw = evaluate_on_background("flrw")
    a = sp.Function("a")(sp.Symbol("t", real=True))
    print(f"  rho = {sp.expand(flrw['T_cov'][0, 0])}")
    print(f"  p = T_11/a^2 = {sp.expand(flrw['T_cov'][1, 1] / a**2)}")

    print("\n6. სტატუსი")
    print("  - Minkowski/FLRW/Bianchi I diagonal ფონებზე იდენტობა შემოწმდა სრული L-ით")
    print("  - Schwarzschild diagonal ფონზე — reduced Y+I1 smoke-test")
    print("  - off-diagonal ვარიაცია smoke-test-ით შემოწმდა")
    print("  - sanity check-მა აჩვენა, რომ იდენტობა ამ ფონზე მცდარ T-კოეფიციენტს იჭერს")
    print("  - generic non-diagonal ფონებზე სრული proof ჯერ ღია")
    print("  - შემდეგი ნაბიჯი: generic non-diagonal ფონების proof და სრული perturbation audit")


# ===================== merged from p01_core.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 23: RFG-ის ჩასმა Horndeski/DHOST/ESS ფარგლებში
================================================================================

რეფერენცია: NOTATION.md, phase22, p08_cmb.py

მიზანი (STRATEGY.md ეტაპი II §1):
- RFG-ის ცხადი მაპირება Horndeski-ის G_2, G_3, G_4, G_5 ფუნქციებზე
- DHOST გავრცობის გადახედვა — საიდან მოვა G_3, G_4, G_5
- ESS (Effective Solid State, Endlich-Nicolis-Wang) ფარგლი I_k სოლიდისთვის
- Bellini-Sawicki α_K, α_B, α_M, α_T სრული გადათვლა

დასკვნა:
- RFG დღეს მხოლოდ k-essence (G_2(X) only) Horndeski-ის ქვეჯგუფშია
- I_k სოლიდი არ ჯდება სუფთა Horndeski-ში — საჭიროა ESS გავრცობა
- DHOST-ის ფუნქციები G_4(X), G_5(X) ცარიელ ცდად რჩება

ცენტრალური მაპირება (Bellini-Sawicki კონვენცია X = -½ g^μν ∂_μΦ ∂_νΦ):
    Y = -2X    (NOTATION.md § Horndeski/EFT Map)

RFG ლაგრანჟიანი (NOTATION-ის Y კონვენციით):
    L = c_Y·Y + c_Y2·Y^2 + c_I1·I_1 + c_I1sq·I_1^2 + c_I2·I_2 + c_I3·I_3 + c_YI1·Y·I_1

Horndeski მაპირება:
    G_2(X) = c_Y·(-2X) + c_Y2·(-2X)^2 = -2·c_Y·X + 4·c_Y2·X^2
    G_3 = 0       (no kinetic mixing yet)
    G_4 = M_Pl^2/2  (no Planck mass running yet)
    G_5 = 0       (no Gauss-Bonnet/derivative coupling yet)

I_k სოლიდი ESS framework-ში არის ცალკე ფაქტორი:
    L_solid = c_I1·I_1 + c_I1sq·I_1^2 + c_I2·I_2 + c_I3·I_3 + c_YI1·Y·I_1
    ეს არ ჯდება Horndeski-ის G_i-ში — გადადის ESS-ის ცალკე სტრესის სექტორად.
"""

import sympy as sp


# ============================================================================
# Horndeski მაპირება
# ============================================================================


def rfg_to_horndeski():
    """
    RFG-ის Y-სექტორის ცხადი ჩასმა Horndeski G_2(X)-ში.
    """
    X = sp.Symbol("X", real=True)
    c_Y, c_Y2 = sp.symbols("c_Y c_Y2", real=True)
    M_Pl = sp.Symbol("M_Pl", positive=True)

    # Y = -2X კონვერსია (NOTATION § Horndeski/EFT Map)
    Y_in_X = -2 * X

    # G_2(X) = c_Y·Y + c_Y2·Y^2 -> X-ის ფუნქციაა
    G_2 = c_Y * Y_in_X + c_Y2 * Y_in_X**2
    G_3 = sp.Integer(0)
    G_4 = M_Pl**2 / 2
    G_5 = sp.Integer(0)

    return {
        "G_2": sp.expand(G_2),
        "G_3": G_3,
        "G_4": G_4,
        "G_5": G_5,
        "X_def": "X = -1/2 * g^μν * ∂_μΦ * ∂_νΦ",
        "Y_to_X": "Y = -2X",
    }


# ============================================================================
# Bellini-Sawicki α პარამეტრები
# ============================================================================


def bellini_sawicki_alphas():
    """
    α_K, α_B, α_M, α_T ცხადი ფორმულები phase23-ის Horndeski მაპირებიდან.

    α_T = 2X(G_{4,X} - G_{5,φ}) / M_*^2 + ... = 0  (G_4 const, G_5 = 0)
    α_M = (1/H) * d/dt (ln M_*^2) = 0  (M_*^2 = 2G_4 = M_Pl^2 const)
    α_B = 2*(X*G_{3,X}*φ_dot/H ...) / M_*^2 = 0  (G_3 = 0, G_4 const)
    α_K = (2X*G_{2,X} + 4X^2*G_{2,XX} + ...) / (H^2 * M_*^2)
    """
    X = sp.Symbol("X", real=True)
    H = sp.Symbol("H", real=True, positive=True)
    c_Y, c_Y2, M_Pl = sp.symbols("c_Y c_Y2 M_Pl", real=True, positive=True)

    Y_in_X = -2 * X
    G_2 = c_Y * Y_in_X + c_Y2 * Y_in_X**2
    G_2_X = sp.diff(G_2, X)
    G_2_XX = sp.diff(G_2, X, 2)

    M_star_sq = M_Pl**2

    alpha_T = sp.Integer(0)  # G_4_X = 0, G_5 = 0
    alpha_M = sp.Integer(0)  # M_*^2 const
    alpha_B = sp.Integer(0)  # G_3 = 0, G_4 const
    alpha_K = sp.simplify((2 * X * G_2_X + 4 * X**2 * G_2_XX) / (H**2 * M_star_sq))

    return {
        "alpha_T": alpha_T,
        "alpha_M": alpha_M,
        "alpha_B": alpha_B,
        "alpha_K": alpha_K,
        "G_2_X": G_2_X,
        "G_2_XX": G_2_XX,
    }


# ============================================================================
# I_k სოლიდი — ESS framework
# ============================================================================


def ess_solid_sector():
    """
    I_k სოლიდი არ ჯდება სუფთა Horndeski-ში.
    ESS framework (Endlich-Nicolis-Wang 2012, Ballesteros-Bellazzini 2013).

    L_solid = c_I1*I_1 + c_I1sq*I_1^2 + c_I2*I_2 + c_I3*I_3 + c_YI1*Y*I_1
    """
    return {
        "framework": "ESS (Effective Solid State, Endlich-Nicolis-Wang 2012)",
        "structure": "L_solid = c_I1*I_1 + c_I1sq*I_1^2 + c_I2*I_2 + c_I3*I_3 + c_YI1*Y*I_1",
        "horndeski_compatibility": (
            "I_k-ი არ ჯდება ცხადად G_2-G_5-ში. ESS გავრცობა საჭიროა."
        ),
        "mode_count": (
            "scalar (Φ) + 2 transverse vector phonon (I_k) + 1 longitudinal phonon. "
            "ჯამში 4 propagating mode (Horndeski-ის 1-ის ნაცვლად)."
        ),
        "extra_alpha": (
            "minimal ESS solid changes stress/perturbation equations, not M_*^2; "
            "alpha_M remains zero unless a nonminimal curvature coupling is added."
        ),
    }


# ============================================================================
# DHOST გავრცობა (Beyond Horndeski) — open completion targets
# ============================================================================


def dhost_extension():
    """
    DHOST (Degenerate Higher Order Scalar Tensor) framework.
    Crisostomi-Koyama-Tasinato 2016, Langlois-Noui 2016.

    RFG-ის ფესვი დღეს არ მოიცავს DHOST-ის Class I, II, III ფუნქციებს.
    G_4(X), G_5(X) X-დამოკიდებული ვერსიები ჯერ არ არის დაფარული; ქვემოთ
    ჩამოთვლილია კონკრეტული completion targets.
    """
    return [
        "G_4(X) — X-დამოკიდებული coupling, Brans-Dicke-ის ბუნებრივი გავრცობა",
        "G_5(X) — Gauss-Bonnet-ტიპის, BH regularization-ისთვის ბუნებრივი",
        "DHOST Class I (A_1 - A_5 ფუნქციები) — degenerate higher derivative",
        "Beyond Horndeski (BH) — extra mode-ის გარეშე higher derivatives",
    ]


# ============================================================================
# Main
# ============================================================================


if __name__ == "__main__" and _should_run_main_section("horndeski"):
    print("=" * 72)
    print("PHASE 23: RFG-ის ჩასმა Horndeski/DHOST/ESS ფარგლებში")
    print("რეფერენცია: NOTATION.md, phase22, p08_cmb.py")
    print("=" * 72)

    print("\n1. Horndeski მაპირება (Y = -2X კონვერსია)")
    horndeski = rfg_to_horndeski()
    print(f"  X დეფინიცია: {horndeski['X_def']}")
    print(f"  Y → X: {horndeski['Y_to_X']}")
    print(f"  G_2(X) = {horndeski['G_2']}")
    print(f"  G_3 = {horndeski['G_3']}")
    print(f"  G_4 = {horndeski['G_4']}")
    print(f"  G_5 = {horndeski['G_5']}")

    print("\n2. Bellini-Sawicki α პარამეტრები")
    alphas = bellini_sawicki_alphas()
    print(f"  α_T = {alphas['alpha_T']}    (G_4_X = 0, G_5 = 0)")
    print(f"  α_M = {alphas['alpha_M']}    (M_*^2 = M_Pl^2 const)")
    print(f"  α_B = {alphas['alpha_B']}    (G_3 = 0)")
    print(f"  α_K = {alphas['alpha_K']}")
    print(f"  G_2_X = {alphas['G_2_X']}")
    print(f"  G_2_XX = {alphas['G_2_XX']}")

    print("\n3. I_k სოლიდი — ESS framework (Horndeski-ის გავრცობა)")
    ess = ess_solid_sector()
    for key, value in ess.items():
        print(f"  {key:25s}: {value}")

    print("\n4. DHOST გავრცობა — ღია ცდები")
    for i, task in enumerate(dhost_extension(), 1):
        print(f"  {i}. {task}")

    print("\n5. სტატუსი")
    print("  - RFG = k-essence (G_2(X) only) ქვეჯგუფი Horndeski-ში")
    print("  - α_T, α_M, α_B = 0 ფიქსირდება ცხადად")
    print("  - α_K = (-4*c_Y*X + 48*c_Y2*X²) / (H²·M_Pl²) — მიიღება sympy-დან")
    print("  - I_k სოლიდი ESS-ის გავრცობას მოითხოვს (Horndeski არ ფარდდება)")
    print("  - G_3, G_4(X), G_5 — DHOST გავრცობის შემდეგი ნაბიჯი")
    print("  - p08_cmb.py-ის Bellini-Sawicki ცდა ამ მაპირებას ეთანხმება")


# ===================== merged from p01_core.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 24: Principal symbol, hyperbolicity, Cauchy well-posedness
================================================================================

სტატუსი:
Strategy 3 / M2-ის შესრულება. ეს ფაილი აღარ არის სქელეტი: იგი აგებს
quadratic principal Lagrangian-ს perturbation-ებისთვის და ამოწმებს
characteristic determinant-ს.

მთავარი ობიექტი:
    L2_principal = A*dPhi_dot^2 + B*pi_L_dot^2
                   + C*dPhi_z^2 + D*pi_L_z^2
                   + M*(dPhi_dot*pi_L_z + pi_L_dot*dPhi_z)

Plane wave-ით exp(i(kz - omega t)):
    det M(s) = (A*s + C)*(B*s + D) - M^2*s = 0,
    s = omega^2/k^2.

თუ mixing M=0:
    s_phi = -C/A,   s_L = -D/B.

FLRW-ზე s არის comoving speed squared; physical speed is a^2*s.
Schwarzschild-ზე principal symbol locally orthonormal frame-ში იგივეა, ხოლო
coordinate radial relation არის omega_coord^2 = f(r)^2*c_local^2*k_r^2.
"""

import sympy as sp

# merged import removed: from p01_core import analyze_sound_speeds


def quadratic_principal_coefficients(scale_factor=None):
    """
    Build principal coefficients around FLRW/Minkowski background.

    scale_factor=None means symbolic a. scale_factor=1 gives Minkowski.
    Perturbations propagate along comoving z. The coefficients are the closed
    FLRW scaling of the second-order expansion; at a=1 they are verified
    against p01_core.analyze_sound_speeds().
    """
    a = sp.Symbol("a", positive=True) if scale_factor is None else sp.sympify(scale_factor)

    dPhi_dot, dPhi_z = sp.symbols("dPhi_dot dPhi_z", real=True)
    pi1_dot, pi3_dot = sp.symbols("pi1_dot pi3_dot", real=True)
    pi1_z, pi3_z = sp.symbols("pi1_z pi3_z", real=True)
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )

    A = sp.simplify(c_Y + 6 * c_Y2 + 3 * c_YI1 / a**2)
    B_long = sp.simplify(-c_I1 - c_YI1 - (6 * c_I1sq + 2 * c_I2) / a**2 - c_I3 / a**4)
    C = sp.simplify(-(c_Y + 2 * c_Y2) / a**2 - 3 * c_YI1 / a**4)
    D = sp.simplify((c_I1 + c_YI1) / a**2 + (10 * c_I1sq + 2 * c_I2) / a**4 + c_I3 / a**6)
    M_mix = sp.simplify(4 * c_YI1 / a**2)

    K_T = B_long
    C_T = sp.simplify((c_I1 + c_YI1) / a**2 + (6 * c_I1sq + c_I2) / a**4)
    cs2_T_comoving = sp.simplify(-C_T / K_T)
    L2 = sp.simplify(
        A * dPhi_dot**2
        + B_long * pi3_dot**2
        + C * dPhi_z**2
        + D * pi3_z**2
        + M_mix * (dPhi_dot * pi3_z + pi3_dot * dPhi_z)
        + K_T * pi1_dot**2
        + C_T * pi1_z**2
    )

    return {
        "a": a,
        "L2": L2,
        "A": A,
        "B_long": B_long,
        "C": C,
        "D": D,
        "M_mix": M_mix,
        "K_T": K_T,
        "C_T": C_T,
        "cs2_T_comoving": cs2_T_comoving,
        "cs2_T_physical": sp.simplify(a**2 * cs2_T_comoving),
    }


def characteristic_polynomial(coeffs, solve_symbolic=False):
    """det M(s)=0 for mixed phase-longitudinal sector."""
    s = sp.Symbol("s", real=True)
    A = coeffs["A"]
    B = coeffs["B_long"]
    C = coeffs["C"]
    D = coeffs["D"]
    M = coeffs["M_mix"]
    det = sp.factor(sp.simplify((A * s + C) * (B * s + D) - M**2 * s))
    roots = sp.solve(det, s) if solve_symbolic else None
    return s, det, roots


def mixed_mode_stability_conditions(coeffs):
    """
    Necessary and sufficient algebraic conditions for the 2x2 principal
    mixed phase-longitudinal sector on the chosen homogeneous background.

    For
        det M(s) = p2*s^2 + p1*s + p0,
    with s=omega^2/k^2, the two mixed speeds are real and positive iff:
        p2 > 0, p1 < 0, p0 > 0, discriminant >= 0.

    This closes the local 2x2 principal-symbol test.  It does not close the
    full curved-background/global hyperbolicity problem.
    """
    A = coeffs["A"]
    B = coeffs["B_long"]
    C = coeffs["C"]
    D = coeffs["D"]
    M = coeffs["M_mix"]
    K_T = coeffs["K_T"]
    C_T = coeffs["C_T"]

    p2 = sp.factor(sp.simplify(A * B))
    p1 = sp.factor(sp.simplify(A * D + B * C - M**2))
    p0 = sp.factor(sp.simplify(C * D))
    discriminant = sp.factor(sp.simplify(p1**2 - 4 * p2 * p0))

    return {
        "mixed_polynomial_coefficients": {
            "p2": p2,
            "p1": p1,
            "p0": p0,
            "discriminant": discriminant,
        },
        "no_ghost_required": {
            "A": A,
            "B_long": B,
            "K_T": K_T,
            "conditions": "A>0, B_long>0, K_T>0",
        },
        "mixed_speed_required": {
            "conditions": "p2>0, p1<0, p0>0, discriminant>=0",
            "meaning": "two real positive mixed eigen-speeds s=omega^2/k^2",
        },
        "transverse_required": {
            "C_T": C_T,
            "condition": "with K_T>0 require C_T<0 so c_T^2=-C_T/K_T>0",
        },
        "scope": "local homogeneous principal-symbol closure; not a full curved/global proof",
    }


def mixed_mode_numeric_condition_check():
    """Compare the algebraic mixed-mode criteria against the existing numeric cases."""
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )
    coeffs, _s, _det, _roots = minkowski_principal_symbol()
    conditions = mixed_mode_stability_conditions(coeffs)
    poly = conditions["mixed_polynomial_coefficients"]

    rows = []
    for row in numeric_hyperbolicity_cases():
        name = row["name"]
        if name == "stable_decoupled":
            subs = {c_Y: 1.0, c_Y2: 0.10, c_I1: -2.0, c_I1sq: 0.10, c_I2: -0.10, c_I3: 0.0, c_YI1: 0.0}
        elif name == "stable_mixed_small":
            subs = {c_Y: 1.0, c_Y2: 0.10, c_I1: -2.0, c_I1sq: 0.10, c_I2: -0.10, c_I3: 0.0, c_YI1: 0.05}
        elif name == "ghost_fail_phase":
            subs = {c_Y: -1.0, c_Y2: 0.05, c_I1: -2.0, c_I1sq: 0.10, c_I2: -0.10, c_I3: 0.0, c_YI1: 0.0}
        elif name == "gradient_fail_solid":
            subs = {c_Y: 1.0, c_Y2: 0.10, c_I1: -0.20, c_I1sq: 0.00, c_I2: 0.50, c_I3: 0.0, c_YI1: 0.0}
        else:
            continue

        values = {
            "p2": float(sp.N(poly["p2"].subs(subs))),
            "p1": float(sp.N(poly["p1"].subs(subs))),
            "p0": float(sp.N(poly["p0"].subs(subs))),
            "discriminant": float(sp.N(poly["discriminant"].subs(subs))),
        }
        algebraic_pass = (
            values["p2"] > 0
            and values["p1"] < 0
            and values["p0"] > 0
            and values["discriminant"] >= 0
        )
        rows.append({
            "name": name,
            **values,
            "mixed_algebraic_status": "PASS" if algebraic_pass else "FAIL",
            "root_status": row["mixed_roots_status"],
        })
    return rows


def minkowski_principal_symbol():
    coeffs = quadratic_principal_coefficients(scale_factor=1)
    s, det, roots = characteristic_polynomial(coeffs)
    return coeffs, s, det, roots


def flrw_principal_symbol():
    coeffs = quadratic_principal_coefficients(scale_factor=None)
    s, det, roots = characteristic_polynomial(coeffs)
    physical_det = sp.factor(sp.simplify(det.subs(s, sp.Symbol("c_phys2", real=True) / coeffs["a"]**2)))
    return coeffs, s, det, roots, physical_det


def decoupled_mixed_roots(coeffs):
    """Closed roots when c_YI1=0 removes phase-longitudinal mixing."""
    c_YI1 = sp.Symbol("c_YI1", real=True)
    A = coeffs["A"].subs(c_YI1, 0)
    B = coeffs["B_long"].subs(c_YI1, 0)
    C = coeffs["C"].subs(c_YI1, 0)
    D = coeffs["D"].subs(c_YI1, 0)
    return sp.simplify(-C / A), sp.simplify(-D / B)


def schwarzschild_local_symbol():
    """
    Local orthonormal principal symbol outside horizon.

    The local Cauchy problem is governed by the same determinant as Minkowski.
    Coordinate radial propagation is redshifted:
        omega_coord^2 = f(r)^2 * c_local^2 * k_r^2.
    """
    r, r_s, c_local2, k_r = sp.symbols("r r_s c_local2 k_r", positive=True)
    omega = sp.Symbol("omega", real=True)
    f = 1 - r_s / r
    radial_dispersion = sp.Eq(omega**2, sp.simplify(f**2 * c_local2 * k_r**2))
    return {
        "f": f,
        "local_det": minkowski_principal_symbol()[2],
        "radial_coordinate_dispersion": radial_dispersion,
        "horizon_note": "coordinate speed -> 0 as f -> 0, local characteristic remains finite",
    }


def verify_against_phase1():
    """
    Compare phase24 coefficients with p01_core.analyze_sound_speeds().
    We compare coefficients, not phase1's printed characteristic convention.
    """
    phase24 = minkowski_principal_symbol()[0]
    cs2_T_phase1, _eq_phase1, coeffs_phase1, _roots_phase1 = analyze_sound_speeds()
    checks = {
        "A_difference": sp.simplify(phase24["A"] - coeffs_phase1["A"]),
        "B_long_difference": sp.simplify(phase24["B_long"] - coeffs_phase1["B_pi3"]),
        "C_difference": sp.simplify(phase24["C"] - coeffs_phase1["C"]),
        "D_difference": sp.simplify(phase24["D"] - coeffs_phase1["D"]),
        "M_mix_difference": sp.simplify(phase24["M_mix"] - coeffs_phase1["M_mix"]),
        "transverse_speed_difference": sp.simplify(phase24["cs2_T_comoving"] - cs2_T_phase1),
    }
    return checks


def numeric_hyperbolicity_cases():
    """
    Numeric smoke-test. PASS means:
        kinetic coefficients are positive,
        mixed roots are real and positive,
        transverse speed is positive.
    """
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )
    coeffs, s, det, _ = minkowski_principal_symbol()

    cases = [
        {
            "name": "stable_decoupled",
            c_Y: 1.0,
            c_Y2: 0.10,
            c_I1: -2.0,
            c_I1sq: 0.10,
            c_I2: -0.10,
            c_I3: 0.0,
            c_YI1: 0.0,
        },
        {
            "name": "stable_mixed_small",
            c_Y: 1.0,
            c_Y2: 0.10,
            c_I1: -2.0,
            c_I1sq: 0.10,
            c_I2: -0.10,
            c_I3: 0.0,
            c_YI1: 0.05,
        },
        {
            "name": "ghost_fail_phase",
            c_Y: -1.0,
            c_Y2: 0.05,
            c_I1: -2.0,
            c_I1sq: 0.10,
            c_I2: -0.10,
            c_I3: 0.0,
            c_YI1: 0.0,
        },
        {
            "name": "gradient_fail_solid",
            c_Y: 1.0,
            c_Y2: 0.10,
            c_I1: -0.20,
            c_I1sq: 0.00,
            c_I2: 0.50,
            c_I3: 0.0,
            c_YI1: 0.0,
        },
    ]

    rows = []
    for case in cases:
        name = case["name"]
        subs = {key: value for key, value in case.items() if key != "name"}
        kinetic_values = [
            float(sp.N(coeffs["A"].subs(subs))),
            float(sp.N(coeffs["B_long"].subs(subs))),
            float(sp.N(coeffs["K_T"].subs(subs))),
        ]
        roots = [complex(root) for root in sp.nroots(sp.Poly(det.subs(subs), s), n=20, maxsteps=100)]
        roots_real = [root.real for root in roots if abs(root.imag) < 1e-7]
        transverse = float(sp.N(coeffs["cs2_T_comoving"].subs(subs)))

        kinetic_pass = all(value > 0 for value in kinetic_values)
        roots_real_pass = len(roots_real) == 2
        roots_positive = roots_real_pass and all(value > 0 for value in roots_real)
        transverse_pass = transverse > 0
        overall = kinetic_pass and roots_positive and transverse_pass

        rows.append(
            {
                "name": name,
                "kinetic_values": kinetic_values,
                "mixed_roots": roots_real,
                "transverse_cs2": transverse,
                "kinetic": "PASS" if kinetic_pass else "FAIL",
                "mixed_roots_status": "PASS" if roots_positive else "FAIL",
                "transverse_status": "PASS" if transverse_pass else "FAIL",
                "overall": "PASS" if overall else "FAIL",
            }
        )
    return rows


def status_assessment():
    return {
        "minkowski": "det M(s)=0 computed; mixed-mode algebraic positivity criteria are explicit",
        "flrw": "comoving det M(s)=0 computed; physical speed is a^2*s; same algebraic criteria apply after scaling",
        "schwarzschild": "local orthonormal determinant equals Minkowski; coordinate radial redshift added",
        "remaining": "global hyperbolicity and full curved-background perturbation system remain open",
    }


def p01_proof_gap_register():
    """Direct list of p01 axiom boundaries and weaknesses that must not be hidden."""
    return [
        {
            "gap": "foundational_measurability_axiom_boundary",
            "current_status": (
                "Y=1, B=delta is explicitly marked as the effective normalized "
                "background representation of the Intuitive_Theory.md §0 "
                "measurability axiom"
            ),
            "risk": (
                "the only risk is mislabeling this axiom as a p01-derived "
                "polynomial minimum; the homogeneous modulus-extremum route "
                "would give K_pi=0 if incorrectly imposed"
            ),
            "next_step": (
                "keep it declared as a foundational axiom and strengthen only "
                "the effective consistency tests around it"
            ),
        },
        {
            "gap": "mixed_mode_stability",
            "current_status": "local 2x2 principal-symbol algebraic criteria are now explicit",
            "risk": "this closes the local homogeneous symbol, not the full curved/global Cauchy problem",
            "next_step": "extend to curved background perturbations and prove global hyperbolicity conditions",
        },
        {
            "gap": "lorentz_invariance",
            "current_status": "one x-boost background-stress condition plus PPN comparison",
            "risk": "all boost directions and perturbation sector are not audited",
            "next_step": "derive the full boost-direction tensor condition and perturbation-sector Lorentz audit",
        },
        {
            "gap": "generic_noether_identity",
            "current_status": "diagonal backgrounds plus off-diagonal coefficient sanity-check",
            "risk": "generic non-diagonal metric proof is not written",
            "next_step": "prove the Noether identity for generic symmetric metric variables",
        },
        {
            "gap": "dirac_bergmann_closure",
            "current_status": "candidate DOF count only",
            "risk": "second-class bracket non-degeneracy/anomaly closure is open",
            "next_step": "construct the constraint matrix and check non-degeneracy/anomaly closure",
        },
    ]


def article_core_theorem():
    """
    Article-facing theorem ledger for the p01 core.

    This function is the clean bridge from p01_core.py into the Georgian
    article draft: it exposes only the action/sign/stability facts that are
    actually supported inside this file.
    """
    Y, I1, I2, I3 = init_variables()
    L_poly = get_polynomial_lagrangian(Y, I1, I2, I3)
    K_Phi_c, K_pi_c = analyze_lorentz_constrained_stability()
    horndeski_map = rfg_to_horndeski()
    alphas = bellini_sawicki_alphas()
    coeffs_m, _s, _det, _roots = minkowski_principal_symbol()
    mixed_conditions = mixed_mode_stability_conditions(coeffs_m)

    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    return {
        "article_use": "core action, sign convention, and necessary local stability gate",
        "postulate_boundary": (
            "Y=1, B^{AB}=delta^{AB} is the normalized effective background, "
            "not a polynomial-minimum derivation."
        ),
        "invariants": {
            "Y": "g^{mu nu} d_mu Phi d_nu Phi",
            "B_AB": "-g^{mu nu} d_mu phi^A d_nu phi^B",
            "I1": "tr(B)",
            "I2": "1/2*(I1^2-tr(B^2))",
            "I3": "det(B)",
        },
        "polynomial_lagrangian": L_poly,
        "sign_bridge": {
            "Y_to_X": horndeski_map["Y_to_X"],
            "c_X": "-2*c_Y^(Y)",
            "c_X2": "4*c_Y2^(Y)",
        },
        "horndeski_Y_sector": {
            "G_2": horndeski_map["G_2"],
            "G_3": horndeski_map["G_3"],
            "G_4": horndeski_map["G_4"],
            "G_5": horndeski_map["G_5"],
            "alpha_T": alphas["alpha_T"],
            "alpha_M": alphas["alpha_M"],
            "alpha_B": alphas["alpha_B"],
            "alpha_K": alphas["alpha_K"],
        },
        "lorentz_branch_relations": {
            "c_I1sq": sp.Eq(sp.Symbol("c_I1sq", real=True), sp.Symbol("c_Y2", real=True)),
            "c_I1": sp.Eq(
                sp.Symbol("c_I1", real=True),
                sp.Symbol("c_Y", real=True)
                - 4 * sp.Symbol("c_Y2", real=True)
                + 2 * sp.Symbol("c_YI1", real=True)
                - 2 * sp.Symbol("c_I2", real=True)
                - sp.Symbol("c_I3", real=True),
            ),
        },
        "necessary_no_ghost_window": {
            "K_PhiPhi_after_relations": K_Phi_c,
            "K_pipi_after_relations": K_pi_c,
            "article_window": [
                sp.Gt(c_Y2, 0),
                sp.Gt(c_Y + 3 * c_YI1, -6 * c_Y2),
                sp.Lt(c_Y + 3 * c_YI1, -2 * c_Y2),
            ],
        },
        "mixed_mode_gate": mixed_conditions,
        "article_status": {
            "action": "CLOSED_MINIMAL_POLYNOMIAL",
            "sign_convention": "CLOSED_Y_TO_X_BRIDGE",
            "no_ghost": "NECESSARY_LOCAL_WINDOW",
            "mixed_modes": "LOCAL_PRINCIPAL_SYMBOL_CRITERIA_AVAILABLE",
            "global_stability": "SEPARATE_PROOF_TARGET",
            "dof_count": "CANDIDATE_LEDGER_ONLY",
        },
    }


def compact_det_label():
    return "(A*s + C)*(B*s + D) - M_mix**2*s = 0"


if __name__ == "__main__" and _should_run_main_section("hyperbolicity"):
    print("=" * 72)
    print("PHASE 24: Principal symbol + hyperbolicity")
    print("=" * 72)

    print("\n1. Minkowski principal symbol")
    coeffs_m, s, det_m, roots_m = minkowski_principal_symbol()
    for key in ["A", "B_long", "C", "D", "M_mix", "K_T", "C_T"]:
        print(f"  {key:10s}: {coeffs_m[key]}")
    print(f"  det M(s): {compact_det_label()}")
    print(f"  expanded degree in s: {sp.degree(det_m, s)}")
    print(f"  decoupled roots (c_YI1=0): {decoupled_mixed_roots(coeffs_m)}")
    print(f"  transverse c_s^2: {coeffs_m['cs2_T_comoving']}")
    mixed_conditions = mixed_mode_stability_conditions(coeffs_m)
    print("\n1b. Mixed-mode algebraic stability criterion")
    for key, value in mixed_conditions["mixed_polynomial_coefficients"].items():
        print(f"  {key:12s}: {value}")
    print("  mixed condition:", mixed_conditions["mixed_speed_required"]["conditions"])
    print("  no-ghost condition:", mixed_conditions["no_ghost_required"]["conditions"])
    print("  transverse condition:", mixed_conditions["transverse_required"]["condition"])
    print("  scope:", mixed_conditions["scope"])

    print("\n2. FLRW principal symbol")
    coeffs_f, s_f, det_f, roots_f, det_phys = flrw_principal_symbol()
    print(f"  a: {coeffs_f['a']}")
    for key in ["A", "B_long", "C", "D", "M_mix"]:
        print(f"  {key:10s}: {coeffs_f[key]}")
    print(f"  det_comoving M(s): {compact_det_label()}")
    print(f"  expanded degree in s: {sp.degree(det_f, s_f)}")
    print("  physical substitution: s = c_phys^2/a^2")
    print(f"  transverse c_phys^2: {coeffs_f['cs2_T_physical']}")

    print("\n3. Schwarzschild local principal symbol")
    sch = schwarzschild_local_symbol()
    print(f"  f                             : {sch['f']}")
    print(f"  local_det                     : {compact_det_label()}")
    print(f"  radial_coordinate_dispersion  : {sch['radial_coordinate_dispersion']}")
    print(f"  horizon_note                  : {sch['horizon_note']}")

    print("\n4. phase1 coefficient verification")
    for key, value in verify_against_phase1().items():
        print(f"  {key:30s}: {value}")

    print("\n5. Numeric hyperbolicity smoke-test")
    for row in numeric_hyperbolicity_cases():
        roots_text = ", ".join(f"{value:.4g}" for value in row["mixed_roots"])
        kinetic_text = ", ".join(f"{value:.4g}" for value in row["kinetic_values"])
        print(f"  {row['name']:22s}: overall={row['overall']}")
        print(f"    K values: {kinetic_text} -> {row['kinetic']}")
        print(f"    mixed roots: {roots_text} -> {row['mixed_roots_status']}")
        print(f"    transverse c_s^2: {row['transverse_cs2']:.4g} -> {row['transverse_status']}")

    print("\n5b. Numeric check against algebraic mixed-mode criterion")
    for row in mixed_mode_numeric_condition_check():
        print(
            f"  {row['name']:22s}: p2={row['p2']:.4g}, p1={row['p1']:.4g}, "
            f"p0={row['p0']:.4g}, disc={row['discriminant']:.4g} -> "
            f"{row['mixed_algebraic_status']} (roots: {row['root_status']})"
        )

    print("\n6. Status")
    for key, value in status_assessment().items():
        print(f"  {key:14s}: {value}")


# ===================== merged from p01_core.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 25: EFT-of-Dark-Energy — necessary no-ghost window + Bellini-Sawicki α-ები
================================================================================

სტატუსი:
Strategy 3 / X4+M3-ის შესრულება.

ამ ფაილის მიზანი:
1. phase21/phase23-ის ghost კონფლიქტის მკაცრი დახურვა:
       α_K = (-4*c_Y*X + 48*c_Y2*X^2)/(H^2*M_Pl^2)
   და α_K > 0 მოთხოვნის ცალკე ჩვენება X-სქემაში.

2. Y-სქემაში ფაზური background-dependent no-ghost აუცილებელი პირობა:
       K_Phi = q00 * (c_Y + 6*c_Y2*Y0 + c_YI1*I1_bg) > 0

   FLRW normalized background:
       K_Phi = c_Y + 6*c_Y2 + 3*c_YI1/a^2 > 0
       a=1 -> c_Y + 6*c_Y2 + 3*c_YI1 > 0

3. Schwarzschild local/static ფონზე იგივე პირობის smoke-test.
   ადგილობრივი ორთონორმალური ჩარჩო უბრუნდება Minkowski/FLRW local პირობას;
   coordinate Phi=t smoke-test ცალკე იბეჭდება, რადგან კოორდინატული ნორმალიზაცია
   ფიზიკური no-ghost პირობა არ არის.

4. Solid sector აღარ იკარგება perturbation ledger-ში. მინიმალური coupling-ის
   შემთხვევაში იგი არ ცვლის M_*^2-ს და alpha_M-ს; alpha_B/alpha_K ტიპის
   ეფექტური წვლილი საჭიროებს ESS perturbation derivation-ს. სრული CMB fit
   მაინც phase21/hi_class ამოცანად რჩება.
"""

import sympy as sp


def horndeski_y_sector_alphas():
    """
    Bellini-Sawicki X scheme, Y = -2X.

    Pure Y-sector:
        G2(X) = -2*c_Y*X + 4*c_Y2*X^2
        alpha_K = (2X*G2_X + 4X^2*G2_XX)/(H^2*M_Pl^2)
                = (-4*c_Y*X + 48*c_Y2*X^2)/(H^2*M_Pl^2)

    This block does not impose a new sign convention.  It translates the
    Y-scheme coefficient into the standard X variable.  Since Y=-2X, the
    coefficient multiplying X carries the opposite sign from the coefficient
    multiplying Y.

    If I1 is treated as a fixed background value, c_YI1*Y*I1 contributes
    by c_Y -> c_Y + c_YI1*I1_bg in G2.
    """
    X, I1_bg = sp.symbols("X I1_bg", real=True)
    H, M_Pl = sp.symbols("H M_Pl", positive=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    y_in_x = -2 * X
    G2_pure = c_Y * y_in_x + c_Y2 * y_in_x**2
    G2_with_i1 = (c_Y + c_YI1 * I1_bg) * y_in_x + c_Y2 * y_in_x**2

    def alpha_k(G2):
        G2_X = sp.diff(G2, X)
        G2_XX = sp.diff(G2, X, 2)
        return sp.simplify((2 * X * G2_X + 4 * X**2 * G2_XX) / (H**2 * M_Pl**2))

    return {
        "G2_pure": sp.expand(G2_pure),
        "alpha_T_Y_sector": sp.Integer(0),
        "alpha_M_Y_sector": sp.Integer(0),
        "alpha_B_Y_sector": sp.Integer(0),
        "alpha_K_pure": alpha_k(G2_pure),
        "alpha_K_with_I1_background": alpha_k(G2_with_i1),
        "Y_to_X_sign_bridge": (
            "Y=-2X. A positive Y-scheme phase coefficient appears with the "
            "opposite sign in the X-scheme G2_X coefficient; compare kinetic "
            "positivity only after choosing one scheme."
        ),
        "ghost_rule_X_scheme": (
            "require alpha_K_total > 0 in the X scheme; this is the translated "
            "form of the Y-scheme no-ghost window, not an independent c_Y sign rule."
        ),
    }


def y_scheme_no_ghost_conditions():
    """
    General Y-scheme quadratic coefficient for phase perturbation.

    Let Phi = Phi_bg + pi and Y = q00*(1 + pi_dot)^2 on a static time-like
    background. The pi_dot^2 coefficient is:

        K_Phi = q00 * (c_Y + 6*c_Y2*Y0 + c_YI1*I1_bg)

    q00 > 0 outside horizons, so the bracket controls the sign.
    """
    a, r, r_s, theta = sp.symbols("a r r_s theta", positive=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)
    f = 1 - r_s / r

    K_general = sp.Symbol("q00", positive=True) * (
        c_Y + 6 * c_Y2 * sp.Symbol("Y0", positive=True) + c_YI1 * sp.Symbol("I1_bg", positive=True)
    )

    K_flrw = sp.simplify(c_Y + 6 * c_Y2 + 3 * c_YI1 / a**2)
    K_flrw_today = sp.simplify(K_flrw.subs(a, 1))
    flrw_cYI1_bound = sp.solve_univariate_inequality(K_flrw_today > 0, c_YI1)
    flrw_cY2_bound = sp.solve_univariate_inequality(K_flrw_today > 0, c_Y2)

    # Schwarzschild, local orthonormal static frame: Y0=1, I1=3.
    K_schw_local = K_flrw_today

    # Coordinate smoke-test for Phi=t and solid labels (r, theta, phi), outside r>r_s.
    I1_schw_coord = sp.simplify(f + 1 / r**2 + 1 / (r**2 * sp.sin(theta) ** 2))
    K_schw_coord_bracket = sp.simplify(c_Y + 6 * c_Y2 / f + c_YI1 * I1_schw_coord)
    K_schw_coord = sp.simplify(K_schw_coord_bracket / f)
    K_schw_coord_equator = sp.simplify(K_schw_coord.subs(theta, sp.pi / 2))

    return {
        "K_general": K_general,
        "K_FLRW": K_flrw,
        "K_FLRW_today": K_flrw_today,
        "FLRW_today_c_YI1_bound": flrw_cYI1_bound,
        "FLRW_today_c_Y2_bound": flrw_cY2_bound,
        "K_Schwarzschild_local": K_schw_local,
        "I1_Schwarzschild_coordinate": I1_schw_coord,
        "K_Schwarzschild_coordinate_equator": K_schw_coord_equator,
    }


def sign_window_sweep():
    """
    Numeric smoke-test for c_Y2 > 0 and c_YI1 window.

    These are not fitted values; they are small examples that show the inequality
    catches pass/fail branches.
    """
    a_value = 1.0
    r_value = 10.0
    r_s_value = 2.0
    theta_value = sp.pi / 2

    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)
    a, r, r_s, theta = sp.symbols("a r r_s theta", positive=True)
    f = 1 - r_s / r

    K_flrw = c_Y + 6 * c_Y2 + 3 * c_YI1 / a**2
    K_schw_local = c_Y + 6 * c_Y2 + 3 * c_YI1
    I1_schw = f + 1 / r**2 + 1 / (r**2 * sp.sin(theta) ** 2)
    K_schw = (c_Y + 6 * c_Y2 / f + c_YI1 * I1_schw) / f

    cases = [
        {"name": "healthy_example", "c_Y": 1.0, "c_Y2": 0.10, "c_YI1": 0.0},
        {"name": "mixed_term_too_negative", "c_Y": 1.0, "c_Y2": 0.10, "c_YI1": -1.0},
        {"name": "c_Y2_negative_risky", "c_Y": 1.0, "c_Y2": -0.30, "c_YI1": 0.0},
    ]

    rows = []
    for case in cases:
        subs = {
            c_Y: case["c_Y"],
            c_Y2: case["c_Y2"],
            c_YI1: case["c_YI1"],
            a: a_value,
            r: r_value,
            r_s: r_s_value,
            theta: theta_value,
        }
        flrw_value = float(K_flrw.subs(subs))
        schw_local_value = float(K_schw_local.subs(subs))
        schw_value = float(K_schw.subs(subs))
        rows.append(
            {
                **case,
                "K_FLRW_today": flrw_value,
                "FLRW_status": "PASS" if flrw_value > 0 else "FAIL",
                "K_Schwarzschild_local": schw_local_value,
                "Schwarzschild_local_status": "PASS" if schw_local_value > 0 else "FAIL",
                "K_Schwarzschild_coord": schw_value,
                "Schwarzschild_coord_status": "PASS" if schw_value > 0 else "FAIL",
            }
        )
    return rows


def ess_solid_alpha_bookkeeping():
    """
    ESS/Ballesteros-Bellazzini style solid-sector bookkeeping.

    Minimal solid matter does not run the Planck mass by itself.  Therefore
    alpha_M and alpha_T stay zero in this minimal branch; possible alpha_B/alpha_K
    stress-sector terms remain explicit open functions until the ESS perturbation
    derivation is written.
    """
    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    H = sp.Function("H")(t)
    M_Pl = sp.Symbol("M_Pl", positive=True)
    X = sp.Symbol("X", real=True)
    I1_bg = sp.Symbol("I1_bg", positive=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    alpha_B_solid = sp.Function("alpha_B_solid")(a)
    alpha_K_solid = sp.Function("alpha_K_solid")(a)
    delta_M2_nonminimal = sp.Function("delta_M2_nonminimal")(a)

    M_eff_sq = M_Pl**2
    alpha_K_y_i1 = sp.simplify(
        (-4 * X * (c_Y + c_YI1 * I1_bg) + 48 * c_Y2 * X**2) / (H**2 * M_eff_sq)
    )
    alpha_M_nonminimal_candidate = sp.simplify(
        sp.diff(sp.log(M_Pl**2 + delta_M2_nonminimal), t) / H
    )

    return {
        "M_eff_sq": M_eff_sq,
        "alpha_K_total": alpha_K_y_i1 + alpha_K_solid,
        "alpha_B_total": alpha_B_solid,
        "alpha_M_total_minimal": sp.Integer(0),
        "alpha_T_total_minimal": sp.Integer(0),
        "alpha_M_nonminimal_candidate": alpha_M_nonminimal_candidate,
        "GW170817_filter": "minimal branch gives alpha_T=0; nonminimal tensor terms must satisfy |alpha_T| < O(1e-15)",
        "solid_note": "alpha_B_solid/alpha_K_solid require ESS perturbation derivation before CLASS/hi_class fit",
    }


def observational_filters():
    return {
        "alpha_T": "|alpha_T| < O(1e-15) from GW170817/GRB170817A; minimal branch has alpha_T=0",
        "alpha_K": "alpha_K_total > 0 is necessary; full kinetic matrix/eigenvalues still required",
        "alpha_M_alpha_B": "minimal alpha_M=0; alpha_B_solid must be derived before CMB/LSS/BAO fit",
        "DESI_link": "static Lambda_eff is not enough for w(z); dynamic alpha-sector is needed",
    }


def class_camb_interface_open():
    return [
        "export alpha_K(a), alpha_B(a), alpha_M(a), alpha_T(a) arrays after ESS closure",
        "choose ESS closure for alpha_B_solid(a), alpha_K_solid(a); keep alpha_M=0 unless nonminimal coupling is added",
        "run hi_class/CLASS Planck 2018 likelihood",
        "add BAO/LSS/DESI likelihoods after background H(a) is fixed",
    ]


if __name__ == "__main__" and _should_run_main_section("eft"):
    print("=" * 72)
    print("PHASE 25: EFT-of-Dark-Energy — no-ghost + alpha sweep")
    print("=" * 72)

    print("\n1. Horndeski/Y-sector Bellini-Sawicki alphas")
    alphas = horndeski_y_sector_alphas()
    for key, value in alphas.items():
        print(f"  {key:30s}: {value}")

    print("\n2. Y-scheme phase no-ghost necessary conditions")
    conditions = y_scheme_no_ghost_conditions()
    for key, value in conditions.items():
        print(f"  {key:34s}: {value}")

    print("\n3. c_Y2 / c_YI1 sign-window smoke-test")
    for row in sign_window_sweep():
        print(
            f"  {row['name']:24s}: "
            f"c_Y={row['c_Y']:+.2f}, c_Y2={row['c_Y2']:+.2f}, c_YI1={row['c_YI1']:+.2f} | "
            f"FLRW K={row['K_FLRW_today']:+.3f} {row['FLRW_status']} | "
            f"Schw local K={row['K_Schwarzschild_local']:+.3f} {row['Schwarzschild_local_status']} | "
            f"coord smoke K={row['K_Schwarzschild_coord']:+.3f} {row['Schwarzschild_coord_status']}"
        )
    print("  note: Schwarzschild coord smoke is not the physical ghost verdict; local K controls the sign.")

    print("\n4. ESS solid-sector alpha bookkeeping")
    ess = ess_solid_alpha_bookkeeping()
    for key, value in ess.items():
        print(f"  {key:22s}: {value}")

    print("\n5. Observational filters")
    for key, value in observational_filters().items():
        print(f"  {key:18s}: {value}")

    print("\n6. CLASS/hi_class open interface")
    for i, task in enumerate(class_camb_interface_open(), 1):
        print(f"  {i}. {task}")

    print("\n7. Status")
    print("  - Strategy 3 X4: background-dependent phase no-ghost window is now explicit.")
    print("  - Strategy 3 M3: alpha_K formula is written explicitly and solid-sector deltas are visible.")
    print("  - Full ESS perturbation derivation and Planck chi^2 fit remain phase21/hi_class work.")


# ===================== merged from p01_core.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 34: ნორმალიზებული სუპერსოლიდური ფონი — ლოურენცის background-stress შემოწმება
================================================================================

სტატუსი:
ეს ბლოკი წარმოადგენს Priority A / X1-ის ფორმალური შემოწმების სამუშაო ბლოკს.
ამოცანა: შევამოწმოთ, რა ალგებრული პირობა სჭირდება ნორმალიზებულ
სუპერსოლიდურ ფონს (სადაც φ^A = x^A), რომ x-boost-ის შემდეგ მისი
background stress დარჩეს იზოტროპული.

თუ background-stress დონეზე ფონი ლოურენც-სიმეტრიულად იკითხება, მაშინ ბუსტირებულ ათვლის სისტემაშიც მისი 
ენერგია-იმპულსის ტენზორი უნდა დარჩეს T_μν ∝ η_μν ფორმის. კერძოდ:
1. არ უნდა გაჩნდეს იმპულსის ნაკადი: T_01 = 0
2. არ უნდა გაჩნდეს სივრცული ანიზოტროპია: T_11 - T_22 = 0

ამ კოდში SymPy-ს მეშვეობით ვასრულებთ ფონური ველების (Φ = t, φ^A = x^A) 
ბუსტს v სიჩქარით, ვითვლით სრულ T_μν-ს და გამოგვაქვს ის ზუსტი ალგებრული 
პირობა კოეფიციენტებზე, რომელიც ანულებს T_01-ს.

შედეგი: x-boost background stress-ის დონეზე ნორმალიზებული სუპერსოლიდური ფონი ლოურენც-სიმეტრიულად იკითხება
მხოლოდ მაშინ, თუ სრულდება კონკრეტული კონსტრეინტი. სრული perturbation-sector
Lorentz audit ცალკე დასახური რჩება.
"""

import sympy as sp

def analyze_lorentz_boost():
    v = sp.Symbol('v', real=True)
    gamma = 1 / sp.sqrt(1 - v**2)
    t, x, y, z = sp.symbols('t x y z', real=True)

    # ველების ბუსტი x ღერძის გასწვრივ
    Phi = gamma * (t - v * x)
    phi1 = gamma * (x - v * t)
    phi2 = y
    phi3 = z

    d_Phi = [sp.diff(Phi, c) for c in (t, x, y, z)]
    d_phi = [[sp.diff(p, c) for c in (t, x, y, z)] for p in (phi1, phi2, phi3)]

    q00, q11, q22, q33 = sp.symbols('q00 q11 q22 q33', real=True)
    q01, q02, q03, q12, q13, q23 = sp.symbols('q01 q02 q03 q12 q13 q23', real=True)
    
    g_inv = sp.Matrix([
        [q00, q01, q02, q03],
        [q01, q11, q12, q13],
        [q02, q12, q22, q23],
        [q03, q13, q23, q33]
    ])

    Y = sp.simplify(sum(g_inv[i,j]*d_Phi[i]*d_Phi[j] for i in range(4) for j in range(4)))
    
    B = sp.zeros(3, 3)
    for A in range(3):
        for B_idx in range(3):
            B[A, B_idx] = sp.simplify(sum(-g_inv[i,j]*d_phi[A][i]*d_phi[B_idx][j] for i in range(4) for j in range(4)))

    I1 = sp.simplify(B.trace())
    I2 = sp.simplify(sp.Rational(1, 2) * (I1**2 - (B*B).trace()))
    I3 = sp.simplify(B.det())

    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols('c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1', real=True)
    L = c_Y*Y + c_Y2*Y**2 + c_I1*I1 + c_I1sq*I1**2 + c_I2*I2 + c_I3*I3 + c_YI1*Y*I1

    # T_01 გამოთვლა (off-diagonal)
    T01 = sp.diff(L, q01) 
    
    subs_minkowski = {
        q00: 1, q11: -1, q22: -1, q33: -1,
        q01: 0, q02: 0, q03: 0, q12: 0, q13: 0, q23: 0
    }
    
    T01_eval = sp.simplify(T01.subs(subs_minkowski))
    
    # T01 უნდა იყოს 0. ვაკვირდებით, რომ ის პროპორციულია (gamma^2 * v)
    # ამოვიღოთ კოეფიციენტი
    lorentz_constraint = sp.simplify(T01_eval / (-2 * gamma**2 * v))
    
    # შევამოწმოთ T11 - T22 (ანიზოტროპია)
    T11 = 2*sp.diff(L, q11) + L
    T22 = 2*sp.diff(L, q22) + L
    T11_eval = sp.simplify(T11.subs(subs_minkowski))
    T22_eval = sp.simplify(T22.subs(subs_minkowski))
    
    anisotropy = sp.simplify(T11_eval - T22_eval)
    aniso_constraint = sp.simplify(anisotropy / (2 * gamma**2 * v**2))
    
    return lorentz_constraint, aniso_constraint

def compare_with_ppn():
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols('c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1', real=True)
    
    # ლოურენცის პირობა T01 = 0-დან
    lorentz_req = c_Y/2 + c_Y2 + c_YI1 - c_I1/2 - 3*c_I1sq - c_I2 - c_I3/2
    lorentz_req = sp.simplify(lorentz_req * 2) # ვაორმაგებთ სიმარტივისთვის
    
    # PPN gamma=1 კონსტრეინტი შემოდის მზის-სისტემის სექტორიდან (p03_solar.py).
    ppn_gamma_req = c_Y + 4*c_Y2 + 2*c_YI1 - c_I1 - 8*c_I1sq - 2*c_I2 - c_I3
    
    diff = sp.simplify(ppn_gamma_req - lorentz_req)
    
    return lorentz_req, ppn_gamma_req, diff


if __name__ == "__main__" and _should_run_main_section("lorentz"):
    print("=" * 72)
    print("PHASE 34: ლოურენცის background-stress შემოწმება")
    print("=" * 72)

    lorentz_constr, aniso_constr = analyze_lorentz_boost()
    
    print("\n1. იმპულსის ნაკადი ბუსტირებულ ნორმალიზებულ ფონში (T_01)")
    print(f"  T_01 ∝ -2 * γ² * v * [ {lorentz_constr} ]")
    print(f"  background-stress სიმეტრიის პირობა T_01 = 0 ითხოვს, რომ ფრჩხილი განულდეს.")
    
    print("\n2. სივრცული ანიზოტროპია ბუსტირებულ ნორმალიზებულ ფონში (T_11 - T_22)")
    print(f"  T_11 - T_22 ∝ 2 * γ² * v² * [ {aniso_constr} ]")
    print("  ანიზოტროპიის განულება იგივე ალგებრულ პირობას ითხოვს.")
    
    print("\n3. PPN γ=1 შედარება (p03_solar.py)")
    lor_req, ppn_req, diff = compare_with_ppn()
    print(f"  Lorentz პირობა: {lor_req} = 0")
    print(f"  PPN γ=1 პირობა: {ppn_req} = 0")
    print(f"  ამ ორი პირობის სხვაობა: {diff} = 0  =>  c_Y2 = c_I1sq")
    
    print("\n4. ფიზიკური დასკვნა")
    print("  RFG-ის ნორმალიზებული სუპერსოლიდური ფონი ზოგად შემთხვევაში არ არის")
    print("  background-stress დონეზე boost-სიმეტრიული.")
    print("  თუ კოეფიციენტები აკმაყოფილებს მიღებულ კონსტრეინტს, ფონის სტრეს-ტენზორი")
    print("  x-boost background-stress დონეზე რჩება T_μν ∝ η_μν.")
    print("  სრული ლოურენცის claim საჭიროებს ყველა boost direction-ისა და perturbation sector-ის audit-ს.")
    print("  დამატებით, p03-ის PPN γ=1 პირობასთან ერთად მიიღება ფაზური და")
    print("  ელასტიური კვადრატული სიხისტეების ტოლობა: c_Y2 = c_I1sq.")


# ===================== OLD BACKBONE INTEGRATION =====================

"""
STAGE A1: OLD backbone -> RFG core candidate ledger

Purpose:
    This block records the valuable mathematical backbone of the old files

        OLD/1. ISPG_FieldEquations.tex
        OLD/2. ISPG_EnergyMomentum.tex
        OLD/3. ISPG_EmergentGeometry.tex
        OLD/11. ISPG_Stability.tex

    inside the new working core.  It is intentionally a candidate ledger:
    the detailed stress, Bianchi, Horndeski/EFT, principal-symbol and Lorentz
    calculations already live above in this same p01_core.py file.

Status:
    Represented as working RFG-core bookkeeping, not closed theory text.
    The Dirac-Bergmann second-class bracket closure, anomaly audit, and
    matter-channel EFT stability remain open.  The DOF count below is a
    candidate target until that proof is written.
"""


def old_variational_backbone_ledger():
    """
    Compact one-metric variational ledger inherited from OLD/1 and OLD/2.

    The key point is the distinction between:
      - direct Euler-Lagrange variation of the scalar: Box(phi)=0 under
        minimal matter coupling;
      - the matter-sourced scalar equation used in the reduced sector, which
        is a trace/bi-conformal consistency relation, not a second independent
        scalar Euler-Lagrange equation.
    """
    G, c = sp.symbols("G c", positive=True)
    kappa = sp.Symbol("kappa_core", real=True)
    box_phi_el, box_phi_trace, T_matter = sp.symbols("Box_phi_EL Box_phi_trace T_matter", real=True)
    dphi_sq = sp.Symbol("(d_phi)^2", real=True)

    return {
        "single_metric": "Matter couples minimally to one physical metric g_mn; no Einstein/Jordan frame split.",
        "matter_stress_definition": "delta S_m = -1/2 int sqrt(-g) T_mn^(m) delta g^mn",
        "scalar_tensor": "T_mn^(phi) = d_m phi d_n phi - 1/2 g_mn (d phi)^2",
        "scalar_trace": sp.Eq(sp.Symbol("T_phi"), -dphi_sq),
        "metric_equation": "G_mn = 8*pi*G*T_mn^(m) + kappa*T_mn^(phi)",
        "kappa": kappa,
        "legacy_negative_kappa": "rejected; RFG core does not use a wrong-sign scalar sector",
        "direct_scalar_EL_unsourced": sp.Eq(box_phi_el, 0),
        "reduced_trace_relation_not_independent_EL": sp.Eq(
            box_phi_trace, -8 * sp.pi * G * T_matter / c**4
        ),
        "divergence_identity": "nabla^m T_mn^(phi) = Box(phi) * d_n phi",
        "conservation_loop": (
            "Bianchi + minimal matter coupling close the system on-shell; "
            "the sourced scalar relation is a reduced-sector consistency relation, "
            "not a simultaneous second Euler-Lagrange equation for the same scalar."
        ),
    }


def old_operational_geometry_ledger():
    """
    Operational metric derivation inherited from OLD/3.

    The old theory's useful content is not a separate postulate set anymore:
    it becomes the RFG interpretation of why the diagonal metric is
    bi-conformal and why gamma_PPN=1 follows structurally.
    """
    phi = sp.Symbol("phi", real=True)
    c, dt, d_sigma = sp.symbols("c dt d_sigma", real=True, positive=True)
    d_tau = sp.exp(phi / 2) * dt
    d_ell = sp.exp(-phi / 2) * d_sigma
    ds2 = sp.exp(phi) * c**2 * dt**2 - sp.exp(-phi) * d_sigma**2

    return {
        "pressure_potential": "phi = log(P_stat/P_max), with phi=0 at asymptotic vacuum normalization.",
        "clock_law": sp.Eq(sp.Symbol("d_tau"), d_tau),
        "rod_law": sp.Eq(sp.Symbol("d_ell"), d_ell),
        "biconformal_metric": sp.Eq(sp.Symbol("ds2"), ds2),
        "weak_field_identification": "Phi_N = Psi = c^2*phi/2",
        "ppn_gamma": sp.Eq(sp.Symbol("gamma_PPN"), 1),
        "stationary_extension": "rotation lives in g_0i/A_i; the diagonal pressure metric is not altered.",
        "one_metric_rule": "No disformal matter metric is introduced; all tests are in the same physical metric.",
    }


def old_constraint_dof_count():
    """
    Dirac-Bergmann bookkeeping inherited from OLD/11.

    The count is a candidate structural stability argument, not yet the full
    bracket closure proof.  It is valuable only as a target ledger until the
    Dirac matrix non-degeneracy/anomaly closure is written.
    """
    phase_dim = 14  # h_ij, pi^ij, phi, p_phi
    first_class = 4  # H and H_i
    second_class_candidate = 4  # two biconformal constraints + two secondary constraints
    dof_unconstrained = sp.Rational(phase_dim - 2 * first_class, 2)
    dof_constrained_candidate = sp.Rational(
        phase_dim - 2 * first_class - second_class_candidate, 2
    )

    return {
        "phase_space_dimension": phase_dim,
        "first_class_constraints": first_class,
        "second_class_constraints_candidate": second_class_candidate,
        "unconstrained_dof": dof_unconstrained,
        "constrained_dof_candidate": dof_constrained_candidate,
        "formula": "N_DOF = (dim Gamma - 2*n_first - n_second)/2",
        "physical_reading": (
            "If the second-class closure is non-degenerate, the biconformal "
            "background removes the independent wrong-sign partner channel."
        ),
        "scope_warning": (
            "Full second-class bracket non-degeneracy/anomaly closure remains "
            "a technical proof target; matter-mediated UV stability still needs "
            "the EFT completion already tracked in phase25."
        ),
    }


def old_stability_gate():
    """
    Consolidated stability bookkeeping and limits.

    This is the logical gate for keeping old stability material in the new
    RFG work file without promoting it to final theory.
    """
    return {
        "hyperbolicity": "covered above by phase24 principal-symbol and well-posedness checks",
        "no_ghost_window": "phase25 covers the phase-sector necessary window; full kinetic-gradient matrix remains required",
        "lorentz_background_stress": "phase34 covers background stress under one boost; perturbation-sector Lorentz audit remains required",
        "constraint_reduction": old_constraint_dof_count(),
        "no_fifth_force": (
            "one-metric minimal coupling means matter follows g_mn geodesics; "
            "there is no independent matter-frame scalar force in the core theory"
        ),
        "energy_caveat": (
            "wrong-sign energy channels are not claimed closed until the "
            "Dirac-bracket closure and EFT matter-channel audit are complete"
        ),
        "old_files_represented_as_candidate_ledgers": [
            "OLD/1. ISPG_FieldEquations.tex",
            "OLD/2. ISPG_EnergyMomentum.tex",
            "OLD/3. ISPG_EmergentGeometry.tex",
            "OLD/11. ISPG_Stability.tex",
        ],
    }


def stage_a1_old_backbone_status():
    return {
        "variational": old_variational_backbone_ledger(),
        "operational_geometry": old_operational_geometry_ledger(),
        "constraint_dof": old_constraint_dof_count(),
        "stability_gate": old_stability_gate(),
        "integration_status": "Stage A.1 represented in p01_core.py as candidate/open-proof bookkeeping",
    }


if __name__ == "__main__" and _should_run_main_section("old"):
    print("=" * 72)
    print("STAGE A1: OLD backbone -> RFG core candidate ledger")
    print("=" * 72)

    status = stage_a1_old_backbone_status()

    print("\n1. Variational backbone")
    for key, value in status["variational"].items():
        print(f"  {key:26s}: {value}")

    print("\n2. Operational geometry")
    for key, value in status["operational_geometry"].items():
        print(f"  {key:26s}: {value}")

    print("\n3. Constraint DOF count")
    for key, value in status["constraint_dof"].items():
        print(f"  {key:26s}: {value}")

    print("\n4. Stability gate")
    gate = status["stability_gate"]
    for key, value in gate.items():
        print(f"  {key:26s}: {value}")

    print("\n5. Candidate ledger verdict")
    print("  - OLD field-equation, stress, emergent-geometry and stability material")
    print("    is represented inside p01_core.py as working bookkeeping, not final theory.")
    print("  - Remaining proof targets: full second-class bracket closure/anomaly audit")
    print("    and matter-channel EFT stability.")


if __name__ == "__main__" and _should_run_main_section("audit"):
    print("=" * 72)
    print("P01 axiom-boundary / proof-gap register")
    print("=" * 72)

    print("\n1. Foundational axiom effective background")
    vacuum_audit = normalized_substrate_background_audit()
    for key, value in vacuum_audit.items():
        if key == "axiom_bridge":
            print(f"  {key:34s}:")
            for bridge_key, bridge_value in value.items():
                print(f"    {bridge_key:32s}: {bridge_value}")
            continue
        print(f"  {key:34s}: {value}")

    print("\n1b. Polynomial-minimum no-go diagnostic")
    derivation_attempt = homogeneous_vacuum_derivation_attempt()
    for key in [
        "phase_stationarity_L_y0",
        "solid_scale_stationarity_L_b0",
        "zero_stress_solution",
        "phase_stationary_solution",
        "K_phi_phase_stationary",
        "K_pi_phase_stationary",
        "obstruction",
        "guardrail",
    ]:
        print(f"  {key:34s}: {derivation_attempt[key]}")

    print("\n2. Local mixed-mode stability gate")
    coeffs_m, _s, _det, _roots = minkowski_principal_symbol()
    mixed_conditions = mixed_mode_stability_conditions(coeffs_m)
    for key, value in mixed_conditions["mixed_polynomial_coefficients"].items():
        print(f"  {key:34s}: {value}")
    print(f"  conditions                        : {mixed_conditions['mixed_speed_required']['conditions']}")
    print(f"  scope                             : {mixed_conditions['scope']}")

    print("\n3. Axiom boundary and remaining proof gaps")
    for item in p01_proof_gap_register():
        print(f"\n  - {item['gap']}")
        print(f"    current_status: {item['current_status']}")
        print(f"    risk          : {item['risk']}")
        print(f"    next_step     : {item['next_step']}")


if __name__ == "__main__" and not _requested_main_sections():
    print("p01_core.py loaded. Run a section explicitly:")
    print("  python RFG\\work\\p01_core.py base|spherical|moduli|stress|horndeski|hyperbolicity|eft|lorentz|old|audit")

