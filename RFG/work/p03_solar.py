# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
p03: solar-system / weak-field ledger.

Status:
- 1PN geometry check is conditional compatibility, not a full RFG exterior proof.
- Mercury and 1PN Shapiro are sanity checks once gamma=beta=1 is derived.
- 2PN Shapiro/light-bending is a candidate discriminator until the active RFG
  exterior optical index is derived.
- Frame dragging is blocked until the rotating solution and preferred-frame
  PPN parameters alpha_1, alpha_2, alpha_3 are derived.
"""

import sympy as sp
from p01_core import get_polynomial_lagrangian


def coeff_U(expr, U, n):
    """Coefficient of U^n in a symbolic series."""
    return sp.simplify(sp.diff(expr, U, n).subs(U, 0) / sp.factorial(n))

def analyze_ppn():
    r, GM = sp.symbols('r GM', real=True, positive=True)
    U_expr = GM / r
    gamma, beta, a2 = sp.symbols('gamma beta a2', real=True)
    kappa = sp.Symbol('kappa', real=True)
    
    # Standard Schwarzschild-like PPN coordinates
    # B არის g_tt, A არის -g_rr
    A = 1 + 2*gamma*U_expr + a2*U_expr**2
    B = 1 - 2*U_expr + 2*(beta - 1)*U_expr**2

    # აინშტაინის ტენზორი G^mu_nu
    G_tt = -sp.diff(A, r) / (r * A**2) + (1/A - 1)/r**2
    G_rr = sp.diff(B, r) / (r * A * B) + (1/A - 1)/r**2
    G_thth = sp.diff(B, r, 2)/(2*A*B) - sp.diff(B, r)**2/(4*A*B**2) - sp.diff(A, r)*sp.diff(B, r)/(4*A**2*B) + sp.diff(B, r)/(2*r*A*B) - sp.diff(A, r)/(2*r*A**2)

    # U-ზე დაყვანა (GM = U * r)
    U = sp.Symbol('U', real=True, positive=True)
    G_tt_U = sp.simplify(G_tt.subs(GM, U * r))
    G_rr_U = sp.simplify(G_rr.subs(GM, U * r))
    G_thth_U = sp.simplify(G_thth.subs(GM, U * r))

    # Scale by r^2 რათა გახდეს უგანზომილებო პოლინომი U-ში
    G_tt_scaled = sp.simplify(sp.series(G_tt_U * r**2, U, 0, 3).removeO())
    G_rr_scaled = sp.simplify(sp.series(G_rr_U * r**2, U, 0, 3).removeO())
    G_thth_scaled = sp.simplify(sp.series(G_thth_U * r**2, U, 0, 3).removeO())

    # ინვარიანტები (სწორი იდენტიფიკაციით: Y=g^tt=1/B, I1=-g^rr-g^thth-g^phiphi=1/A+2)
    Y = 1/B
    I1 = 2 + 1/A
    I2 = 1 + 2/A
    I3 = 1/A
    
    Y_s, I1_s, I2_s, I3_s = sp.symbols('Y I1 I2 I3', real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    
    L_eval = L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_Y = sp.diff(L_poly, Y_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I1 = sp.diff(L_poly, I1_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I2 = sp.diff(L_poly, I2_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I3 = sp.diff(L_poly, I3_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    
    T_tt = 2 * L_Y / B - L_eval
    T_rr = 2 * (L_I1 / A + 2 * L_I2 / A + L_I3 / A) - L_eval
    T_thth = 2 * (L_I1 + L_I2 * (1 + 1/A) + L_I3 / A) - L_eval

    T_tt_U = sp.simplify(T_tt.subs(GM, U * r))
    T_rr_U = sp.simplify(T_rr.subs(GM, U * r))
    T_thth_U = sp.simplify(T_thth.subs(GM, U * r))
    
    T_tt_series = sp.simplify(sp.series(T_tt_U, U, 0, 3).removeO())
    T_rr_series = sp.simplify(sp.series(T_rr_U, U, 0, 3).removeO())
    T_thth_series = sp.simplify(sp.series(T_thth_U, U, 0, 3).removeO())

    return U, gamma, beta, a2, G_tt_scaled, G_rr_scaled, G_thth_scaled, T_tt_series, T_rr_series, T_thth_series


def ppn_geometry_gate():
    """
    Geometry-only Schwarzschild-like PPN compatibility check.

    This shows what the vacuum Einstein-tensor part would require. It does not
    close the RFG exterior solution until the RFG stress constraints vanish.
    """
    res = analyze_ppn()
    U, gamma, beta, a2, G_tt_s, G_rr_s, G_thth_s, _, _, _ = res

    G_tt_O1 = coeff_U(G_tt_s, U, 1)
    G_rr_O1 = coeff_U(G_rr_s, U, 1)
    G_thth_O1 = coeff_U(G_thth_s, U, 1)
    G_tt_O2 = coeff_U(G_tt_s, U, 2)
    G_rr_O2 = coeff_U(G_rr_s, U, 2)
    G_thth_O2 = coeff_U(G_thth_s, U, 2)

    geometry_subs = {gamma: 1, a2: 4, beta: 1}

    return {
        "status": "CONDITIONAL_GEOMETRIC_VACUUM_COMPATIBILITY",
        "coordinate_system": "Schwarzschild-like / areal-radius PPN",
        "gamma_condition": [
            sp.Eq(G_rr_O1, 0),
            sp.Eq(G_thth_O1, 0),
        ],
        "a2_condition_after_gamma_1": sp.Eq(sp.simplify(G_tt_O2.subs(gamma, 1)), 0),
        "beta_condition_after_gamma_1_a2_4": [
            sp.Eq(sp.simplify(G_rr_O2.subs({gamma: 1, a2: 4})), 0),
            sp.Eq(sp.simplify(G_thth_O2.subs({gamma: 1, a2: 4})), 0),
        ],
        "residuals_after_geometry_conditions": [
            sp.simplify(G_tt_s.subs(geometry_subs)),
            sp.simplify(G_rr_s.subs(geometry_subs)),
            sp.simplify(G_thth_s.subs(geometry_subs)),
        ],
        "not_claimed": "full RFG exterior solution; stress constraints remain open",
    }


def weak_field_stress_constraint_gate():
    """
    RFG stress constraints through O(U^2) on the GR-like geometry branch.

    These are blockers: gamma=beta=1 is not an RFG exterior proof until these
    coefficients are either derived to vanish or shown to be suppressed.
    """
    res = analyze_ppn()
    U, gamma, beta, a2, _, _, _, T_tt_s, T_rr_s, T_thth_s = res
    geometry_subs = {gamma: 1, a2: 4, beta: 1}

    components = {
        "T^t_t": sp.simplify(T_tt_s.subs(geometry_subs)),
        "T^r_r": sp.simplify(T_rr_s.subs(geometry_subs)),
        "T^theta_theta": sp.simplify(T_thth_s.subs(geometry_subs)),
    }
    coeffs = {
        comp: {
            "O(U^0)": coeff_U(expr, U, 0),
            "O(U^1)": coeff_U(expr, U, 1),
            "O(U^2)": coeff_U(expr, U, 2),
        }
        for comp, expr in components.items()
    }
    equations = {
        comp: [sp.Eq(value, 0) for value in orders.values()]
        for comp, orders in coeffs.items()
    }

    variables = sorted(
        list({
            symbol
            for orders in coeffs.values()
            for value in orders.values()
            for symbol in value.free_symbols
        }),
        key=lambda symbol: symbol.name,
    )
    eqs_o0 = [orders["O(U^0)"] for orders in coeffs.values()]
    eqs_o1 = [orders["O(U^1)"] for orders in coeffs.values()]
    eqs_o2 = [orders["O(U^2)"] for orders in coeffs.values()]
    leading_solution = sp.solve(eqs_o0 + eqs_o1, variables, dict=True)
    strict_solution = sp.solve(eqs_o0 + eqs_o1 + eqs_o2, variables, dict=True)
    leading_subs = leading_solution[0] if leading_solution else {}
    o2_after_leading = {
        comp: sp.simplify(orders["O(U^2)"].subs(leading_subs))
        for comp, orders in coeffs.items()
    }

    return {
        "status": "PARTIAL_1PN_STRESS_CLOSURE_WITH_2PN_OBSTRUCTION",
        "branch": "gamma=1, beta=1, a2=4 geometry inserted",
        "stress_series_after_geometry_conditions": components,
        "coefficients_to_vanish": coeffs,
        "equations_to_solve": equations,
        "leading_O0_O1_solution": leading_solution,
        "O2_residual_after_leading_solution": o2_after_leading,
        "strict_O0_O1_O2_solution": strict_solution,
        "warning": (
            "O(U^0) and O(U^1) have a nontrivial closure branch, but adding "
            "O(U^2) forces the trivial all-zero coefficient solution"
        ),
    }


def solar_1pn_closure_branch():
    """
    Nontrivial stress closure through O(U^1) on the GR-like 1PN branch.

    This strengthens the solar ledger: 1PN compatibility has a concrete
    coefficient branch. It does not close 2PN; the remaining residual is the
    next prediction/obstruction target.
    """
    c_Y2, c_YI1 = sp.symbols("c_Y2 c_YI1", real=True)
    branch = {
        sp.Symbol("c_I1", real=True): 4 * c_Y2 + 2 * c_YI1,
        sp.Symbol("c_I1sq", real=True): c_Y2,
        sp.Symbol("c_I2", real=True): -10 * c_Y2 - 3 * c_YI1,
        sp.Symbol("c_I3", real=True): 8 * c_Y2 + 4 * c_YI1,
        sp.Symbol("c_Y", real=True): -4 * c_Y2 - 2 * c_YI1,
    }
    K_phi = sp.simplify(branch[sp.Symbol("c_Y", real=True)] + 6 * c_Y2 + 3 * c_YI1)
    K_pi = sp.simplify(
        -branch[sp.Symbol("c_I1", real=True)]
        - 6 * branch[sp.Symbol("c_I1sq", real=True)]
        - 2 * branch[sp.Symbol("c_I2", real=True)]
        - branch[sp.Symbol("c_I3", real=True)]
        - c_YI1
    )

    return {
        "status": "NONTRIVIAL_1PN_STRESS_CLOSURE_BRANCH",
        "branch": branch,
        "free_parameters": [c_Y2, c_YI1],
        "phase_no_ghost_prefactor": K_phi,
        "solid_no_ghost_prefactor": K_pi,
        "healthy_window_needed": [
            sp.Gt(K_phi, 0),
            sp.Gt(K_pi, 0),
        ],
        "O2_residual_on_this_branch": {
            "T^t_t O(U^2)": 16 * c_Y2,
            "T^r_r O(U^2)": 16 * c_Y2,
            "T^theta_theta O(U^2)": 8 * c_YI1,
        },
        "interpretation": (
            "1PN Solar-System compatibility can be supported by a nontrivial "
            "coefficient branch; exact GR-like 2PN stress-free closure cannot "
            "be kept unless c_Y2=c_YI1=0, which collapses this branch"
        ),
    }

if __name__ == "__main__":
    res = analyze_ppn()
    U, gamma, beta, a2, G_tt_s, G_rr_s, G_thth_s, T_tt_s, T_rr_s, T_thth_s = res
    
    print("--- PPN ექსპანსია და აინშტაინის განტოლებები ---")
    print("გამოყენებულია სტანდარტული სფერული (Schwarzschild-like) PPN კოორდინატები:")
    print("g_tt = B(r) = 1 - 2U + 2(beta-1)U^2")
    print("g_rr = -A(r) = -(1 + 2*gamma*U + a2*U^2)")
    print("g_thth = -r^2\n")

    print("განზომილებათა აცდენა:")
    print("G^mu_nu ტენზორის კომპონენტები შეიცავენ 1/r^2 ფაქტორს. რადგან U = GM/r, 1/r^2 = U^2/(GM)^2.")
    print("ამიტომ G_scaled = G * r^2 იწყება O(U) რიგით, რაც ნიშნავს G ~ U/r^2 ~ 1/r^3.")
    print("T^mu_nu სტრეს-ტენზორი იწყება O(1) და O(U) რიგით, ანუ T ~ 1 + 1/r.")

    geometry_gate = ppn_geometry_gate()
    stress_gate = weak_field_stress_constraint_gate()

    G_tt_O1 = coeff_U(G_tt_s, U, 1)
    G_rr_O1 = coeff_U(G_rr_s, U, 1)
    G_thth_O1 = coeff_U(G_thth_s, U, 1)

    print("\n--- O(U) რიგის გეომეტრიული ნაწილი (პროპორციულია 1/r^3-ის) ---")
    print("G^t_t (O(U)) =", G_tt_O1)
    print("G^r_r (O(U)) =", G_rr_O1)
    print("G^th_th (O(U)) =", G_thth_O1)
    print("რადგან T^mu_nu-ში 1/r^3 წევრები არ არის (ისინი O(U^3)-ზე იწყება), ეს გეომეტრიული")
    print("წევრები დამოუკიდებლად უნდა განულდეს ვაკუუმში:")
    print(f"G^r_r = 0  =>  {G_rr_O1} = 0  =>  gamma = 1")
    print(f"G^th_th = 0 =>  {G_thth_O1} = 0  =>  gamma = 1")

    G_tt_O2 = coeff_U(G_tt_s, U, 2)
    G_rr_O2 = coeff_U(G_rr_s, U, 2)
    G_thth_O2 = coeff_U(G_thth_s, U, 2)

    G_tt_O2_g1 = sp.simplify(G_tt_O2.subs(gamma, 1))
    G_rr_O2_g1 = sp.simplify(G_rr_O2.subs(gamma, 1))
    G_thth_O2_g1 = sp.simplify(G_thth_O2.subs(gamma, 1))

    print("\n--- O(U^2) რიგის გეომეტრიული ნაწილი (gamma=1 ჩასმით, პროპორციულია 1/r^4-ის) ---")
    print("G^t_t (O(U^2)) =", G_tt_O2_g1)
    print("G^r_r (O(U^2)) =", G_rr_O2_g1)
    print("G^th_th (O(U^2)) =", G_thth_O2_g1)

    print("ვაკუუმში 1/r^4 წევრებიც უნდა განულდეს:")
    print(f"G^t_t = 0  =>  {G_tt_O2_g1} = 0  =>  a2 = 4")

    G_rr_O2_g1_a2 = sp.simplify(G_rr_O2_g1.subs(a2, 4))
    G_thth_O2_g1_a2 = sp.simplify(G_thth_O2_g1.subs(a2, 4))
    print(f"G^r_r (a2=4 ჩასმით) = {G_rr_O2_g1_a2}")
    print(f"G^th_th (a2=4 ჩასმით) = {G_thth_O2_g1_a2}")
    print(f"ორივე განტოლება იძლევა თავსებად პირობას: {G_rr_O2_g1_a2} = 0 (ან {G_thth_O2_g1_a2} = 0)  =>  beta = 1")

    print("\n--- სუპერსოლიდის T_mn constraint gate O(U^2)-მდე ---")
    print("სტატუსი:", stress_gate["status"])
    print("ჩასმული გეომეტრიული შტო:", stress_gate["branch"])
    for comp, coeffs in stress_gate["coefficients_to_vanish"].items():
        print(comp, "coefficients:", coeffs)
    print("O(U^0)+O(U^1) nontrivial solution:", stress_gate["leading_O0_O1_solution"])
    print("O(U^2) residual after that solution:", stress_gate["O2_residual_after_leading_solution"])
    print("strict O(U^0)+O(U^1)+O(U^2) solution:", stress_gate["strict_O0_O1_O2_solution"])
    print("დასკვნა:", stress_gate["warning"])
    
    print("\n--- აგენტთა საბჭოს დასკვნები ---")
    print("1. A და B კოდში სწორადაა იდენტიფიცირებული: B=g_tt, A=-g_rr. სტრეინები B^{AB} ითვლება -1/g_ij-ით (სწორია).")
    print("2. A=1/B წინასწარ აღარ იდება. გამოყვანილია დამოუკიდებელი a2 პარამეტრით G^t_t=0 პირობიდან.")
    print("3. beta=1 დგინდება G^r_r=0 და G^th_th=0 განტოლებებიდან a2=4 ჩასმის შემდეგ.")
    print("4. G^th_th კომპონენტიც დაემატა; corrected O(U), O(U^2) კოეფიციენტებით სისტემა უკვე შეიძლება ერთობლივად შემოწმდეს.")
    print("5. კოორდინატები ცხადად გამოცხადდა როგორც Standard Schwarzschild PPN (არა isotropic).")
    print("6. 1PN stress closure არსებობს, მაგრამ 2PN exact-GR stress-free closure ტრივიალიზდება.")
    print("geometry gate:", geometry_gate["status"])
    print("1PN closure branch:", solar_1pn_closure_branch()["status"])

# ===================== CONSOLIDATED PHASE SECTIONS =====================


# ===================== merged from p03_solar.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
დაკვირვებითი sanity-check.
ეს ბლოკი დამოუკიდებლად არ ამტკიცებს RFG მეტრიკას. ის ამოწმებს, რომ
თუ gamma=1 და beta=1 უკვე გამოყვანილია, მერკურის GR/PPN წევრი სწორად ბრუნდება.
"""
import math

def calculate_mercury_precession():
    # გრავიტაციული და ფიზიკური მუდმივები (SI ერთეულებში)
    G_val = 6.67430e-11
    M_sun = 1.98847e30
    c_val = 299792458.0
    
    # მერკურის ორბიტის პარამეტრები
    a_mercury = 57.90905e9   # მეტრი (დიდი ნახევარღერძი)
    e_mercury = 0.205630     # ექსცენტრისიტეტი
    
    # კეპლერის პერიოდი T = 2π√(a³/GM), ფარული ემპირიული შენატანის თავიდან ასაცილებლად
    T_mercury_sec = 2 * math.pi * math.sqrt(a_mercury**3 / (G_val * M_sun))
    T_mercury_days = T_mercury_sec / (24.0 * 3600.0)
    days_per_century = 36525.0
    
    # პირობითი PPN პარამეტრები; მათი RFG-დან გამოყვანა stress gate-ს ელოდება.
    gamma_val = 1.0
    beta_val = 1.0
    
    # თეორიული PPN ფაქტორი: (2 + 2*gamma - beta) / 3
    ppn_factor = (2 + 2 * gamma_val - beta_val) / 3.0
    
    # პრეცესია თითო ორბიტაზე (რადიანებში)
    delta_phi_rad = (6 * math.pi * G_val * M_sun) / (c_val**2 * a_mercury * (1 - e_mercury**2)) * ppn_factor
    
    # გადაყვანა არკწამებში თითო საუკუნეზე
    orbits_per_century = days_per_century / T_mercury_days
    rad_to_arcsec = (180.0 / math.pi) * 3600.0
    precession_arcsec_per_century = delta_phi_rad * orbits_per_century * rad_to_arcsec
    
    return ppn_factor, precession_arcsec_per_century


def mercury_precession_gate():
    """Mercury perihelion gate once gamma=beta=1 is independently derived."""
    ppn_factor, precession = calculate_mercury_precession()
    return {
        "status": "MATCHES_STANDARD_GR_TERM_IF_GAMMA_BETA_1",
        "ppn_factor": ppn_factor,
        "precession_arcsec_per_century": precession,
        "reference_value_arcsec_per_century": 42.98,
        "not_claimed": (
            "independent RFG proof of gamma=beta=1; this is a conditional "
            "solar-system sanity check"
        ),
    }

if __name__ == "__main__":
    ppn_factor, precession = calculate_mercury_precession()
    print("--- მერკურის პერიჰელიონის პრეცესია: conditional sanity-check ---")
    print(f"PPN ფაქტორი ((2 + 2*gamma - beta)/3): {ppn_factor}")
    print(f"გამოთვლილი პრეცესია: {precession:.2f} არკწამი/საუკუნეში")
    print("სტანდარტული GR/PPN მნიშვნელობა: 42.98 არკწამი/საუკუნეში")
    assert abs(precession - 42.98) < 0.1, f"ცდომილება დიდია: პრეცესია = {precession}"
    print("სტატუსი:", mercury_precession_gate()["status"])


# ===================== merged from p03_solar.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
Shapiro Time Delay (1PN) Sanity-Check.
ეს ფაილი არ გამოჰყავს gamma-ს, არამედ იყენებს p03_solar.py-ში
მიღებულ პირობით შედეგს (gamma=1), რათა შეამოწმოს 1PN ფორმა და Cassini gate.
"""
import sympy as sp

def calculate_shapiro_delay():
    x, b, G, M, c_sym = sp.symbols('x b G M c', real=True, positive=True)
    gamma_sym = sp.Symbol('gamma', real=True)

    # მანძილი გრავიტაციულ ცენტრამდე (r) ტრაექტორიის გასწვრივ (სადაც b არის impact parameter)
    r = sp.sqrt(x**2 + b**2)
    U = G * M / (c_sym**2 * r)

    # ეფექტური რეფრაქციული ინდექსი PPN მეტრიკაში სინათლისთვის (1PN მიახლოება)
    # n(r) = c / c_coord ≈ 1 + (1 + gamma) * U
    n_r = 1 + (1 + gamma_sym) * U

    # დაყოვნების ნაწილი: Delta_n = n(r) - 1
    delta_n = (1 + gamma_sym) * U

    # დროის ინტეგრალი dt = (dx / c) * n(x)
    # \Delta t = \int_{-x_0}^{x_1} (delta_n / c) dx
    x0, x1 = sp.symbols('x0 x1', real=True, positive=True)

    delay_integrand = delta_n / c_sym
    
    # ინტეგრაცია (გვაძლევს asinh(x/b), რაც ლოგარითმში გადადის)
    # asinh(x/b) = ln(x/b + sqrt((x/b)^2 + 1)) = ln((x + sqrt(x^2 + b^2))/b)
    integral_res = (1 + gamma_sym) * G * M / c_sym**3 * sp.ln(x + sp.sqrt(x**2 + b**2))
    
    # ინტეგრალის საზღვრები: -x0-დან x1-მდე
    log_term = sp.ln((x1 + sp.sqrt(x1**2 + b**2)) * (x0 + sp.sqrt(x0**2 + b**2)) / b**2)
    coef = (1 + gamma_sym) * G * M / c_sym**3
    # ვიყენებთ Mul(..., evaluate=False) რათა coef არ შევიდეს log-ის შიგნით (base**coef)
    delta_t_general = sp.Mul(coef, log_term, evaluate=False)
    
    return delta_n, delta_t_general, gamma_sym


def cassini_gamma_gate(gamma_value=1.0):
    """
    Cassini 1PN gamma gate.

    Bertotti, Iess, Tortora (Nature 2003): gamma - 1 = (2.1 +/- 2.3)e-5.
    This gate is passed only after gamma is independently derived in the RFG
    exterior branch.
    """
    central = 2.1e-5
    sigma = 2.3e-5
    conservative_bound = abs(central) + sigma
    deviation = abs(gamma_value - 1.0)
    return {
        "status": "PASS_IF_GAMMA_DERIVED" if deviation <= conservative_bound else "FAIL",
        "gamma_value": gamma_value,
        "abs_gamma_minus_1": deviation,
        "cassini_central_gamma_minus_1": central,
        "cassini_sigma": sigma,
        "conservative_bound": conservative_bound,
        "source": "Bertotti, Iess, Tortora, Nature 425, 374-376 (2003)",
        "not_claimed": "Cassini pass before gamma is derived from RFG stress/exterior equations",
    }

if __name__ == "__main__":
    delta_n, dt_gen, gamma_sym = calculate_shapiro_delay()
    print("--- Shapiro Time Delay (1PN) ---")
    print("ეფექტური რეფრაქციის დანამატი (Delta n):", delta_n)
    print("ზოგადი 1PN დაყოვნება:", dt_gen)
    print("RFG/GR დაყოვნება (gamma = 1 პირობაში):", dt_gen.subs(gamma_sym, 1))
    print("Cassini gate:", cassini_gamma_gate())

    print("\n--- აგენტთა საბჭოს შენიშვნები / მათემატიკური იდენტობა ---")
    print("1. საზღვრების (x1 და -x0) ჩასმისას ვიღებთ: ln(x1 + sqrt(x1^2+b^2)) - ln(-x0 + sqrt(x0^2+b^2))")
    print("2. მნიშვნელი გარდაიქმნება იდენტობით: (-x0 + sqrt(x0^2+b^2)) = b^2 / (x0 + sqrt(x0^2+b^2))")
    print("3. ეს გვაძლევს ფიზიკურად გამჭვირვალე ფორმას: ln[(x1+l1)(x0+l0)/b^2].")
    print("4. gamma=1 პარამეტრი პირობითია; Cassini pass ძალაშია მხოლოდ RFG exterior proof-ის შემდეგ.")


# ===================== merged from p03_solar.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 14: 2PN Shapiro and light-bending discriminator.

ძველი ISPG Appendix 5/6-ის კანდიდატი discriminator ახალ RFG ენაზე:
1PN დონეზე RFG და GR ემთხვევა, რადგან gamma=beta=1.
განსხვავება იწყება 2PN რიგში, თუ RFG-ის exterior optical index დაიხურა როგორც

    n_RFG = exp(2 r_g/r),

ხოლო GR-ის isotropic Schwarzschild optical index არის rational ფორმა:

    n_GR = (1+r_g/(2r))^3/(1-r_g/(2r)).

ამ ansatz-ის closed differential შედეგებია:
    Delta t_2PN^(RFG-GR) = (r_g^2/(c b)) * (pi/4),
    Delta theta_RFG = 2 r_s/b + pi*r_s^2/b^2 + O(r_s^3/b^3),
    Delta theta_2PN^RFG / Delta theta_2PN^GR = 16/15.
"""

import sympy as sp


def calculate_shapiro_2pn_discriminator():
    """
    Candidate RFG-GR 2PN Shapiro difference.

    Straight reference path:
        r_0(z)^2 = b^2 + z^2.

    Expansions:
        n_RFG = 1 + 2r_g/r + 2r_g^2/r^2 + ...
        n_GR  = 1 + 2r_g/r + (7/4)r_g^2/r^2 + ...

    Shared 1PN bending/path terms cancel in the RFG-GR difference, so the
    differential coefficient is controlled by alpha_RFG-alpha_GR = 1/4.
    """
    z, b, r_g, c, z_1, z_2 = sp.symbols(
        'z b r_g c z_1 z_2',
        real=True,
        positive=True,
    )
    eps = sp.Symbol('eps', real=True, positive=True)
    r = sp.sqrt(b**2 + z**2)

    n_rfg = sp.exp(2 * eps * r_g / r)
    n_gr = (1 + eps * r_g / (2 * r))**3 / (1 - eps * r_g / (2 * r))

    alpha_rfg = sp.simplify(
        sp.series(n_rfg, eps, 0, 3).coeff(eps, 2) * r**2 / r_g**2
    )
    alpha_gr = sp.simplify(
        sp.series(n_gr, eps, 0, 3).coeff(eps, 2) * r**2 / r_g**2
    )
    delta_alpha = sp.simplify(alpha_rfg - alpha_gr)

    master_integral_infinite = sp.integrate(1 / (b**2 + z**2), (z, -sp.oo, sp.oo))
    master_integral_finite = sp.integrate(1 / (b**2 + z**2), (z, -z_1, z_2))
    delta_t_infinite = sp.simplify(delta_alpha * r_g**2 * master_integral_infinite / c)
    delta_t_finite = sp.simplify(delta_alpha * r_g**2 * master_integral_finite / c)
    dimensionless_delta_b = sp.simplify(delta_t_infinite * c * b / r_g**2)

    return {
        "status": "CANDIDATE_2PN_DISCRIMINATOR",
        "input_status": "BLOCKED_UNTIL_RFG_EXTERIOR_OPTICAL_INDEX_DERIVED",
        "n_RFG": sp.Eq(sp.Symbol('n_RFG'), n_rfg.subs(eps, 1)),
        "n_GR_isotropic": sp.Eq(sp.Symbol('n_GR'), n_gr.subs(eps, 1)),
        "alpha_RFG": alpha_rfg,
        "alpha_GR": alpha_gr,
        "delta_alpha": delta_alpha,
        "master_integral_infinite": master_integral_infinite,
        "master_integral_finite": master_integral_finite,
        "Delta_t_2PN_RFG_minus_GR": sp.Eq(sp.Symbol('Delta_t'), delta_t_infinite),
        "finite_endpoint_Delta_t": sp.Eq(sp.Symbol('Delta_t_finite'), delta_t_finite),
        "Delta_B": sp.Eq(sp.Symbol('Delta_B'), dimensionless_delta_b),
        "candidate_result": "if n_RFG=exp(2*r_g/r), then Delta_B=pi/4",
        "coordinate_warning": (
            "this block is isotropic optical-index language; it must be bridged "
            "to the areal-radius PPN block before theory export"
        ),
    }


def calculate_light_deflection_2pn_discriminator():
    """
    Old Appendix 5 candidate recovered under the exponential optical-index ansatz.

    RFG:
        Delta theta = 2 r_s/b + pi r_s^2/b^2 + ...
    GR:
        Delta theta = 2 r_s/b + (15pi/16) r_s^2/b^2 + ...
    """
    b, r_s = sp.symbols('b r_s', real=True, positive=True)

    theta_1pn = 2 * r_s / b
    theta_2pn_rfg = sp.pi * r_s**2 / b**2
    theta_2pn_gr = sp.Rational(15, 16) * sp.pi * r_s**2 / b**2
    ratio = sp.simplify(theta_2pn_rfg / theta_2pn_gr)
    delta = sp.simplify(theta_2pn_rfg - theta_2pn_gr)

    return {
        "status": "CANDIDATE_2PN_BENDING_DISCRIMINATOR",
        "input_status": "HARD_CODED_UNTIL_RAY_EQUATION_DERIVED_FROM_RFG_METRIC",
        "theta_1PN_shared": theta_1pn,
        "theta_2PN_RFG": theta_2pn_rfg,
        "theta_2PN_GR": theta_2pn_gr,
        "theta_total_RFG": theta_1pn + theta_2pn_rfg,
        "theta_total_GR": theta_1pn + theta_2pn_gr,
        "RFG_over_GR_2PN_ratio": ratio,
        "RFG_2PN_enhancement_percent": sp.N((ratio - 1) * 100, 8),
        "Delta_theta_2PN_RFG_minus_GR": delta,
        "candidate_result": (
            "if the exponential optical-index branch holds, RFG has a 16/15 "
            "enhancement of the GR 2PN bending term"
        ),
    }


if __name__ == "__main__":
    print("--- Shapiro Time Delay (2PN): RFG-GR candidate discriminator ---")
    shapiro = calculate_shapiro_2pn_discriminator()
    for key, value in shapiro.items():
        print(f"{key:34s}: {value}")

    print("\n--- Light deflection (2PN): exponential optical index ---")
    bending = calculate_light_deflection_2pn_discriminator()
    for key, value in bending.items():
        print(f"{key:34s}: {value}")

    print("\n--- აგენტთა საბჭოს შენიშვნები / ტექნიკური შეზღუდვები ---")
    print("1. 1PN დონეზე RFG და GR ემთხვევა: gamma=beta=1.")
    print("2. 2PN Shapiro-ში shared bent-ray/path terms ქრება RFG-GR სხვაობაში.")
    print("3. candidate operational discriminator არის Delta_B=pi/4.")
    print("4. 2PN light bending-ში candidate კოეფიციენტი არის 16/15-ჯერ დიდი GR-ის 2PN წევრზე.")
    print("5. theory export დაბლოკილია, სანამ optical index და null geodesics RFG metric-იდან არ გამოვა.")


# ===================== merged from p03_solar.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 30: Lense-Thirring frame dragging — Gravity Probe B
================================================================================

რეფერენცია: p03_solar.py frame-dragging/PPN სამუშაო ბლოკი

დაკვირვება — Gravity Probe B (Everitt et al. 2011):
- gyroscope on satellite, Earth-orbit
- geodetic precession: 6601.8 ± 18.3 mas/yr (GR: 6606.1)
- frame-dragging (Lense-Thirring): 37.2 ± 7.2 mas/yr (GR: 39.2)

GR Lense-Thirring formula:
Ω_LT = (G/c²r³) · [3(S·r̂)r̂ - S]

RFG-ის ცდა:
- bi-conformal scalar — gravitomagnetic g_0i sector
- phase8 PPN γ=1: geodetic + leading 1.5PN gravitomagnetic sector is conditional
- preferred-frame/vector-PPN proof remains the tightening task
"""

import math


GP_B = {
    "geodetic_precession_obs": 6601.8,  # mas/yr, Everitt et al. 2011
    "geodetic_precession_err": 18.3,
    "geodetic_GR": 6606.1,  # GR prediction
    "Lense_Thirring_obs": 37.2,  # mas/yr
    "Lense_Thirring_err": 7.2,
    "Lense_Thirring_GR": 39.2,
}


def gr_lense_thirring_formula():
    """GR Lense-Thirring formula summary."""
    return {
        "formula": "Ω_LT = (G/c²r³) · [3(S·r̂)r̂ - S]",
        "Earth_S": "Earth angular momentum I·ω",
        "satellite_orbit": "GP-B 642 km altitude polar orbit",
        "GR_prediction_mas_yr": GP_B["Lense_Thirring_GR"],
    }


def rfg_gravitomagnetic_open():
    """RFG bi-conformal gravitomagnetic sector — tightening tasks."""
    return [
        "Leading 1.5PN Lense-Thirring is conditional on the one-metric minimal-coupling GR sector.",
        "Full stationary rotating bi-conformal solution should derive the same g_0i coefficient.",
        "BLOCKER: PPN preferred-frame parameters (α_1, α_2, α_3) must be derived as zero.",
        "MOND rotational bridge must remain inert in the Solar System: Z_rot≈a0/g << 1.",
        "LARES/LARES-2/GINGER are observational context only, not RFG proof.",
    ]


def lageos_lares_comparison():
    """LAGEOS + LARES — improved frame-dragging measurements."""
    return {
        "status": "OBSERVATIONAL_CONTEXT_NOT_USED_AS_RFG_PROOF",
        "LAGEOS_I_II_2011": "reported Lense-Thirring tests; systematics debated",
        "LARES_2016": "reported improved frame-dragging tests; source-tag before export",
        "LARES-2": "targeted/reported high precision; do not use as proof here",
        "GINGER": "Earth-based ring-laser program; future/context only",
    }


def rfg_predictions():
    """RFG-ის ცდა Lense-Thirring-ისთვის."""
    return {
        "PPN_gamma_1PN": "γ=1 branch — geodetic precession matches GR only after stress gate closes",
        "Lense_Thirring_RFG": "conditional leading 1.5PN: Ω_LT = GR under one-metric minimal coupling",
        "MOND_rotational_slot": "Z_rot≈a0/g, so Solar-System correction is <10^-8 to 10^-11",
        "preferred_frame_PPN": "BLOCKED_UNTIL_ALPHA_I_ZERO_DERIVED",
        "current_status": "conditional checklist; rotating RFG solution remains open",
    }


def frame_dragging_gate():
    """Strict gate for Solar-System rotating-sector claims."""
    return {
        "status": "BLOCKED_UNTIL_ROTATING_SOLUTION_AND_ALPHA_I_ZERO_DERIVED",
        "GP_B_reference": GP_B,
        "required_before_claim": [
            "derive stationary rotating RFG exterior g_0i",
            "derive alpha_1=alpha_2=alpha_3=0",
            "show MOND/rotational bridge is inert at Solar-System accelerations",
        ],
        "not_claimed": "RFG independent frame-dragging pass",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 30: Lense-Thirring frame dragging — Gravity Probe B")
    print("რეფერენცია: Everitt 2011, p03_solar frame-dragging/PPN ბლოკი")
    print("=" * 72)

    print("\n1. დაკვირვება (Gravity Probe B)")
    for key, val in GP_B.items():
        print(f"  {key:30s}: {val}")

    print("\n2. GR Lense-Thirring formula")
    gr = gr_lense_thirring_formula()
    for key, val in gr.items():
        print(f"  {key:25s}: {val}")

    print("\n3. RFG-ის გასამკაცრებელი ნაბიჯები")
    for i, task in enumerate(rfg_gravitomagnetic_open(), 1):
        print(f"  {i}. {task}")

    print("\n4. LAGEOS/LARES გადახედვა")
    for key, val in lageos_lares_comparison().items():
        print(f"  {key:25s}: {val}")

    print("\n5. RFG predictions")
    for key, val in rfg_predictions().items():
        print(f"  {key:25s}: {val}")

    print("\n6. სტატუსი")
    print("  - GR L-T 39.2 mas/yr vs GP-B 37.2±7.2 — within 1σ")
    print("  - RFG leading 1.5PN frame-dragging is conditional on one-metric minimal coupling.")
    print("  - preferred-frame α_1, α_2, α_3 derivation is a blocker.")
    print("  - gate:", frame_dragging_gate()["status"])


# =============================================================================
# STAGE C1: OLD solar-system precision gate
# =============================================================================

def stage_c1_old_solar_precision_status():
    """Deletion-gate marker for OLD/4--OLD/8."""
    return {
        "OLD_4_PPN": {
            "status": "migrated_as_conditional_geometry_check",
            "target": "p03_solar.py",
            "core": "Schwarzschild-like weak-field geometry gives gamma=1, beta=1 if RFG stress constraints close",
        },
        "OLD_5_light_deflection": {
            "status": "migrated_as_candidate_discriminator",
            "target": "p03_solar.py PHASE 14",
            "core": "1PN matches GR if gamma=beta=1; 2PN coefficient is candidate until RFG optical index is derived",
        },
        "OLD_6_shapiro": {
            "status": "migrated_as_candidate_discriminator",
            "target": "p03_solar.py PHASE 14",
            "core": "finite RFG-GR 2PN differential discriminator Delta_B=pi/4 under exponential optical-index ansatz",
        },
        "OLD_7_perihelion": {
            "status": "migrated_as_conditional_sanity_check",
            "target": "p03_solar.py Mercury block",
            "core": "PPN factor (2+2gamma-beta)/3=1 gives Mercury 42.98 arcsec/century",
        },
        "OLD_8_frame_dragging": {
            "status": "migrated_with_open_tightening",
            "target": "p03_solar.py PHASE 30",
            "core": "leading 1.5PN Lense-Thirring inherited under one-metric minimal coupling",
            "open": "derive preferred-frame/vector PPN alpha_1=alpha_2=alpha_3=0 explicitly",
        },
        "rotational_MOND_bridge": {
            "status": "kept_as_speculative_and_inert_in_solar_system",
            "suppression": "Z_rot approximately a0/g, so Earth/Sun frame-dragging correction is << current sensitivity",
        },
    }


def stage_c1_solar_falsification_targets():
    """What remains observationally useful after the 1PN GR match."""
    return {
        "closed_1PN": [
            "light bending, conditional on gamma=1",
            "Shapiro logarithmic delay, conditional on gamma=1",
            "perihelion precession, conditional on gamma=beta=1",
            "geodetic precession, conditional on gamma=1",
        ],
        "precision_discriminators": [
            "candidate 2PN Shapiro finite differential Delta_B=pi/4",
            "candidate 2PN bending coefficient enhancement 16/15",
            "strong-field lensing/timing near Sgr A*, pulsar-BH systems, ngEHT/BHEX",
        ],
        "theory_tightening": [
            "solve RFG weak-field stress constraints through O(U^2)",
            "derive active RFG exterior optical index and coordinate bridge",
            "full stationary rotating solution for g_0i",
            "preferred-frame PPN alpha_i proof",
            "nonperturbative null geodesics for compact-object imaging",
        ],
    }


def solar_system_claim_gate():
    """Top-level status ledger for p03."""
    stress_gate = weak_field_stress_constraint_gate()
    one_pn_branch = solar_1pn_closure_branch()
    return {
        "file_export_status": "NOT_READY_FOR_RFG_THEORY_EXPORT",
        "ppn_geometry": ppn_geometry_gate()["status"],
        "rfg_stress_constraints": stress_gate["status"],
        "solar_1PN_closure_branch": one_pn_branch["status"],
        "strict_2PN_stress_free_branch": (
            "TRIVIAL_ONLY" if stress_gate["strict_O0_O1_O2_solution"] else "CHECK"
        ),
        "mercury": mercury_precession_gate()["status"],
        "cassini_shapiro_1PN": cassini_gamma_gate()["status"],
        "shapiro_2PN": calculate_shapiro_2pn_discriminator()["status"],
        "light_bending_2PN": calculate_light_deflection_2pn_discriminator()["status"],
        "frame_dragging": frame_dragging_gate()["status"],
        "preferred_frame": "BLOCKED_UNTIL_ALPHA_I_ZERO_DERIVED",
        "do_not_claim": [
            "do not claim full Solar-System pass",
            "do not claim exact GR-like 2PN stress-free exterior for nonzero coefficients",
            "do not claim 2PN Delta_B=pi/4 as final until optical index is derived",
            "do not claim frame-dragging pass before rotating solution and alpha_i=0",
        ],
    }


if __name__ == "__main__":
    print("\n" + "=" * 72)
    print("p03 top-level claim gate")
    print("=" * 72)
    for key, value in solar_system_claim_gate().items():
        print(f"{key}: {value}")

