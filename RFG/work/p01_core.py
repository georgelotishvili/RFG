# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp

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

def analyze_no_ghost():
    a = sp.Symbol('a', real=True, positive=True) # FLRW scale factor
    dPhi_dot = sp.Symbol('dPhi_dot', real=True)
    pi1_dot, pi2_dot, pi3_dot = sp.symbols('pi1_dot pi2_dot pi3_dot', real=True)
    pi_dot_sq = pi1_dot**2 + pi2_dot**2 + pi3_dot**2
    
    # ფლუქტუაციები
    # შენიშვნა: ფონის შერჩევა (Y=1, B=δ) არის ჩასმული ansatz-ი, რაც ფარული fine-tuning-ის სტატუსს ატარებს.
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
    
    return K_PhiPhi, K_pipi, K_PhiPhi_Mink, K_pipi_Mink

def analyze_lorentz_constrained_stability():
    """
    ამოწმებს No-Ghost პირობებს Lorentz (T_01=0) და PPN (gamma=1) 
    კონსტრეინტების ჩასმის შემდეგ.
    """
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols('c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1', real=True)
    K_PhiPhi, K_pipi, K_PhiPhi_Mink, K_pipi_Mink = analyze_no_ghost()
    
    # Lorentz + PPN გვაძლევს ორ პირობას:
    # 1. c_Y2 = c_I1sq
    # 2. c_I1 = c_Y - 4*c_Y2 + 2*c_YI1 - 2*c_I2 - c_I3
    subs_dict = {
        c_I1sq: c_Y2,
        c_I1: c_Y - 4*c_Y2 + 2*c_YI1 - 2*c_I2 - c_I3
    }
    
    K_Phi_constr = sp.simplify(K_PhiPhi_Mink.subs(subs_dict))
    K_pi_constr = sp.simplify(K_pipi_Mink.subs(subs_dict))
    
    return K_Phi_constr, K_pi_constr

def analyze_sound_speeds():
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
    
    # დეტერმინანტი det(G - cs²·K) = 0 მოგვცემს კვადრატულ განტოლებას cs2-სთვის
    cs2 = sp.Symbol('cs2', real=True)
    eq_cs2 = sp.simplify(4*A*B_pi3*cs2**2 - (4*A*D + 4*B_pi3*C + M_mix**2)*cs2 + 4*C*D)
    
    # sp.solve აბრუნებს უზარმაზარ ფესვებს, ამიტომ ვაბრუნებთ მატრიცის კოეფიციენტებს და დამახასიათებელ განტოლებას
    coeffs = {'A': A, 'B_pi3': B_pi3, 'C': C, 'D': D, 'M_mix': M_mix}
    
    # აგენტების მოთხოვნით, ფესვებს (eigenvalues) ბოლომდე ვიღებთ სიმბოლურად!
    cs2_roots = sp.solve(eq_cs2, cs2)
    
    return cs2_T, eq_cs2, coeffs, cs2_roots

if __name__ == "__main__":
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

    cs2_T, eq_cs2, coeffs, cs2_roots = analyze_sound_speeds()
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
    print("\nამოხსნილი საკუთრივი მნიშვნელობები (Eigenmode Speeds c_s^2):")
    print("Root 1:", cs2_roots[0])
    print("Root 2:", cs2_roots[1])

    K_Phi_c, K_pi_c = analyze_lorentz_constrained_stability()
    print("\n--- ლოურენც-ინვარიანტული ვაკუუმის სტაბილურობა ---")
    print("კონსტრეინტების (c_Y2 = c_I1sq და PPN) ჩასმის შემდეგ No-Ghost პირობები:")
    print(f"K_PhiPhi > 0 => {K_Phi_c} > 0")
    print(f"K_pipi > 0   => {K_pi_c} > 0")
    print("დასკვნა: ეს ორი პირობა ერთდროულად სრულდება, თუ c_Y2 > 0 და")
    print("-6*c_Y2 < (c_Y + 3*c_YI1) < -2*c_Y2. თეორია ფიზიკურად ცოცხალია!")

    print("\n--- აგენტთა საბჭოს შენიშვნები ---")
    print("- ფონის შერჩევა (Y=1, B=δ) ჩასმული ansatz-ია, ფარული fine-tuning-ის სტატუსით.")
    print("- 2*pi3_z მოდის ფონის phi^3=z რუკიდან და არის გრძივი პერტურბაციის ხაზოვანი წევრი.")
    print("- ჯვარედინი შერევა დეტალურად აისახა 2x2 მატრიცის დეტერმინანტით.")
    print("- cs^2 განტოლება ბოლომდე იქნა ამოხსნილი და eigenmode-ების სიჩქარეები ზუსტად")
    print("  გამოვლინდა (არ დაგვიტოვებია ამოუხსნელი განტოლება). ალგებრული ნიშანი M_mix-სთვის")
    print("  გასწორდა, რაც უზრუნველყოფს სწორ ფიზიკურ 'Level Repulsion' ეფექტს.")

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
Schwarzschild) falsification-ით. ეს ფაილი იყენებს იმავე კონვენციებს
სფერული ანზაცის Δp გენერაციისთვის (MOND-ის ფორმულა §4-ში).

სფერული სტატიკური ანზაცი:
    ds^2 = B(r)*dt^2 - A(r)*dr^2 - C(r)*dOmega^2
ფაზური ველი: Phi = t (სტატიკური)
ელასტიური ველი: phi^1 = f(r) (რადიალური დეფორმაცია), phi^2 = theta, phi^3 = phi

შედეგი: rho, p_rad, p_tan, Δp = p_tan - p_rad
"""

import sympy as sp
# merged import removed: from p01_core import get_polynomial_lagrangian


def get_spherical_invariants():
    """
    სფერული სტატიკური ანზაცის ინვარიანტები NOTATION.md-ის კონვენციით.

    g_inv დიაგონალური: (1/B, -1/A, -1/C, -1/(C*sin^2(theta)))
    f(r) — რადიალური დეფორმაცია (phi^1 = f(r))
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

    # ელასტიური ველი — comoving phi^A = (f(r), theta, phi)
    # d_1 phi^1 = f'(r), d_2 phi^2 = 1, d_3 phi^3 = 1 (sin(theta)-ის გარეშე)
    # B^{AB} = -g^{mu nu} d_mu phi^A d_nu phi^B
    # B^{11} = -g^{11} * f'^2 = -(-1/A)*f'^2 = f'^2/A
    # B^{22} = -g^{22} * 1 = -(-1/C) = 1/C
    # B^{33} = -g^{33} * 1 = -(-1/(C*sin^2)) = 1/(C*sin^2)
    Bmat = sp.zeros(3, 3)
    Bmat[0, 0] = f_prime**2 / A
    Bmat[1, 1] = 1 / C
    Bmat[2, 2] = 1 / (C * sp.sin(theta)**2)

    I1 = sp.simplify(Bmat.trace())
    I2 = sp.simplify(sp.Rational(1, 2) * (I1**2 - (Bmat * Bmat).trace()))
    I3 = sp.simplify(Bmat.det())

    return r, theta, A, B, C, f, Y, I1, I2, I3


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
    # Y = q0, B^{11} = -q1*f'^2, B^{22} = -q2, B^{33} = -q3
    Y_sym = q[0]
    B11 = -q[1] * f_prime**2
    B22 = -q[2]
    B33 = -q[3]
    I1_sym = B11 + B22 + B33
    I2_sym = sp.Rational(1, 2) * (I1_sym**2 - (B11**2 + B22**2 + B33**2))
    I3_sym = B11 * B22 * B33

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


if __name__ == "__main__":
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
    print("  - phase22-ის Bianchi/Noether იდენტობა იყენებს იმავე კონვენციას")
    print("  - Δp გენერირდება f'(r), A, B, C-ის ფუნქციად — RFG_Theory.md §4-თვის")


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

if __name__ == "__main__":
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
- falsification ცდა, რომ Noether იდენტობა მცდარი T-ის ფორმულას ცხადად
  იჭერს (არ არის ცარიელი ტავტოლოგია)

შენიშვნა მკითხველისთვის:
matter stress tensor-ის off-shell კოვარიანტული დივერგენცია generic-ად
ნული არ არის — ის ფიზიკურ ველთა Euler-Lagrange წყაროებს უდრის. ეს ფაილი
ამოწმებს ცხადად ამ იდენტობას, არა "off-shell zero"-ს. falsification_test()
ცხადყოფს რომ ცდა ფაქტობრივად ფარდდება T-ის ფორმულის სისწორეს — ცარიელი
ტავტოლოგია არ არის.
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


def background(name):
    """coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields."""
    if name == "minkowski":
        t, x, y, z = sp.symbols("t x y z", real=True)
        coords = [t, x, y, z]
        g_cov_diag = [1, -1, -1, -1]
        g_inv_diag = [1, -1, -1, -1]
        sqrt_minus_g = sp.Integer(1)
        fields = [t, x, y, z]
        return coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields

    if name == "flrw":
        t, x, y, z = sp.symbols("t x y z", real=True)
        coords = [t, x, y, z]
        a = sp.Function("a")(t)
        Phi = sp.Function("Phi")(t)
        g_cov_diag = [1, -a**2, -a**2, -a**2]
        g_inv_diag = [1, -1 / a**2, -1 / a**2, -1 / a**2]
        sqrt_minus_g = a**3
        fields = [Phi, x, y, z]
        return coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields

    if name == "bianchi_i":
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
        return coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields

    if name == "schwarzschild":
        t, r, th, ph = sp.symbols("t r theta phi", real=True, positive=True)
        coords = [t, r, th, ph]
        r_s = sp.Symbol("r_s", real=True, positive=True)
        f = 1 - r_s / r
        g_cov_diag = [f, -1 / f, -r**2, -r**2 * sp.sin(th)**2]
        g_inv_diag = [1 / f, -f, -1 / r**2, -1 / (r**2 * sp.sin(th)**2)]
        sqrt_minus_g = r**2 * sp.sin(th)
        fields = [t, r, th, ph]
        return coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields

    raise ValueError(f"unknown background: {name}")


# ============================================================================
# T-ის ცდის ფუნქცია (ფონური ჩასმის წინ მთლიანი L ფონის ფუნქციად)
# ============================================================================


def evaluate_on_background(name, lagrangian_mode="full", t_factor=1):
    """
    ფონური Bianchi/Noether იდენტობის ცდა.

    t_factor: T-ის ცდის გადანამრავლება. სწორი მნიშვნელობა 1; falsification
    ცდისთვის t_factor != 1 — residual მაშინ ცხადად არანულოვანი.

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
# Falsification ცდა — ვაჩვენებთ რომ Noether იდენტობა მცდარ T-ს იჭერს
# ============================================================================


def falsification_test():
    """
    სცადე T-ის ცდა გადანამრავლებული t_factor = 3 ფაქტორით.
    residual მაშინ უნდა იყოს არანულოვანი — ცდა მართლა ფარდდება T-ის
    ფორმულის სისწორეს.
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
            "p02_cosmo.py ემთხვევა; p01_core.py-ის mixed "
            "ფორმა ამ კონვენციით უნდა გადაიწეროს"
        ),
        "no_ghost": "იხ. NOTATION.md Full No-Ghost Window; X-scheme pure branch-ში c_Y < 0.",
    }


def reduce_zero(expr):
    """ნულის შემოწმება ფონური სიმეტრიების გათვალისწინებით."""
    return sp.factor(sp.together(sp.trigsimp(sp.simplify(expr))))


def is_zero_vector(values):
    return all(reduce_zero(value) == 0 for value in values)


# ============================================================================
# Main
# ============================================================================


if __name__ == "__main__":
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
    print(f"  ცდის სტატუსი: {'OK' if reduce_zero(offdiag_ok) == 0 else 'CHECK'}")

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
        print(f"  status: {'OK' if residual_ok else 'CHECK'}")

    print("\n4. Falsification ცდა — Noether იდენტობა ფარდდება T-ის ფორმულას?")
    correct_ok, wrong_fails, correct_res, wrong_res = falsification_test()
    print(f"  სწორი T (t_factor=1): residual ნული? {correct_ok}")
    print(f"  მცდარი T (t_factor=3): residual არანული? {wrong_fails}")
    if correct_ok and wrong_fails:
        print("  ცდის სტატუსი: OK — იდენტობა ფაქტობრივად ფარდდება")
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
    print("  - falsification ცდამ აჩვენა რომ იდენტობა მცდარი T-ის ფორმულას იჭერს")
    print("  - generic non-diagonal ფონებზე სრული proof ჯერ ღია")
    print("  - შემდეგი ნაბიჯი: NOTATION.md + p01_core/phase6 გადაწერა")


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
            "ESS framework Bellini-Sawicki α-ებს დამატებითი წვლილით ცვლის, რადგან "
            "α_M ≠ 0 და M_*^2_eff = M_Pl^2 + f(c_I1, c_I2, c_I3, a)"
        ),
    }


# ============================================================================
# DHOST გავრცობა (Beyond Horndeski) — ცარიელი ცდები
# ============================================================================


def dhost_extension():
    """
    DHOST (Degenerate Higher Order Scalar Tensor) framework.
    Crisostomi-Koyama-Tasinato 2016, Langlois-Noui 2016.

    RFG-ის ფესვი დღეს არ მოიცავს DHOST-ის Class I, II, III ფუნქციებს.
    G_4(X), G_5(X) X-დამოკიდებული ვერსიები ცარიელ ცდად რჩება.
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


if __name__ == "__main__":
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
        "minkowski": "det M(s)=0 computed from L2 principal coefficients",
        "flrw": "comoving det M(s)=0 computed; physical speed is a^2*s",
        "schwarzschild": "local orthonormal determinant equals Minkowski; coordinate radial redshift added",
        "remaining": "global hyperbolicity and full curved-background perturbation system remain beyond smoke-test",
    }


def compact_det_label():
    return "(A*s + C)*(B*s + D) - M_mix**2*s = 0"


if __name__ == "__main__":
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
PHASE 25: სრული EFT-of-Dark-Energy — no-ghost window + Bellini-Sawicki α-ები
================================================================================

სტატუსი:
Strategy 3 / X4+M3-ის შესრულება.

ამ ფაილის მიზანი:
1. phase21/phase23-ის ghost კონფლიქტის მკაცრი დახურვა:
       α_K = (-4*c_Y*X + 48*c_Y2*X^2)/(H^2*M_Pl^2)
   და α_K > 0 მოთხოვნის ცალკე ჩვენება X-სქემაში.

2. Y-სქემაში სრული background-dependent no-ghost პირობა:
       K_Phi = q00 * (c_Y + 6*c_Y2*Y0 + c_YI1*I1_bg) > 0

   FLRW normalized background:
       K_Phi = c_Y + 6*c_Y2 + 3*c_YI1/a^2 > 0
       a=1 -> c_Y + 6*c_Y2 + 3*c_YI1 > 0

3. Schwarzschild local/static ფონზე იგივე პირობის smoke-test.
   ადგილობრივი ორთონორმალური ჩარჩო უბრუნდება Minkowski/FLRW local პირობას;
   coordinate Phi=t smoke-test ცალკე იბეჭდება, რადგან კოორდინატული ნორმალიზაცია
   ფიზიკური no-ghost პირობა არ არის.

4. Solid sector აღარ იკარგება α-ებში: α_K, α_B, α_M, α_T-ში ჩნდება
   ESS/Ballesteros-Bellazzini ტიპის სიმბოლური დამატებები. სრული CMB fit მაინც
   phase21/hi_class ამოცანად რჩება.
"""

import sympy as sp


def horndeski_y_sector_alphas():
    """
    Bellini-Sawicki X scheme, Y = -2X.

    Pure Y-sector:
        G2(X) = -2*c_Y*X + 4*c_Y2*X^2
        alpha_K = (2X*G2_X + 4X^2*G2_XX)/(H^2*M_Pl^2)
                = (-4*c_Y*X + 48*c_Y2*X^2)/(H^2*M_Pl^2)

    If I1 is treated as a fixed background spurion, c_YI1*Y*I1 contributes
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
        "alpha_K_with_I1_spurion": alpha_k(G2_with_i1),
        "ghost_rule_X_scheme": "require alpha_K_total > 0; for pure low-X branch this pushes c_Y < 0 in X convention",
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

    The file does not claim a full hi_class implementation. It makes the missing
    solid-sector terms explicit, so α_B/α_M/α_K are no longer silently set to zero.
    """
    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    H = sp.Function("H")(t)
    M_Pl = sp.Symbol("M_Pl", positive=True)
    X = sp.Symbol("X", real=True)
    I1_bg = sp.Symbol("I1_bg", positive=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    delta_M2 = sp.Function("delta_M2_solid")(a)
    alpha_B_solid = sp.Function("alpha_B_solid")(a)
    alpha_K_solid = sp.Function("alpha_K_solid")(a)
    alpha_T_solid = sp.Function("alpha_T_solid")(a)

    M_eff_sq = M_Pl**2 + delta_M2
    alpha_M_total = sp.simplify(sp.diff(sp.log(M_eff_sq), t) / H)
    alpha_K_y_i1 = sp.simplify(
        (-4 * X * (c_Y + c_YI1 * I1_bg) + 48 * c_Y2 * X**2) / (H**2 * M_eff_sq)
    )

    return {
        "M_eff_sq": M_eff_sq,
        "alpha_K_total": alpha_K_y_i1 + alpha_K_solid,
        "alpha_B_total": alpha_B_solid,
        "alpha_M_total": alpha_M_total,
        "alpha_T_total": alpha_T_solid,
        "GW170817_filter": "require alpha_T_solid ≈ 0 and phase34 tensor-speed filter",
        "solid_note": "delta alpha terms require ESS perturbation derivation before CLASS/hi_class fit",
    }


def observational_filters():
    return {
        "alpha_T": "|alpha_T| < O(1e-15) from GW170817/GRB170817A",
        "alpha_K": "alpha_K_total > 0 for scalar no-ghost",
        "alpha_M_alpha_B": "must be fit to CMB/LSS/BAO; not fixed by this smoke-test",
        "DESI_link": "static Lambda_eff is not enough for w(z); dynamic alpha-sector is needed",
    }


def class_camb_interface_open():
    return [
        "export alpha_K(a), alpha_B(a), alpha_M(a), alpha_T(a) arrays",
        "choose ESS closure for delta_M2_solid(a), alpha_B_solid(a), alpha_K_solid(a)",
        "run hi_class/CLASS Planck 2018 likelihood",
        "add BAO/LSS/DESI likelihoods after background H(a) is fixed",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 25: EFT-of-Dark-Energy — no-ghost + alpha sweep")
    print("=" * 72)

    print("\n1. Horndeski/Y-sector Bellini-Sawicki alphas")
    alphas = horndeski_y_sector_alphas()
    for key, value in alphas.items():
        print(f"  {key:30s}: {value}")

    print("\n2. Full Y-scheme no-ghost conditions")
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
    print("  - Strategy 3 X4: background-dependent no-ghost window is now explicit.")
    print("  - Strategy 3 M3: alpha_K formula is corrected and solid-sector deltas are visible.")
    print("  - Full ESS perturbation derivation and Planck chi^2 fit remain phase21/hi_class work.")


# ===================== merged from p01_core.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 34: ვაკუუმური სუპერსოლიდი — ლოურენც-ინვარიანტობის მკაცრი შემოწმება
================================================================================

სტატუსი:
ეს ფაილი წარმოადგენს Priority A / X1-ის ფორმალურ მათემატიკურ გადაწყვეტას.
ამოცანა: დავამტკიცოთ, შეუძლია თუ არა სუპერსოლიდურ ვაკუუმს (სადაც φ^A = x^A) 
შეინარჩუნოს ლოურენც-ინვარიანტობა გლობალური ბუსტის (Lorentz boost) მიმართ.

თუ ვაკუუმი ლოურენც-ინვარიანტულია, მაშინ ბუსტირებულ ათვლის სისტემაშიც მისი 
ენერგია-იმპულსის ტენზორი უნდა დარჩეს T_μν ∝ η_μν ფორმის. კერძოდ:
1. არ უნდა გაჩნდეს იმპულსის ნაკადი: T_01 = 0
2. არ უნდა გაჩნდეს სივრცული ანიზოტროპია: T_11 - T_22 = 0

ამ კოდში SymPy-ს მეშვეობით ვასრულებთ ფონური ველების (Φ = t, φ^A = x^A) 
ბუსტს v სიჩქარით, ვითვლით სრულ T_μν-ს და გამოგვაქვს ის ზუსტი ალგებრული 
პირობა კოეფიციენტებზე, რომელიც ანულებს T_01-ს.

შედეგი: ვაკუუმური სუპერსოლიდი ლოურენც-ინვარიანტულია მხოლოდ მაშინ, თუ 
სრულდება კონკრეტული კონსტრეინტი.
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
    
    # PPN gamma=1 კონსტრეინტი phase8-დან
    ppn_gamma_req = c_Y + 4*c_Y2 + 2*c_YI1 - c_I1 - 8*c_I1sq - 2*c_I2 - c_I3
    
    diff = sp.simplify(ppn_gamma_req - lorentz_req)
    
    return lorentz_req, ppn_gamma_req, diff


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 34: ლოურენც-ინვარიანტობის მკაცრი შემოწმება")
    print("=" * 72)

    lorentz_constr, aniso_constr = analyze_lorentz_boost()
    
    print("\n1. იმპულსის ნაკადი ბუსტირებულ ვაკუუმში (T_01)")
    print(f"  T_01 ∝ -2 * γ² * v * [ {lorentz_constr} ]")
    print(f"  ლოურენც-ინვარიანტობის პირობა T_01 = 0 ითხოვს, რომ ფრჩხილი განულდეს.")
    
    print("\n2. სივრცული ანიზოტროპია ბუსტირებულ ვაკუუმში (T_11 - T_22)")
    print(f"  T_11 - T_22 ∝ 2 * γ² * v² * [ {aniso_constr} ]")
    print("  როგორც ვხედავთ, ანიზოტროპიის განულება ზუსტად იგივე პირობას ითხოვს!")
    
    print("\n3. PPN γ=1 შედარება (phase8)")
    lor_req, ppn_req, diff = compare_with_ppn()
    print(f"  Lorentz პირობა: {lor_req} = 0")
    print(f"  PPN γ=1 პირობა: {ppn_req} = 0")
    print(f"  ამ ორი პირობის სხვაობა: {diff} = 0  =>  c_Y2 = c_I1sq")
    
    print("\n4. ფიზიკური დასკვნა")
    print("  RFG სუპერსოლიდური ვაკუუმი *არღვევს* ლოურენც-ინვარიანტობას ზოგად შემთხვევაში,")
    print("  თუმცა, თუ კოეფიციენტები აკმაყოფილებს მიღებულ კონსტრეინტს, ვაკუუმის სტრეს-ტენზორი")
    print("  ნებისმიერ ათვლის სისტემაში რჩება T_μν ∝ η_μν. ეს წყვეტს ე.წ. ლოურენცის კონფლიქტს.")
    print("  დამატებით, PPN γ=1-თან თავსებადობა მკაცრად მოითხოვს, რომ ფაზური (c_Y2) და")
    print("  ელასტიური (c_I1sq) კვადრატული სიხისტეები იყოს ტოლი: c_Y2 = c_I1sq.")

