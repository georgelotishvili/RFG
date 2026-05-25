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