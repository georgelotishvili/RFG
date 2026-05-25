# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from phase1_action import init_variables, get_polynomial_lagrangian

def analyze_oscillon():
    r = sp.Symbol('r', real=True, positive=True)
    t = sp.Symbol('t', real=True)
    theta = sp.Symbol('theta', real=True) # theta = omega * t
    omega = sp.Symbol('omega', real=True, positive=True)
    
    # ოსცილონის რადიალური პროფილი და მისი გრადიენტი
    Phi0 = sp.Function('Phi0')(r)
    Phi0_prime = sp.diff(Phi0, r)
    Phi1 = sp.Function('Phi1')(r) # მეორე ჰარმონიკა არაწრფივი შერევისთვის
    Phi1_prime = sp.diff(Phi1, r)
    
    # სკალარული ველი ფონის (t) და ოსცილაციის (მინიმუმ 2 ჰარმონიკით) ჩათვლით
    delta_Phi = Phi0 * sp.sin(theta) + Phi1 * sp.sin(3*theta)
    Phi_total = t + delta_Phi
    
    # წარმოებულები (theta-თი ვაწარმოებთ t-ს ნაცვლად)
    Phi_dot = 1 + omega * (Phi0 * sp.cos(theta) + 3 * Phi1 * sp.cos(3*theta))
    Phi_r = Phi0_prime * sp.sin(theta) + Phi1_prime * sp.sin(3*theta)
    
    # ფაზური ინვარიანტი Y მეტრიკით g^00=1, g^rr=-1 (Minkowski)
    Y_eval = sp.expand(Phi_dot**2 - Phi_r**2)
    
    Y, I1, I2, I3 = init_variables()
    L_poly = get_polynomial_lagrangian(Y, I1, I2, I3)
    
    # ენერგიის სიმკვრივე: T^0_0 = 2 * g^00 * (dPhi/dt)^2 * dL/dY - L
    dL_dY = sp.diff(L_poly, Y)
    rho_expr = 2 * Phi_dot**2 * dL_dY - L_poly
    
    # ვსვამთ ელასტიურ ინვარიანტებს Minkowski ფონზე (I1=3, I2=3, I3=1)
    bg_subs = {I1: 3, I2: 3, I3: 1}
    rho_sub = rho_expr.subs(bg_subs).subs(Y, Y_eval)
    
    # ფონური ენერგია (როცა ოსცილაცია არ გვაქვს)
    rho_bg = rho_sub.subs({Phi0: 0, Phi0_prime: 0, Phi1: 0, Phi1_prime: 0})
    
    # ოსცილონის წმინდა ენერგია 
    rho_pert = sp.expand(rho_sub - rho_bg)
    
    # დროის ერთ პერიოდზე გასაშუალოება
    rho_avg = sp.integrate(rho_pert, (theta, 0, 2*sp.pi)) / (2*sp.pi)
    rho_avg = sp.simplify(rho_avg)
    
    # ენერგიის სრული ინტეგრალი
    E_total = sp.Integral(rho_avg * 4 * sp.pi * r**2, r)
    
    # ვირიალური პირობა: dE/domega = 0 რეზონანსული სიხშირის დასაფიქსირებლად
    virial_integrand = sp.simplify(sp.diff(rho_avg * 4 * sp.pi * r**2, omega))
    
    return Phi_total, Y_eval, rho_avg, E_total, virial_integrand

if __name__ == "__main__":
    Phi_total, Y_eval, rho_avg, E_total, virial_integrand = analyze_oscillon()
    
    print("--- ოსცილონის ანალიზი (Phi ველის ვარიაცია) ---")
    print("\nსრული ველი Phi(t,r):")
    print(Phi_total.subs(sp.Symbol('theta', real=True), sp.Symbol('omega', real=True, positive=True) * sp.Symbol('t', real=True)))
    
    print("\nფაზური ინვარიანტი Y:")
    print(Y_eval.subs(sp.Symbol('theta', real=True), sp.Symbol('omega', real=True, positive=True) * sp.Symbol('t', real=True)))
    
    print("\nგასაშუალოებული წმინდა ენერგიის სიმკვრივე <rho_osc>:")
    # შევაგროვოთ Phi0-ის ხარისხების მიხედვით
    Phi0 = sp.Function('Phi0')(sp.Symbol('r', real=True, positive=True))
    Phi0_prime = sp.diff(Phi0, sp.Symbol('r', real=True, positive=True))
    Phi1 = sp.Function('Phi1')(sp.Symbol('r', real=True, positive=True))
    Phi1_prime = sp.diff(Phi1, sp.Symbol('r', real=True, positive=True))
    print(sp.collect(sp.expand(rho_avg), [Phi0**2, Phi1**2, Phi0_prime**2]))
    
    print("\nსრული ენერგიის ინტეგრალი (E):")
    print(E_total)
    
    print("\nვირიალური პირობის ინტეგრანდი (dE/domega = 0):")
    print(virial_integrand)
    
    print("\n--- აგენტთა საბჭოს დასკვნები ---")
    print("1. ოსცილაცია ხდება რეალურ Phi ველზე და Y არის არაწრფივი შედეგი (Y=1+2Φ̇+Φ̇²-Φ'²).")
    print("2. I1, I2, I3 ინარჩუნებენ Minkowski ფონის (3, 3, 1) წვლილს ენერგიის გამოთვლაში.")
    print("3. სასრული ენერგიისთვის აუცილებელია Phi0(r), Phi1(r) და მათი წარმოებულები საკმარისად სწრაფად")
    print("   ქრებოდნენ უსასრულობაში, ხოლო r=0-ზე პროფილი რეგულარული იყოს.")
    print("4. 2 ჰარმონიკის ჩართვამ (Phi0, Phi1) დაადასტურა, რომ ენერგიაში ჩნდება არაწრფივი ჯვარედინი")
    print("   შერევები. c_Y2>0 დადებითად მოქმედებს quartic წევრებზე, მაგრამ სრული ენერგიის")
    print("   პოზიტიურობა მოითხოვს rho_avg-ის სრული გამოსახულების ანალიზს.")
    print("5. omega-ს ფიქსაციის ფორმალური პირობაა dE/domega = 0; რეალური omega-ს მისაღებად")
    print("   საჭიროა პროფილის ამოხსნა და საზღვრული პირობები.")
    print("6. დასკვნა: ეს ფაილი წარმოადგენს ფორმალური ოსცილონის ენერგიის ფუნქციონალის ესკიზს,")
    print("   და ჯერ არ ამტკიცებს finite-energy მდგრად ოსცილონს სრულად.")
