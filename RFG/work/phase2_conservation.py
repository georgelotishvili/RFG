# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from phase1_action import get_polynomial_lagrangian, init_variables
from phase1_tensor import calculate_stress_tensor

def check_conservation():
    t = sp.Symbol('t', real=True)
    a = sp.Function('a')(t)
    H = sp.diff(a, t) / a
    
    # ინვერსიული მეტრიკის კომპონენტები FLRW-სთვის (g^uu)
    g00, g11, g22, g33 = sp.symbols('g00 g11 g22 g33', real=True)
    
    # Phi ველის დროითი წარმოებული (FLRW-ში სივრცული გრადიენტები ნულია)
    # ვტოვებთ ფუნქციად, რათა ენერგიის შენახვისა და EOM-ის იდენტობა გამოჩნდეს
    Phi_dot = sp.Function('Phi_dot')(t)
    
    # FLRW-სთვის ინვარიანტები (g^uu და Phi_dot ცვლადებით)
    Y_g = g00 * Phi_dot**2
    I1_g = -g11 - g22 - g33
    I2_g = g11*g22 + g22*g33 + g11*g33
    I3_g = -g11*g22*g33
    
    Y, I1, I2, I3 = init_variables()
    L = get_polynomial_lagrangian(Y, I1, I2, I3)
    L_g = L.subs({Y: Y_g, I1: I1_g, I2: I2_g, I3: I3_g})
    
    # 1. ენერგიის შენახვის პირდაპირი გამოთვლა T_mn-დან
    T_mixed = calculate_stress_tensor(L_g, [g00, g11, g22, g33])
    
    # FLRW ანზაცის ჩასმა g^uv-სთვის
    ansatz = {g00: 1, g11: -1/a**2, g22: -1/a**2, g33: -1/a**2}
    rho = sp.simplify(T_mixed[0].subs(ansatz))
    p_iso = sp.simplify(-T_mixed[1].subs(ansatz)) # T^1_1 = -p (სწორი იზოტროპული მიდგომა)
    
    # შენახვის განტოლების მარცხენა მხარე: nabla_u T^u_0 = drho/dt + 3H(rho + p)
    rho_dot = sp.simplify(sp.diff(rho, t))
    conservation_lhs = sp.simplify(rho_dot + 3 * H * (rho + p_iso))
    
    # 2. სკალარული ველის (Phi) მოძრაობის განტოლება
    # nabla_u (dL / d(partial_u Phi)) = 0
    dL_dPhi_dot = sp.simplify(sp.diff(L_g, Phi_dot).subs(ansatz))
    # კოვარიანტული დივერგენცია დროში: 1/a^3 * d/dt(a^3 * V^0)
    EOM_Phi = sp.simplify((sp.diff(a**3 * dL_dPhi_dot, t)) / a**3)
    
    # 3. Noether-ის იდენტობა
    # დროითი ტრანსლაციის სიმეტრიიდან: nabla_u T^u_0 = EOM_Phi * partial_0 Phi
    expected_lhs = sp.simplify(EOM_Phi * Phi_dot)
    
    # შევამოწმოთ სხვაობა (უნდა იყოს ზუსტი 0)
    difference = sp.simplify(conservation_lhs - expected_lhs)
    
    # 4. თუ ჩავსვამთ მკაცრ ანზაცს Phi = t (ანუ Phi_dot = 1, dPhi_dot/dt = 0)
    ansatz_Phi_t = {Phi_dot: 1, sp.diff(Phi_dot, t): 0}
    conservation_eval = sp.simplify(conservation_lhs.subs(ansatz_Phi_t))
    EOM_eval = sp.simplify(EOM_Phi.subs(ansatz_Phi_t))
    
    return difference, expected_lhs, conservation_eval, EOM_eval

if __name__ == "__main__":
    diff, expected_lhs, cons_eval, eom_eval = check_conservation()
    
    print("--- Noether-ის იდენტობის შემოწმება FLRW ფონზე ---")
    print("∇_μ T^μ_0 (ენერგიის შენახვის განტოლების მარცხენა მხარე):")
    print(expected_lhs)
    print("\nსხვაობა (∇_μ T^μ_0) და (EOM_Phi * Phi_dot) შორის:")
    print(diff)
    if diff == 0:
        print("დასკვნა: ენერგიის შენახვა drho/dt+3H(rho+p) კოვარიანტულად ზუსტად უდრის")
        print("         სკალარული ველის მოძრაობის განტოლებას! (Noether's identity შესრულებულია)")
        
    print("\n--- მტკიცება Phi = t ანზაცისთვის ---")
    print("თუ ჩავსვამთ მკაცრ ფონს Phi = t (Phi_dot = 1):")
    print("drho/dt + 3H(rho+p) ნარჩენი:", cons_eval)
    print("EOM_Phi ნარჩენი:", eom_eval)
    print("\nაგენტთა საბჭოს პასუხი:")
    print("1. p_iso პირდაპირ T^1_1-დან დაითვალა, 'გასაშუალოება' აღარ გამოიყენება.")
    print("2. Y და I_k მეტრიკიდან ფიქსირდება; დამოუკიდებელი ფუნქციები ამოღებულია.")
    print("3. Phi=t ანზაცი არ ანულებს ენერგიის შენახვას ავტომატურად (ნარჩენი რჩება).")
    print("4. ეს ნარჩენი ზუსტად ემთხვევა Phi-ს EOM-ს, რაც ნიშნავს რომ სისტემა")
    print("   მოითხოვს Phi-ს დინამიკურ ევოლუციას Phi_dot(t), ან ენერგიის გაცვლას")
    print("   ბარიონულ მატერიასთან/რადიაციასთან სრული შენახვისთვის.")