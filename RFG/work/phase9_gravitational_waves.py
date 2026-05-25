# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from phase1_action import get_polynomial_lagrangian

def analyze_gw_full():
    t, z = sp.symbols('t z', real=True)
    h = sp.Function('h')(t, z)
    h_dot = sp.diff(h, t)
    h_z = sp.diff(h, z)
    eps = sp.Symbol('eps', real=True)
    a = sp.Symbol('a', real=True, positive=True)

    # FLRW მეტრიკა TT პერტურბაციით
    g_11 = -a**2 * (1 - eps * h)
    g_22 = -a**2 * (1 + eps * h)
    g_33 = -a**2
    
    # -g დეტერმინანტი
    det_g = a**6 * (1 - eps**2 * h**2)
    sqrt_g = sp.sqrt(det_g)
    
    g_inv_11 = sp.series(1/g_11, eps, 0, 3).removeO()
    g_inv_22 = sp.series(1/g_22, eps, 0, 3).removeO()
    g_inv_33 = 1/g_33

    B11 = -g_inv_11
    B22 = -g_inv_22
    B33 = -g_inv_33

    I1_pert = sp.simplify(B11 + B22 + B33)
    I2_pert = sp.simplify(sp.Rational(1, 2) * (I1_pert**2 - (B11**2 + B22**2 + B33**2)))
    I3_pert = sp.simplify(B11 * B22 * B33)
    Y_pert = 1 

    Y_s, I1_s, I2_s, I3_s = sp.symbols('Y I1 I2 I3', real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    
    L_eval = L_poly.subs({Y_s: Y_pert, I1_s: I1_pert, I2_s: I2_pert, I3_s: I3_pert})
    
    # სრული ლაგრანჟიანის სიმკვრივე (sqrt(-g) * L)
    L_density = sp.series(sqrt_g * L_eval, eps, 0, 3).removeO()
    
    # O(h^2) წევრი (გავყოთ a^3-ზე, რათა სუფთა კოეფიციენტი დაგვრჩეს)
    L_O2 = sp.simplify(L_density.coeff(eps, 2) / a**3)
    
    # შევამოწმოთ კინეტიკური წევრები L_solid-ში
    coeff_h_dot2 = L_O2.coeff(h_dot**2)
    coeff_h_z2 = L_O2.coeff(h_z**2)
    
    # მასის წევრი
    mass_term = sp.simplify(L_O2.coeff(h**2))
    
    return coeff_h_dot2, coeff_h_z2, mass_term

def analyze_polarizations_and_dipole():
    # აგენტთა საბჭოს შენიშვნის მიხედვით:
    # სიმბოლური m_A და m_B დამოუკიდებლობით დიპოლის "ნულობა" ცარიელი ანზაცია.
    # სრული მტკიცებულება მოითხოვს Strong-field სხეულების თვით-გრავიტაციული
    # ენერგიის (Damour-Esposito-Farèse) რეალურ პერტურბაციულ გამოთვლას.
    dipole_radiation = sp.Symbol('dipole_radiation_TBD')
    
    # პოლარიზაციები (ჯერ მხოლოდ qualitative მოლოდინებია)
    tt_modes = 2 # h_+, h_x მოდის აინშტაინ-ჰილბერტის სექტორიდან
    breathing_mode = 'Qualitative Placeholder (TBD)' # არსებობს (კვალი არ ნულდება), მაგრამ ამპლიტუდა ~ d^2(Trace)/dt^2
    longitudinal_mode = 'Qualitative Placeholder (TBD)' # არსებობს ელასტიურობის გამო, მაგრამ c_s < c
    
    return tt_modes, breathing_mode, longitudinal_mode, dipole_radiation

if __name__ == "__main__":
    coeff_h_dot2, coeff_h_z2, mass_term_flrw = analyze_gw_full()
    a = sp.Symbol('a', real=True, positive=True)
    mass_term_mink = sp.simplify(mass_term_flrw.subs(a, 1))
    
    print("--- solid-sector TT check: c_T არ იცვლება ---")
    print(f"L_solid-ის კინეტიკური წევრი (h_dot^2): {coeff_h_dot2}")
    print(f"L_solid-ის გრადიენტული წევრი (h_z^2): {coeff_h_z2}")
    print("დასკვნა: L_solid არ შეიცავს h_dot^2 ან h_z^2 წევრებს.")
    print("c_T^2 = c^2 (1 + delta), სადაც delta = 0 ზუსტად სრულდება!")
    print("\nთუმცა L_solid წარმოქმნის გრავიტონის ეფექტურ მასას (m_g^2 * h^2):")
    print(f"FLRW ფონზე მასის კოეფიციენტი: {mass_term_flrw}")
    print(f"Minkowski ფონზე (a=1): {mass_term_mink}")
    print("\nGW170817-ის და მასიური დისპერსიის ასარიდებლად მოითხოვება ახალი კონსტრეინტი:")
    print(f"{mass_term_mink} = 0")

    print("\n--- აგენტთა საბჭოს შენიშვნები ---")
    print("1. I2/I3-ის გაფართოება გასწორდა; MD ტექსტში შევიდა +1.5*c_I2 + 0.5*c_I3 (ნიშნები გამოსწორდა).")
    print("2. Y_pert=1-ის მიუხედავად, sqrt(-g) ექსპანსია ურევს ფაზურ წევრებს O(h^2) მასაში (სწორია).")
    print("3. დიპოლური რადიაციის ნულობა მონიშნულია როგორც სამუშაო ჰიპოთეზა (TBD), რადგან ის სრულ perturbation/source ანალიზს ითხოვს.")