# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from phase1_action import init_variables, get_polynomial_lagrangian
from phase1_tensor import calculate_stress_tensor

def get_friedmann_equations():
    # დრო და კოსმოლოგიური მასშტაბური ფაქტორი
    t = sp.Symbol('t', real=True)
    a = sp.Function('a')(t)
    
    # ჰაბლის პარამეტრი: H = \dot{a} / a
    H = sp.diff(a, t) / a
    
    # გრავიტაციული კონსტანტა
    kappa = sp.Symbol('kappa', real=True)
    
    # FLRW სტრეს-ტენზორი მკაცრი ვარიაციით
    g00, g11, g22, g33 = sp.symbols('g00 g11 g22 g33', real=True)
    Y_g = g00
    I1_g = -g11 - g22 - g33
    I2_g = g11*g22 + g22*g33 + g11*g33
    I3_g = -g11*g22*g33
    
    Y, I1, I2, I3 = init_variables()
    L_poly = get_polynomial_lagrangian(Y, I1, I2, I3)
    L_g = L_poly.subs({Y: Y_g, I1: I1_g, I2: I2_g, I3: I3_g})
    
    T_mixed = calculate_stress_tensor(L_g, [g00, g11, g22, g33])
    ansatz = {g00: 1, g11: -1/a**2, g22: -1/a**2, g33: -1/a**2}
    
    rho_solid = sp.simplify(T_mixed[0].subs(ansatz))
    p_iso_solid = sp.simplify(-T_mixed[1].subs(ansatz))
    
    # მატერიის და რადიაციის წვლილები
    rho_m0, rho_rad0 = sp.symbols('rho_m0 rho_rad0', real=True)
    rho_m = rho_m0 / a**3
    rho_rad = rho_rad0 / a**4
    
    p_m = 0
    p_rad_matter = rho_rad / 3
    
    # პირველი ფრიდმანის განტოლება სრული შემადგენლობით
    friedmann1 = sp.Eq(3 * H**2, kappa * (rho_solid + rho_m + rho_rad))
    
    # მეორე ფრიდმანის განტოლება (აჩქარების)
    a_ddot = sp.diff(a, t, t)
    friedmann2 = sp.Eq(2 * a_ddot / a + H**2, -kappa * (p_iso_solid + p_m + p_rad_matter))
    
    return friedmann1, friedmann2, a, t, rho_solid, p_iso_solid

if __name__ == "__main__":
    f1, f2, a, t, rho_solid, p_iso_solid = get_friedmann_equations()
    
    print("პირველი ფრიდმანის განტოლება:", f1)
    print("\nმეორე ფრიდმანის განტოლება:", f2)
    
    print("\n--- სუპერსოლიდის ენერგიის სიმკვრივე FLRW ფონზე ---")
    rho_expanded = sp.expand(rho_solid)
    rho_collected = sp.collect(rho_expanded, a)
    print("rho_solid(a) =", rho_collected)
    
    print("\nშენიშვნა: RFG_Theory.md §3-ში მოყვანილ გამარტივებულ ფორმულას აკლია აქ მიღებული")
    print("          -3*c_I2/a^4 და -c_I3/a^6 წევრები, რომლებიც სრული ვარიაციიდან ჩნდება.")
    
    # ემერჯენტული Λ_eff-ის გამოყოფა: ვუშვებთ, რომ a -> უსასრულობაში
    x = sp.Symbol("x", positive=True)
    Lambda_eff = sp.simplify(sp.limit(rho_solid.subs(a, x), x, sp.oo))
    print("\nეფექტური კოსმოლოგიური მუდმივა (a -> infinity):")
    print("Lambda_eff / kappa =", Lambda_eff)
    
    print("\n--- კოსმოლოგიური ისტორია (დომინანტური წევრები) ---")
    print("1. ადრეული ეპოქა (a -> 0): ჭარბობს a^-4 (რადიაცია) და a^-6 (I3 ელასტიურობა).")
    print("2. შუა ეპოქა (a ~ 1): ჭარბობს a^-3 (მატერია).")
    print("3. გვიანი ეპოქა (a -> infinity): ჭარბობს მუდმივი წევრი Lambda_eff -> აჩქარებული გაფართოება.")
    
    print("\n--- აგენტთა საბჭოს შენიშვნები ---")
    print("- Lambda_eff ცხადად გამოიყო (c_Y + 3*c_Y2), თუმცა თეორია არ ხსნის მის ნატურალურობას")
    print("  (მცირე დაკვირვებად მნიშვნელობას); ექსტრემალური fine-tuning შენარჩუნებულია (იხ. §12).")
    print("- rho_m და rho_rad სწორად სკალირდება a^-3 და a^-4 კანონებით.")
    print("- კოსმოლოგიური სისტემა შეიცავს radiation/matter/late-time acceleration-ისთვის საჭირო სკალირებად წევრებს;")
    print("  სრული დინამიკური ამოხსნა და დაკვირვებითი ფიტი ცალკე ეტაპია.")
    print("- ენერგიის შენახვის განტოლების (drho/dt + 3H(rho+p_iso) = 0) შემოწმება ცალკე გადის")
    print("  phase2_conservation.py-ში, სადაც მკაცრი Phi=t ფონის შემთხვევაში ნარჩენი რჩება.")