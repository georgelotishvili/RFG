# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
რეგულარული ცენტრის ანზაცის შემოწმება
სტატუსი: 
ფაილი ამოწმებს, უზრუნველყოფს თუ არა რეგულარული ცენტრალური ანზაცი (A=1+O(r^2))
სასრულ გეომეტრიულ სიმრუდეს და მედიუმის სტრეს-ტენზორს.
ეს არის სასრულ-ენერგიული ბირთვის პირობის შემოწმება და არა შავი ხვრელის 
გლობალური რიცხვობრივი ამოხსნა (matching, horizons, ISCO).
"""
import sympy as sp
from phase1_action import get_polynomial_lagrangian

def analyze_regular_center():
    r = sp.Symbol('r', real=True, positive=True)
    a_2, b_2 = sp.symbols('a_2 b_2', real=True)
    B_0 = sp.Symbol('B_0', real=True, positive=True) # B(0) > 0 აუცილებელია
    
    # რეგულარული ცენტრის ანზაცი (A(0) = 1 აუცილებელია რეგულარულობისთვის)
    A_core = 1 + a_2 * r**2
    B_core = B_0 + b_2 * r**2
    
    # G^mu_nu გეომეტრიული ნაწილები
    G_tt = -sp.diff(A_core, r) / (r * A_core**2) + (1/A_core - 1)/r**2
    G_rr = sp.diff(B_core, r) / (r * A_core * B_core) + (1/A_core - 1)/r**2
    G_thth = sp.diff(B_core, r, 2)/(2*A_core*B_core) - sp.diff(B_core, r)**2/(4*A_core*B_core**2) - sp.diff(A_core, r)*sp.diff(B_core, r)/(4*A_core**2*B_core) + sp.diff(B_core, r)/(2*r*A_core*B_core) - sp.diff(A_core, r)/(2*r*A_core**2)
    
    # რიჩის სკალარი R = -G^\mu_\mu
    R_scalar = -(G_tt + G_rr + 2*G_thth)
    
    # ლიმიტები r -> 0
    G_tt_0 = sp.simplify(sp.limit(G_tt, r, 0))
    G_rr_0 = sp.simplify(sp.limit(G_rr, r, 0))
    G_thth_0 = sp.simplify(sp.limit(G_thth, r, 0))
    R_0 = sp.simplify(sp.limit(R_scalar, r, 0))
    
    # Kretschmann სკალარი (K = R_{abcd}R^{abcd}) რეალური სინგულარობის შესამოწმებლად
    A_p, B_p = sp.diff(A_core, r), sp.diff(B_core, r)
    B_pp = sp.diff(B_p, r)
    term1 = (B_pp/(2*B_core) - B_p**2/(4*B_core**2) - A_p*B_p/(4*A_core*B_core))**2
    term2 = (B_p/(r*B_core))**2
    term3 = (A_p/(r*A_core))**2
    term4 = (1 - 1/A_core)**2 / r**4
    K_scalar = (4 / A_core**2) * term1 + (2 / A_core**2) * term2 + (2 / A_core**2) * term3 + 4 * term4
    K_0 = sp.simplify(sp.limit(K_scalar, r, 0))
    
    # სუპერსოლიდის სტრეს-ტენზორი r -> 0-სას
    # ვუშვებთ ცენტრალური სკალარული მუხტის არარსებობას (Psi'=0)
    Y = 1/B_core
    I1 = 2 + 1/A_core
    I2 = 1 + 2/A_core
    I3 = 1/A_core
    
    Y_s, I1_s, I2_s, I3_s = sp.symbols('Y I1 I2 I3', real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    
    L_eval = L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_Y = sp.diff(L_poly, Y_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I1 = sp.diff(L_poly, I1_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I2 = sp.diff(L_poly, I2_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I3 = sp.diff(L_poly, I3_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    
    T_tt = 2 * L_Y / B_core - L_eval
    T_rr = 2 * (L_I1 / A_core + 2 * L_I2 / A_core + L_I3 / A_core) - L_eval
    
    T_tt_0 = sp.simplify(sp.limit(T_tt, r, 0))
    T_rr_0 = sp.simplify(sp.limit(T_rr, r, 0))
    
    return G_tt_0, G_rr_0, G_thth_0, R_0, K_0, T_tt_0, T_rr_0

if __name__ == "__main__":
    G_t, G_r, G_th, R_0, K_0, T_t, T_r = analyze_regular_center()
    print("--- რეგულარული ცენტრის ანზაცის შემოწმება (r -> 0) ---")
    print(f"G^t_t (გეომეტრიული სიმრუდე ცენტრში) = {G_t}")
    print(f"G^r_r (გეომეტრიული წნევა ცენტრში) = {G_r}")
    print(f"G^th_th (კუთხური სიმრუდე) = {G_th}")
    print(f"R (რიჩის სკალარი ცენტრში) = {R_0}")
    print(f"K (Kretschmann სკალარი ცენტრში) = {K_0}")
    
    print(f"\nT^t_t (ენერგიის სიმკვრივე r=0-ზე) = {T_t}")
    print(f"T^r_r (რადიალური წნევა r=0-ზე) = {T_r}")
    
    print("\n--- აგენტთა საბჭოს დასკვნები ---")
    print("1. რეგულარულობის ანზაცით (A=1+O(r^2), B>0) სიმრუდის ინვარიანტები, მათ შორის")
    print("   უმკაცრესი Kretschmann (K) სკალარი, ცენტრში აბსოლუტურად სასრულია (არა სინგულარული)!")
    print("2. მედიუმის სტრეს-ტენზორი T^t_t და T^r_r ასევე სასრულია.")
    print("3. T^t_t და T^r_r ზოგად შემთხვევაში არ არიან ტოლი, ამიტომ ეს არ არის სუფთა")
    print("   de Sitter-ის ვაკუუმი (w=-1). ტერმინი MD ტექსტში შეიცვალა 'სასრულ-ენერგიული ბირთვით'.")
    print("4. ჰორიზონტების, EHT/LIGO პარამეტრების და ISCO-ს პოვნა მოითხოვს გლობალურ")
    print("   რიცხვობრივ ამოხსნას და matching-ს გარე Schwarzschild-ის მეტრიკასთან.")