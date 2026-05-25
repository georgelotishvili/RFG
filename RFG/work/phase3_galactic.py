# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
მოძველებული (deprecated) საწყისი ანიზოტროპიის ესკიზი.
სრული MOND თეორია, a_0, mu(x) და BTFR იხილეთ phase11_mond_metric.py-სა და phase20_bullet_cluster.py-ში.
"""
import sympy as sp
from phase1_tensor import get_pressures

def galactic_acceleration():
    r = sp.Symbol('r', real=True, positive=True)
    
    # get_pressures აბრუნებს ანიზოტროპიულ წნევებს, რომლებიც r-ის ფუნქციებია 
    rho, p_rad, p_tan, delta_p = get_pressures(r)
    
    # 1. ტესტ-ნაწილაკის რეალური აჩქარება მედიუმის შიდა წონასწორობიდან:
    # a_test = (p_rad' - 2*delta_p/r) / (rho + p_rad)
    p_rad_prime = sp.diff(p_rad, r)
    a_test = sp.simplify((p_rad_prime - 2 * delta_p / r) / (rho + p_rad))
    
    # 2. ნიუტონური ძალა ჩვეულებრივი ბარიონული მასიდან 
    G = sp.Symbol('G', real=True, positive=True)
    M_b = sp.Function('M_b')(r)
    g_N = G * M_b / r**2
    
    return a_test, delta_p, g_N

if __name__ == "__main__":
    a_test, dp, g_N = galactic_acceleration()
    print("სფერული გალაქტიკური ანიზოტროპია (delta_p):")
    print(dp)
    print("\nტესტ-ნაწილაკის ფორმალური ანიზოტროპიული აჩქარება (ესკიზური):")
    print(a_test)
    print("\nნიუტონური ნაწილი:")
    print(g_N)
    print("\nშენიშვნა: სრული MOND ამოხსნა chi-ველის მეხსიერებით იხილეთ phase11_mond_metric.py-სა და phase20_bullet_cluster.py-ში.")