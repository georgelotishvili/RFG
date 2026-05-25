# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
Shapiro Time Delay (2PN) Consistency Check.
ეს ფაილი წარმოადგენს ლოგიკურ sanity-check-ს phase8_weak_field.py-ში მიღებული 
შედეგების (beta=1, a2=4) საფუძველზე. იგი არ წარმოადგენს სრული 2PN მეტრიკის 
დამოუკიდებელ გამოყვანას საველე განტოლებებიდან.
"""
import sympy as sp

def calculate_shapiro_2pn():
    x, b, rs = sp.symbols('x b rs', real=True, positive=True)
    r = sp.sqrt(x**2 + b**2)
    
    # ექსპანსიის პარამეტრი U = GM/rc^2, rs = 2GM/c^2 -> U = rs / (2r)
    # თუმცა სიმარტივისთვის უშუალოდ rs/r-ით გავშალოთ (რაც პროპორციულია U-სი)
    u = rs / r
    
    # GR მეტრიკა (Schwarzschild 2PN-მდე):
    B_GR = 1 - u
    A_GR = 1 + u + u**2
    
    # RFG მეტრიკა (phase8-დან ვიცით, რომ beta=1 და a2=4, ამიტომ A = 1 + u + u**2):
    # შენიშვნა: A = 1/B ექსპანსია 2PN რიგამდე იძლევა ზუსტად 1 + u + u**2.
    B_RFG = 1 - u
    A_RFG = 1 + u + u**2
    
    # სინათლის კოორდინატული სიჩქარე: dt/dx
    # ds^2 = B dt^2 - A dr^2 - r^2 dphi^2 = 0
    # ტრაექტორიაზე r^2 = x^2 + b^2 -> dr = (x/r)dx
    k = x**2 / r**2
    integrand_GR = (1 / sp.sqrt(B_GR)) * sp.sqrt(A_GR * k + b**2 / r**2)
    integrand_RFG = (1 / sp.sqrt(B_RFG)) * sp.sqrt(A_RFG * k + b**2 / r**2)
    
    # ექსპანსია 2PN (O(rs^2)) რიგამდე
    eps = sp.Symbol('eps', real=True, positive=True)
    int_GR_O2 = sp.simplify(sp.series(integrand_GR.subs(rs, eps * rs), eps, 0, 3).coeff(eps, 2))
    int_RFG_O2 = sp.simplify(sp.series(integrand_RFG.subs(rs, eps * rs), eps, 0, 3).coeff(eps, 2))
    
    diff_integrand = sp.simplify(int_RFG_O2 - int_GR_O2)
    return int_GR_O2, int_RFG_O2, diff_integrand

if __name__ == "__main__":
    term_GR, term_RFG, diff = calculate_shapiro_2pn()
    print("--- Shapiro Time Delay (2PN) ---")
    print("GR 2PN ინტეგრანდი:", term_GR)
    print("RFG 2PN ინტეგრანდი:", term_RFG)
    print("სხვაობა (RFG - GR):", diff)
    print("დასკვნა: ინტეგრანდებს შორის სხვაობა ზუსტად 0-ია 2PN რიგში.")

    print("\n--- აგენტთა საბჭოს შენიშვნები / ტექნიკური შეზღუდვები ---")
    print("1. სტატუსი: ეს არის Consistency check (ტრივიალური a-a=0), რადგან A=1/B და beta=1 ჩასმულია წინასწარ.")
    print("   დამოუკიდებელი proof-ისთვის საჭიროა A(r), B(r)-ის ცალკე ამოხსნა O(u^2)-მდე ველის განტოლებებიდან.")
    print("2. ფიზიკური სიზუსტე: რეალური 2PN ტესტი მოითხოვს bent-ray (სხივის გამრუდების) გათვალისწინებას")
    print("   და სასრულ (finite) საზღვრებს. Straight-line ინტეგრანდი 2PN რიგში ფიზიკურად არასრულია.")
    print("3. BepiColombo ტესტადობა: რადგან RFG ზუსტად ემთხვევა GR-ს 2PN-მდე, მათი გარჩევა მზის სისტემაში")
    print("   გართულებულია. განმასხვავებელი ეფექტები უნდა ვეძებოთ 3PN რიგში ან ძლიერ ველებში (Sgr A*, პულსარები).")