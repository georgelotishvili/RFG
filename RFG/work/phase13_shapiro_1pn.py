# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
Shapiro Time Delay (1PN) Sanity-Check.
ეს ფაილი არ გამოჰყავს gamma-ს, არამედ იყენებს phase8_weak_field.py-ში
მიღებულ შედეგს (gamma=1), რათა შეამოწმოს თავსებადობა Cassini-ს ექსპერიმენტთან.
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

if __name__ == "__main__":
    delta_n, dt_gen, gamma_sym = calculate_shapiro_delay()
    print("--- Shapiro Time Delay (1PN) ---")
    print("ეფექტური რეფრაქციის დანამატი (Delta n):", delta_n)
    print("ზოგადი 1PN დაყოვნება:", dt_gen)
    print("RFG/GR დაყოვნება (gamma = 1 პირობაში):", dt_gen.subs(gamma_sym, 1))

    print("\n--- აგენტთა საბჭოს შენიშვნები / მათემატიკური იდენტობა ---")
    print("1. საზღვრების (x1 და -x0) ჩასმისას ვიღებთ: ln(x1 + sqrt(x1^2+b^2)) - ln(-x0 + sqrt(x0^2+b^2))")
    print("2. მნიშვნელი გარდაიქმნება იდენტობით: (-x0 + sqrt(x0^2+b^2)) = b^2 / (x0 + sqrt(x0^2+b^2))")
    print("3. ეს გვაძლევს ფიზიკურად გამჭვირვალე ფორმას: ln[(x1+l1)(x0+l0)/b^2].")
    print("4. gamma=1 პარამეტრი მოდის phase8_weak_field.py-ს შედეგიდან (Cassini-სთან თავსებადობა - sanity check).")