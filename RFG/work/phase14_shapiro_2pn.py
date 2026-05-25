# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 14: 2PN Shapiro and light-bending discriminator.

ძველი ISPG Appendix 5/6-ის პროგნოზი ახალ RFG ენაზე:
1PN დონეზე RFG და GR ემთხვევა, რადგან gamma=beta=1.
განსხვავება იწყება 2PN რიგში, რადგან RFG-ის optical index არის

    n_RFG = exp(2 r_g/r),

ხოლო GR-ის isotropic Schwarzschild optical index არის rational ფორმა:

    n_GR = (1+r_g/(2r))^3/(1-r_g/(2r)).

ამის closed differential შედეგებია:
    Delta t_2PN^(RFG-GR) = (r_g^2/(c b)) * (pi/4),
    Delta theta_RFG = 2 r_s/b + pi*r_s^2/b^2 + O(r_s^3/b^3),
    Delta theta_2PN^RFG / Delta theta_2PN^GR = 16/15.
"""

import sympy as sp


def calculate_shapiro_2pn_discriminator():
    """
    Closed RFG-GR 2PN Shapiro difference.

    Straight reference path:
        r_0(z)^2 = b^2 + z^2.

    Expansions:
        n_RFG = 1 + 2r_g/r + 2r_g^2/r^2 + ...
        n_GR  = 1 + 2r_g/r + (7/4)r_g^2/r^2 + ...

    Shared 1PN bending/path terms cancel in the RFG-GR difference, so the
    differential coefficient is controlled by alpha_RFG-alpha_GR = 1/4.
    """
    z, b, r_g, c, z_1, z_2 = sp.symbols(
        'z b r_g c z_1 z_2',
        real=True,
        positive=True,
    )
    eps = sp.Symbol('eps', real=True, positive=True)
    r = sp.sqrt(b**2 + z**2)

    n_rfg = sp.exp(2 * eps * r_g / r)
    n_gr = (1 + eps * r_g / (2 * r))**3 / (1 - eps * r_g / (2 * r))

    alpha_rfg = sp.simplify(
        sp.series(n_rfg, eps, 0, 3).coeff(eps, 2) * r**2 / r_g**2
    )
    alpha_gr = sp.simplify(
        sp.series(n_gr, eps, 0, 3).coeff(eps, 2) * r**2 / r_g**2
    )
    delta_alpha = sp.simplify(alpha_rfg - alpha_gr)

    master_integral_infinite = sp.integrate(1 / (b**2 + z**2), (z, -sp.oo, sp.oo))
    master_integral_finite = sp.integrate(1 / (b**2 + z**2), (z, -z_1, z_2))
    delta_t_infinite = sp.simplify(delta_alpha * r_g**2 * master_integral_infinite / c)
    delta_t_finite = sp.simplify(delta_alpha * r_g**2 * master_integral_finite / c)
    dimensionless_delta_b = sp.simplify(delta_t_infinite * c * b / r_g**2)

    return {
        "n_RFG": sp.Eq(sp.Symbol('n_RFG'), n_rfg.subs(eps, 1)),
        "n_GR_isotropic": sp.Eq(sp.Symbol('n_GR'), n_gr.subs(eps, 1)),
        "alpha_RFG": alpha_rfg,
        "alpha_GR": alpha_gr,
        "delta_alpha": delta_alpha,
        "master_integral_infinite": master_integral_infinite,
        "master_integral_finite": master_integral_finite,
        "Delta_t_2PN_RFG_minus_GR": sp.Eq(sp.Symbol('Delta_t'), delta_t_infinite),
        "finite_endpoint_Delta_t": sp.Eq(sp.Symbol('Delta_t_finite'), delta_t_finite),
        "Delta_B": sp.Eq(sp.Symbol('Delta_B'), dimensionless_delta_b),
        "prediction": "closed old-theory discriminator recovered: Delta_B=pi/4",
    }


def calculate_light_deflection_2pn_discriminator():
    """
    Old Appendix 5 prediction recovered in RFG.

    RFG:
        Delta theta = 2 r_s/b + pi r_s^2/b^2 + ...
    GR:
        Delta theta = 2 r_s/b + (15pi/16) r_s^2/b^2 + ...
    """
    b, r_s = sp.symbols('b r_s', real=True, positive=True)

    theta_1pn = 2 * r_s / b
    theta_2pn_rfg = sp.pi * r_s**2 / b**2
    theta_2pn_gr = sp.Rational(15, 16) * sp.pi * r_s**2 / b**2
    ratio = sp.simplify(theta_2pn_rfg / theta_2pn_gr)
    delta = sp.simplify(theta_2pn_rfg - theta_2pn_gr)

    return {
        "theta_1PN_shared": theta_1pn,
        "theta_2PN_RFG": theta_2pn_rfg,
        "theta_2PN_GR": theta_2pn_gr,
        "theta_total_RFG": theta_1pn + theta_2pn_rfg,
        "theta_total_GR": theta_1pn + theta_2pn_gr,
        "RFG_over_GR_2PN_ratio": ratio,
        "RFG_2PN_enhancement_percent": sp.N((ratio - 1) * 100, 8),
        "Delta_theta_2PN_RFG_minus_GR": delta,
        "prediction": "RFG has a 16/15 = 6.67% enhancement of the 2PN bending term.",
    }


if __name__ == "__main__":
    print("--- Shapiro Time Delay (2PN): RFG-GR closed discriminator ---")
    shapiro = calculate_shapiro_2pn_discriminator()
    for key, value in shapiro.items():
        print(f"{key:34s}: {value}")

    print("\n--- Light deflection (2PN): exponential optical index ---")
    bending = calculate_light_deflection_2pn_discriminator()
    for key, value in bending.items():
        print(f"{key:34s}: {value}")

    print("\n--- აგენტთა საბჭოს შენიშვნები / ტექნიკური შეზღუდვები ---")
    print("1. 1PN დონეზე RFG და GR ემთხვევა: gamma=beta=1.")
    print("2. 2PN Shapiro-ში shared bent-ray/path terms ქრება RFG-GR სხვაობაში.")
    print("3. დარჩენილი operational discriminator არის Delta_B=pi/4.")
    print("4. 2PN light bending-ში RFG-ის კოეფიციენტი არის 16/15-ჯერ დიდი GR-ის 2PN წევრზე.")
    print("5. პრაქტიკული ტესტის ფანჯარა არის strong-field timing/lensing: pulsar-BH, Sgr A*, ngEHT/BHEX.")
