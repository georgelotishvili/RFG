# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 1 (tensor): სუპერსოლიდის სტრესი სფერული ანზაცით
================================================================================

რეფერენცია: NOTATION.md, phase22_full_stress_tensor.py

ეს ფაილი იყენებს NOTATION.md-ის აქტიურ კონვენციას:
- სიგნატურა (+---)
- B^{AB} = -g^{mu nu} * d_mu phi^A * d_nu phi^B
- T_{mu nu} = 2*dL/dg^{mu nu} - g_{mu nu}*L (unique symmetric variables)
- off-diagonal: ფაქტორი 1 (არა 2), რადგან g^{mn}=g^{nm}

phase22 ფარდდება Bianchi/Noether იდენტობას სამ ფონზე (Minkowski, FLRW,
Schwarzschild) falsification-ით. ეს ფაილი იყენებს იმავე კონვენციებს
სფერული ანზაცის Δp გენერაციისთვის (MOND-ის ფორმულა §4-ში).

სფერული სტატიკური ანზაცი:
    ds^2 = B(r)*dt^2 - A(r)*dr^2 - C(r)*dOmega^2
ფაზური ველი: Phi = t (სტატიკური)
ელასტიური ველი: phi^1 = f(r) (რადიალური დეფორმაცია), phi^2 = theta, phi^3 = phi

შედეგი: rho, p_rad, p_tan, Δp = p_tan - p_rad
"""

import sympy as sp
from phase1_action import get_polynomial_lagrangian


def get_spherical_invariants():
    """
    სფერული სტატიკური ანზაცის ინვარიანტები NOTATION.md-ის კონვენციით.

    g_inv დიაგონალური: (1/B, -1/A, -1/C, -1/(C*sin^2(theta)))
    f(r) — რადიალური დეფორმაცია (phi^1 = f(r))
    """
    r = sp.Symbol("r", real=True, positive=True)
    theta = sp.Symbol("theta", real=True, positive=True)
    A = sp.Function("A")(r)
    B = sp.Function("B")(r)
    C = sp.Function("C")(r)
    f = sp.Function("f")(r)
    f_prime = sp.diff(f, r)

    # ფაზური ველი — სტატიკური Phi = t
    # ე.ი. d_0 Phi = 1, d_i Phi = 0
    # Y = g^{00} * (d_0 Phi)^2 = 1/B
    Y = 1 / B

    # ელასტიური ველი — comoving phi^A = (f(r), theta, phi)
    # d_1 phi^1 = f'(r), d_2 phi^2 = 1, d_3 phi^3 = 1 (sin(theta)-ის გარეშე)
    # B^{AB} = -g^{mu nu} d_mu phi^A d_nu phi^B
    # B^{11} = -g^{11} * f'^2 = -(-1/A)*f'^2 = f'^2/A
    # B^{22} = -g^{22} * 1 = -(-1/C) = 1/C
    # B^{33} = -g^{33} * 1 = -(-1/(C*sin^2)) = 1/(C*sin^2)
    Bmat = sp.zeros(3, 3)
    Bmat[0, 0] = f_prime**2 / A
    Bmat[1, 1] = 1 / C
    Bmat[2, 2] = 1 / (C * sp.sin(theta)**2)

    I1 = sp.simplify(Bmat.trace())
    I2 = sp.simplify(sp.Rational(1, 2) * (I1**2 - (Bmat * Bmat).trace()))
    I3 = sp.simplify(Bmat.det())

    return r, theta, A, B, C, f, Y, I1, I2, I3


def get_lagrangian_spherical():
    """L = poly(Y, I1, I2, I3) სფერული ანზაცის ცვლადებში."""
    r, theta, A, B, C, f, Y, I1, I2, I3 = get_spherical_invariants()
    Y_s, I1_s, I2_s, I3_s = sp.symbols("Y I1 I2 I3", real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    L = sp.simplify(L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3}))
    return r, theta, A, B, C, f, L, (Y, I1, I2, I3)


def get_stress_tensor():
    """
    T_{mu nu} გენერალური q-ცვლადებით (NOTATION.md-ის კონვენცია), შემდეგ
    სფერული ჩასმა. T-ის სქელეტი — დიაგონალური მხოლოდ (off-diagonal სფერულ
    სტატიკურ ანზაცზე ნულია).
    """
    r, theta, A, B, C, f, _, _ = get_lagrangian_spherical()
    f_prime = sp.diff(f, r)

    # q = (g^00, g^11, g^22, g^33) სიმბოლური
    q = sp.symbols("q0 q1 q2 q3", real=True, nonzero=True)
    # Y = q0, B^{11} = -q1*f'^2, B^{22} = -q2, B^{33} = -q3
    Y_sym = q[0]
    B11 = -q[1] * f_prime**2
    B22 = -q[2]
    B33 = -q[3]
    I1_sym = B11 + B22 + B33
    I2_sym = sp.Rational(1, 2) * (I1_sym**2 - (B11**2 + B22**2 + B33**2))
    I3_sym = B11 * B22 * B33

    Y_s, I1_s, I2_s, I3_s = sp.symbols("Ys I1s I2s I3s", real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    L_sym = L_poly.subs({Y_s: Y_sym, I1_s: I1_sym, I2_s: I2_sym, I3_s: I3_sym})

    # T_{mu mu} = 2 * dL/dq[mu] - L/q[mu]    (q[mu] = g^{mu mu}, g_{mu mu} = 1/q[mu])
    T_cov_general = [
        2 * sp.diff(L_sym, q[mu]) - L_sym / q[mu]
        for mu in range(4)
    ]

    # სფერული ჩასმა
    subs_sph = {
        q[0]: 1 / B,
        q[1]: -1 / A,
        q[2]: -1 / C,
        q[3]: -1 / (C * sp.sin(theta)**2),
    }

    T_cov = [sp.simplify(expr.subs(subs_sph)) for expr in T_cov_general]

    return r, theta, A, B, C, f, T_cov


def get_pressures():
    """
    rho, p_rad, p_tan, Δp სფერული ანზაციდან NOTATION.md-ის კონვენციით.

    rho = T^{0}_{0} = g^{00} * T_{00} = (1/B) * T_{00}
    p_rad = -T^{1}_{1} = -g^{11} * T_{11} = -(-1/A) * T_{11} = T_{11}/A
    p_tan = -T^{2}_{2} = -g^{22} * T_{22} = T_{22}/C
    Δp = p_tan - p_rad
    """
    r, theta, A, B, C, f, T_cov = get_stress_tensor()

    rho = sp.simplify(T_cov[0] / B)
    p_rad = sp.simplify(T_cov[1] / A)
    p_tan = sp.simplify(T_cov[2] / C)
    delta_p = sp.simplify(p_tan - p_rad)

    return r, theta, A, B, C, f, rho, p_rad, p_tan, delta_p


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 1 (tensor): სუპერსოლიდის სტრესი სფერული ანზაცით")
    print("რეფერენცია: NOTATION.md, phase22")
    print("=" * 72)

    r, theta, A, B, C, f, rho, p_rad, p_tan, delta_p = get_pressures()

    c_I1, c_I1sq, c_YI1, c_I2, c_I3 = sp.symbols(
        "c_I1 c_I1sq c_YI1 c_I2 c_I3", real=True
    )

    print("\n1. ენერგიის სიმკვრივე (rho):")
    print(sp.simplify(rho))

    print("\n2. გრძივი წნევა (p_rad):")
    print(sp.simplify(p_rad))

    print("\n3. განივი წნევა (p_tan):")
    print(sp.simplify(p_tan))

    print("\n4. ანიზოტროპია (Δp = p_tan - p_rad):")
    delta_p_expanded = sp.expand(delta_p)
    print(sp.collect(delta_p_expanded, [c_I1, c_I1sq, c_YI1, c_I2, c_I3]))

    print("\n5. სტატუსი:")
    print("  - კონვენცია: NOTATION.md-ის აქტიური ფორმა (T_mn = 2*dL/dg^mn - g_mn*L)")
    print("  - სიგნატურა: (+---)")
    print("  - phase22-ის Bianchi/Noether იდენტობა იყენებს იმავე კონვენციას")
    print("  - Δp გენერირდება f'(r), A, B, C-ის ფუნქციად — RFG_Theory.md §4-თვის")
