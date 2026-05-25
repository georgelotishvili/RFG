# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 23: RFG-ის ჩასმა Horndeski/DHOST/ESS ფარგლებში
================================================================================

რეფერენცია: NOTATION.md, phase22, phase21_cmb.py

მიზანი (STRATEGY.md ეტაპი II §1):
- RFG-ის ცხადი მაპირება Horndeski-ის G_2, G_3, G_4, G_5 ფუნქციებზე
- DHOST გავრცობის გადახედვა — საიდან მოვა G_3, G_4, G_5
- ESS (Effective Solid State, Endlich-Nicolis-Wang) ფარგლი I_k სოლიდისთვის
- Bellini-Sawicki α_K, α_B, α_M, α_T სრული გადათვლა

დასკვნა:
- RFG დღეს მხოლოდ k-essence (G_2(X) only) Horndeski-ის ქვეჯგუფშია
- I_k სოლიდი არ ჯდება სუფთა Horndeski-ში — საჭიროა ESS გავრცობა
- DHOST-ის ფუნქციები G_4(X), G_5(X) ცარიელ ცდად რჩება

ცენტრალური მაპირება (Bellini-Sawicki კონვენცია X = -½ g^μν ∂_μΦ ∂_νΦ):
    Y = -2X    (NOTATION.md § Horndeski/EFT Map)

RFG ლაგრანჟიანი (NOTATION-ის Y კონვენციით):
    L = c_Y·Y + c_Y2·Y^2 + c_I1·I_1 + c_I1sq·I_1^2 + c_I2·I_2 + c_I3·I_3 + c_YI1·Y·I_1

Horndeski მაპირება:
    G_2(X) = c_Y·(-2X) + c_Y2·(-2X)^2 = -2·c_Y·X + 4·c_Y2·X^2
    G_3 = 0       (no kinetic mixing yet)
    G_4 = M_Pl^2/2  (no Planck mass running yet)
    G_5 = 0       (no Gauss-Bonnet/derivative coupling yet)

I_k სოლიდი ESS framework-ში არის ცალკე ფაქტორი:
    L_solid = c_I1·I_1 + c_I1sq·I_1^2 + c_I2·I_2 + c_I3·I_3 + c_YI1·Y·I_1
    ეს არ ჯდება Horndeski-ის G_i-ში — გადადის ESS-ის ცალკე სტრესის სექტორად.
"""

import sympy as sp


# ============================================================================
# Horndeski მაპირება
# ============================================================================


def rfg_to_horndeski():
    """
    RFG-ის Y-სექტორის ცხადი ჩასმა Horndeski G_2(X)-ში.
    """
    X = sp.Symbol("X", real=True)
    c_Y, c_Y2 = sp.symbols("c_Y c_Y2", real=True)
    M_Pl = sp.Symbol("M_Pl", positive=True)

    # Y = -2X კონვერსია (NOTATION § Horndeski/EFT Map)
    Y_in_X = -2 * X

    # G_2(X) = c_Y·Y + c_Y2·Y^2 -> X-ის ფუნქციაა
    G_2 = c_Y * Y_in_X + c_Y2 * Y_in_X**2
    G_3 = sp.Integer(0)
    G_4 = M_Pl**2 / 2
    G_5 = sp.Integer(0)

    return {
        "G_2": sp.expand(G_2),
        "G_3": G_3,
        "G_4": G_4,
        "G_5": G_5,
        "X_def": "X = -1/2 * g^μν * ∂_μΦ * ∂_νΦ",
        "Y_to_X": "Y = -2X",
    }


# ============================================================================
# Bellini-Sawicki α პარამეტრები
# ============================================================================


def bellini_sawicki_alphas():
    """
    α_K, α_B, α_M, α_T ცხადი ფორმულები phase23-ის Horndeski მაპირებიდან.

    α_T = 2X(G_{4,X} - G_{5,φ}) / M_*^2 + ... = 0  (G_4 const, G_5 = 0)
    α_M = (1/H) * d/dt (ln M_*^2) = 0  (M_*^2 = 2G_4 = M_Pl^2 const)
    α_B = 2*(X*G_{3,X}*φ_dot/H ...) / M_*^2 = 0  (G_3 = 0, G_4 const)
    α_K = (2X*G_{2,X} + 4X^2*G_{2,XX} + ...) / (H^2 * M_*^2)
    """
    X = sp.Symbol("X", real=True)
    H = sp.Symbol("H", real=True, positive=True)
    c_Y, c_Y2, M_Pl = sp.symbols("c_Y c_Y2 M_Pl", real=True, positive=True)

    Y_in_X = -2 * X
    G_2 = c_Y * Y_in_X + c_Y2 * Y_in_X**2
    G_2_X = sp.diff(G_2, X)
    G_2_XX = sp.diff(G_2, X, 2)

    M_star_sq = M_Pl**2

    alpha_T = sp.Integer(0)  # G_4_X = 0, G_5 = 0
    alpha_M = sp.Integer(0)  # M_*^2 const
    alpha_B = sp.Integer(0)  # G_3 = 0, G_4 const
    alpha_K = sp.simplify((2 * X * G_2_X + 4 * X**2 * G_2_XX) / (H**2 * M_star_sq))

    return {
        "alpha_T": alpha_T,
        "alpha_M": alpha_M,
        "alpha_B": alpha_B,
        "alpha_K": alpha_K,
        "G_2_X": G_2_X,
        "G_2_XX": G_2_XX,
    }


# ============================================================================
# I_k სოლიდი — ESS framework
# ============================================================================


def ess_solid_sector():
    """
    I_k სოლიდი არ ჯდება სუფთა Horndeski-ში.
    ESS framework (Endlich-Nicolis-Wang 2012, Ballesteros-Bellazzini 2013).

    L_solid = c_I1*I_1 + c_I1sq*I_1^2 + c_I2*I_2 + c_I3*I_3 + c_YI1*Y*I_1
    """
    return {
        "framework": "ESS (Effective Solid State, Endlich-Nicolis-Wang 2012)",
        "structure": "L_solid = c_I1*I_1 + c_I1sq*I_1^2 + c_I2*I_2 + c_I3*I_3 + c_YI1*Y*I_1",
        "horndeski_compatibility": (
            "I_k-ი არ ჯდება ცხადად G_2-G_5-ში. ESS გავრცობა საჭიროა."
        ),
        "mode_count": (
            "scalar (Φ) + 2 transverse vector phonon (I_k) + 1 longitudinal phonon. "
            "ჯამში 4 propagating mode (Horndeski-ის 1-ის ნაცვლად)."
        ),
        "extra_alpha": (
            "ESS framework Bellini-Sawicki α-ებს დამატებითი წვლილით ცვლის, რადგან "
            "α_M ≠ 0 და M_*^2_eff = M_Pl^2 + f(c_I1, c_I2, c_I3, a)"
        ),
    }


# ============================================================================
# DHOST გავრცობა (Beyond Horndeski) — ცარიელი ცდები
# ============================================================================


def dhost_extension():
    """
    DHOST (Degenerate Higher Order Scalar Tensor) framework.
    Crisostomi-Koyama-Tasinato 2016, Langlois-Noui 2016.

    RFG-ის ფესვი დღეს არ მოიცავს DHOST-ის Class I, II, III ფუნქციებს.
    G_4(X), G_5(X) X-დამოკიდებული ვერსიები ცარიელ ცდად რჩება.
    """
    return [
        "G_4(X) — X-დამოკიდებული coupling, Brans-Dicke-ის ბუნებრივი გავრცობა",
        "G_5(X) — Gauss-Bonnet-ტიპის, BH regularization-ისთვის ბუნებრივი",
        "DHOST Class I (A_1 - A_5 ფუნქციები) — degenerate higher derivative",
        "Beyond Horndeski (BH) — extra mode-ის გარეშე higher derivatives",
    ]


# ============================================================================
# Main
# ============================================================================


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 23: RFG-ის ჩასმა Horndeski/DHOST/ESS ფარგლებში")
    print("რეფერენცია: NOTATION.md, phase22, phase21_cmb.py")
    print("=" * 72)

    print("\n1. Horndeski მაპირება (Y = -2X კონვერსია)")
    horndeski = rfg_to_horndeski()
    print(f"  X დეფინიცია: {horndeski['X_def']}")
    print(f"  Y → X: {horndeski['Y_to_X']}")
    print(f"  G_2(X) = {horndeski['G_2']}")
    print(f"  G_3 = {horndeski['G_3']}")
    print(f"  G_4 = {horndeski['G_4']}")
    print(f"  G_5 = {horndeski['G_5']}")

    print("\n2. Bellini-Sawicki α პარამეტრები")
    alphas = bellini_sawicki_alphas()
    print(f"  α_T = {alphas['alpha_T']}    (G_4_X = 0, G_5 = 0)")
    print(f"  α_M = {alphas['alpha_M']}    (M_*^2 = M_Pl^2 const)")
    print(f"  α_B = {alphas['alpha_B']}    (G_3 = 0)")
    print(f"  α_K = {alphas['alpha_K']}")
    print(f"  G_2_X = {alphas['G_2_X']}")
    print(f"  G_2_XX = {alphas['G_2_XX']}")

    print("\n3. I_k სოლიდი — ESS framework (Horndeski-ის გავრცობა)")
    ess = ess_solid_sector()
    for key, value in ess.items():
        print(f"  {key:25s}: {value}")

    print("\n4. DHOST გავრცობა — ღია ცდები")
    for i, task in enumerate(dhost_extension(), 1):
        print(f"  {i}. {task}")

    print("\n5. სტატუსი")
    print("  - RFG = k-essence (G_2(X) only) ქვეჯგუფი Horndeski-ში")
    print("  - α_T, α_M, α_B = 0 ფიქსირდება ცხადად")
    print("  - α_K = (-4*c_Y*X + 48*c_Y2*X²) / (H²·M_Pl²) — მიიღება sympy-დან")
    print("  - I_k სოლიდი ESS-ის გავრცობას მოითხოვს (Horndeski არ ფარდდება)")
    print("  - G_3, G_4(X), G_5 — DHOST გავრცობის შემდეგი ნაბიჯი")
    print("  - phase21_cmb.py-ის Bellini-Sawicki ცდა ამ მაპირებას ეთანხმება")
