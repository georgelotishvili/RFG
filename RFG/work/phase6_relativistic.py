# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 6 (relativistic): FLRW რელატივისტური სტრესი
================================================================================

რეფერენცია: NOTATION.md, phase22_full_stress_tensor.py

ეს ფაილი იყენებს NOTATION.md-ის აქტიურ კონვენციას:
- სიგნატურა (+---)
- B^{AB} = -g^{mu nu} * d_mu phi^A * d_nu phi^B
- T_{mu nu} = 2*dL/dg^{mu nu} - g_{mu nu}*L

ეს არის FLRW ფონური ცდა, რომელიც phase22-ის ბირთვით (Bianchi/Noether
იდენტობა) გადადის. შედეგი — rho(a), p_iso(a) — შედარდება RFG_Theory.md §3-ის
ფორმულას.

FLRW ანზაცი:
    ds^2 = dt^2 - a(t)^2 * (dx^2 + dy^2 + dz^2)
    g^{mu nu} = diag(1, -1/a^2, -1/a^2, -1/a^2)
    Phi = Phi(t), phi^A = (x, y, z)
"""

import sympy as sp
from phase22_full_stress_tensor import evaluate_on_background, reduce_zero


def get_flrw_pressures():
    """
    FLRW ფონური ჩასმა phase22-ის evaluate_on_background-ით.

    აბრუნებს:
    - rho = T_{00}
    - p_iso = T_{11} / a^2
    """
    result = evaluate_on_background("flrw", lagrangian_mode="full")

    a = sp.Function("a")(sp.Symbol("t", real=True))
    rho = sp.simplify(result["T_cov"][0, 0])
    p_iso = sp.simplify(result["T_cov"][1, 1] / a**2)

    return rho, p_iso, a, result


def compare_with_theory(rho):
    """
    RFG_Theory.md §3:57 ფორმულა:
        rho = -3*c_I1/a^2 - 9*c_I1sq/a^4 + c_Y + 3*c_Y2 + 3*c_YI1/a^2

    ეს თეორიის ტექსტში გამოტოვებულია c_I2 და c_I3 წევრებისთვის. phase22 სრულ
    ფორმას ბრუნდება — შესწორება გადასატანია RFG_Theory.md-ში.
    """
    c_Y, c_Y2 = sp.symbols("c_Y c_Y2", real=True)
    c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )

    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    Phi = sp.Function("Phi")(t)
    Phi_dot = sp.diff(Phi, t)

    rho_full_expected = (
        -3 * c_I1 / a**2
        - 9 * c_I1sq / a**4
        - 3 * c_I2 / a**4
        - c_I3 / a**6
        + c_Y * Phi_dot**2
        + 3 * c_Y2 * Phi_dot**4
        + 3 * c_YI1 * Phi_dot**2 / a**2
    )

    rho_md_text = (
        -3 * c_I1 / a**2
        - 9 * c_I1sq / a**4
        + c_Y * Phi_dot**2
        + 3 * c_Y2 * Phi_dot**4
        + 3 * c_YI1 * Phi_dot**2 / a**2
    )

    full_match = sp.simplify(rho - rho_full_expected) == 0
    md_text_match = sp.simplify(rho - rho_md_text) == 0
    missing_terms = sp.simplify(rho_full_expected - rho_md_text)

    return full_match, md_text_match, missing_terms


def check_bianchi_residual(result):
    """phase22-ის Bianchi residual FLRW-ზე უნდა იყოს [0,0,0,0]."""
    residual = [reduce_zero(value) for value in result["residual"]]
    is_ok = all(value == 0 for value in residual)
    return is_ok, residual


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 6 (relativistic): FLRW რელატივისტური სტრესი")
    print("რეფერენცია: NOTATION.md, phase22")
    print("=" * 72)

    rho, p_iso, a, result = get_flrw_pressures()

    print("\n1. FLRW ფონური ჩასმა (NOTATION.md კონვენცია)")
    print(f"  invariants (Y, I1, I2, I3): {result['invariants']}")

    print("\n2. ენერგია-იმპულსის ტენზორი")
    print(f"  rho = T_00 = {sp.expand(rho)}")
    print(f"  p_iso = T_11/a^2 = {sp.expand(p_iso)}")

    print("\n3. შედარება RFG_Theory.md §3:57-ის ფორმულასთან")
    full_match, md_text_match, missing = compare_with_theory(rho)
    print(f"  სრული ფორმა (c_I2 და c_I3 ჩათვლით) ემთხვევა: {full_match}")
    print(f"  ძველი §3 ფორმა (c_I2, c_I3-ის გარეშე) ემთხვევა: {md_text_match}")
    print(f"  RFG_Theory.md §3-ში დასამატებელი წევრები: {sp.expand(missing)}")

    print("\n4. Bianchi residual phase22-ის ცდიდან")
    bianchi_ok, residual = check_bianchi_residual(result)
    print(f"  residual vector: {residual}")
    print(f"  Bianchi/Noether იდენტობა: {'OK' if bianchi_ok else 'CHECK'}")

    print("\n5. სტატუსი")
    print("  - კონვენცია: NOTATION.md-ის აქტიური ფორმა")
    print("  - სიგნატურა: (+---)")
    print("  - phase22-ის Bianchi იდენტობა FLRW-ზე სრულდება")
    print("  - rho-ში c_I2, c_I3 წევრები არსებობს — RFG_Theory.md §3-ის შესწორება")
    print("  - p_iso §3-ის ფორმასთან — შემოწმდეს ცალკე")
