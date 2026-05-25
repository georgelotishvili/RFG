# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 11: MOND სექტორი — Bekenstein-Milgrom AQUAL ლაგრანჟიანი
================================================================================

რეფერენცია: NOTATION.md, phase22_full_stress_tensor.py, RFG_Theory.md § 4

სტატუსი (გულახდილი):
- ეს ფაილი იყენებს Bekenstein-Milgrom AQUAL ფარგლს χ ველისთვის.
- L_χ-ის სრული ფორმა მითითებულია ნორმალიზაციით და dimensional analysis-ით.
- mu(x) = x/(1+x) არის ფენომენოლოგიური არჩევანი (Famaey-Binney 2005).
- a_0 = cH_0/(2π) რიცხობრივი დამთხვევაა; Λ_eff-დან გამოყვანა ცალკე ნაბიჯია.
- BTFR-ის "გამოყვანა" deep-MOND ანზაცის ჩასმით — ალგებრული შედეგია.
- SPARC ~175 გალაქტიკის χ² fit ცარიელია — ცალკე ცდა (phase33-ის კანდიდატი).

ცენტრალური ფარგლები:
1. χ-ის სრული AQUAL ლაგრანჟიანი:
       L_χ = -(a_0^2 / (8πG)) * F(y),    y = |∇χ|^2 / a_0^2
   Euler-Lagrange:
       ∇ · (μ(x) * ∇χ) = 4πG * ρ_b,    μ(x) = dF/dy(x^2)
2. სფერული ფონური: dF/dy * χ' = G*M_b/r^2 = g_N
3. ნიუტონური ლიმიტი (y >> 1, μ → 1): χ' = g_N
4. Deep MOND (y << 1, μ ≈ x): χ' = √(a_0 * g_N)
5. ეფექტური აჩქარება: g_eff = g_N + g_χ
6. BTFR deep-MOND ლიმიტში: M_b = v^4 / (a_0 * G)
"""

import sympy as sp
import math


def aqual_lagrangian_explicit():
    """
    Bekenstein-Milgrom AQUAL ლაგრანჟიანი ცხადი ნორმალიზაციით.

    L_χ = -(a_0^2 / (8πG)) * F(y),   y = (∇χ)^2 / a_0^2
    """
    a0, G_grav = sp.symbols("a0 G", real=True, positive=True)
    y = sp.Symbol("y", real=True, positive=True)
    F = sp.Function("F")(y)

    # ლაგრანჟიანის ფორმა (ფაქტორი 8πG სტანდარტული Bekenstein-Milgrom-ისთვის)
    L_chi = -(a0**2 / (8 * sp.pi * G_grav)) * F

    # Euler-Lagrange: ∇·(μ(x)·∇χ) = 4πG·ρ_b
    # μ(x) = dF/dy(x^2),  x = √y = |∇χ|/a_0
    mu = sp.Function("mu")(y)
    EL_lhs = "div(mu(|∇χ|/a_0) * ∇χ)"
    EL_rhs = "4 * pi * G * rho_b"

    return L_chi, mu, EL_lhs, EL_rhs


def derive_mu_choice():
    """
    μ(x) ფენომენოლოგიური არჩევანი (Famaey-Binney 2005).
    ეს არ არის RFG-დან გამოყვანა — ფენომენოლოგიური არჩევანია.

    Deep-MOND consistency check (μ → x as x → 0):
    F(y) = ∫ μ(√y) dy with μ = x/(1+x)
    leading order y → 0: F(y) → (2/3) y^(3/2)  — standard AQUAL deep-MOND form.
    """
    x = sp.Symbol("x", real=True, positive=True)
    y = sp.Symbol("y", real=True, positive=True)

    mu_simple = x / (1 + x)
    dF_dy_form = mu_simple.subs(x, sp.sqrt(y))
    F_y = sp.simplify(y - 2 * sp.sqrt(y) + 2 * sp.log(1 + sp.sqrt(y)))

    # ლიმიტი y >> 1 (ნიუტონური)
    mu_newton = sp.limit(mu_simple, x, sp.oo)

    # ლიმიტი x << 1 (deep MOND): μ ≈ x
    mu_deep = sp.series(mu_simple, x, 0, 2).removeO()

    # Deep-MOND F(y) leading order verification
    F_deep_leading = sp.series(F_y, y, 0, 3).removeO()
    F_deep_target = sp.Rational(2, 3) * y ** sp.Rational(3, 2)
    deep_mond_limit = sp.simplify(sp.limit(F_y / F_deep_target, y, 0, dir="+"))

    return (
        mu_simple,
        F_y,
        mu_newton,
        mu_deep,
        F_deep_leading,
        F_deep_target,
        deep_mond_limit,
    )


def spherical_field_equation():
    """
    სფერული ფონური AQUAL განტოლება.
    ∇·(μ(x)·∇χ) = 4πG·ρ_b
    სფერული სიმეტრიით → μ(|χ'|/a_0) · χ' = G·M_b(r) / r^2 = g_N
    """
    a0, g_N, chi_p, G_grav, M_b, r = sp.symbols(
        "a0 g_N chi_p G M_b r", real=True, positive=True
    )

    # ნიუტონური აჩქარება
    g_N_expr = G_grav * M_b / r**2

    # AQUAL ფონური განტოლება μ-ით, სადაც x = chi_p / a_0
    x_val = chi_p / a0
    mu_at_x = x_val / (1 + x_val)
    eq_full = sp.Eq(mu_at_x * chi_p, g_N)

    # Deep MOND (x << 1, μ ≈ x):
    eq_deep = sp.Eq(x_val * chi_p, g_N)
    chi_p_deep_sols = sp.solve(eq_deep, chi_p)
    chi_p_deep = chi_p_deep_sols[0] if chi_p_deep_sols else None

    # ნიუტონური (x >> 1, μ ≈ 1):
    eq_newton = sp.Eq(chi_p, g_N)
    chi_p_newton = sp.solve(eq_newton, chi_p)[0]

    return g_N_expr, eq_full, chi_p_deep, chi_p_newton


def btfr_derivation():
    """
    BTFR (Baryonic Tully-Fisher Relation): M_b ∝ v^4
    Deep MOND ლიმიტში g_chi = √(a_0 g_N), ორბიტა v^2/r = g_chi.

    შენიშვნა: ეს არის ალგებრული ჩასმა, არა RFG-დან გამოყვანა. 
    g_chi = √(a_0·g_N)-ის ჩასმა v^2/r-ში აუცილებლად იძლევა v^4 ∝ M_b.
    """
    v, r, M_b, G_grav, a0 = sp.symbols(
        "v r M_b G a0", real=True, positive=True
    )
    g_N = G_grav * M_b / r**2
    g_chi = sp.sqrt(a0 * g_N)
    g_eff = g_N + g_chi

    # ორბიტა — გარე რეგიონში g_N << g_chi → v^2/r = g_chi
    orbit_eq = sp.Eq(v**2 / r, g_chi)
    M_b_BTFR = sp.solve(orbit_eq, M_b)[0]

    return g_eff, g_chi, M_b_BTFR


def a0_cosmological_scale():
    """
    a_0 = c·H_0 / (2π) რიცხობრივი დამთხვევა.
    NB: ეს არ არის Λ_eff-დან გამოყვანა — ეს არის ემპირიული scaling.
    Λ_eff → a_0 ხიდი ცარიელია; ცალკე ცდა საჭიროა.
    """
    c_si = 299792458.0
    mpc_to_m = 3.086e22
    a0_target = 1.2e-10  # m/s^2 (Milgrom)

    H0_values = {
        "SH0ES 73": 73000.0 / mpc_to_m,
        "Planck 67.4": 67400.0 / mpc_to_m,
        "median 70": 70000.0 / mpc_to_m,
    }

    results = []
    for label, H0 in H0_values.items():
        a0_predicted = c_si * H0 / (2 * math.pi)
        err = abs(a0_predicted - a0_target) / a0_target * 100
        results.append((label, a0_predicted, err))

    return results, a0_target


def mu_derivation_audit():
    """
    Strategy 3 / M5 audit.

    RFG-ის ფესვიდან μ(x)-ის გამოყვანა ნიშნავს არა მხოლოდ ისეთი F(y)-ის პოვნას,
    რომლის derivative იძლევა μ-ს, არამედ χ-მედიუმის მიკროფიზიკური მოქმედებიდან
    ამ F(y)-ის აუცილებლობის ჩვენებას.
    """
    mu_simple, F_y, _, _, _, F_deep_target, deep_mond_limit = derive_mu_choice()
    return {
        "candidate_mu": mu_simple,
        "candidate_F_y": F_y,
        "deep_mond_target": F_deep_target,
        "deep_mond_consistency": deep_mond_limit,
        "rfg_microphysical_derivation": "FAIL",
        "reason": (
            "phase11 reconstructs a valid AQUAL F(y) for the chosen Famaey-Binney "
            "mu, but no RFG action or vortex-memory dynamics forces this function."
        ),
        "status": "phenomenological AQUAL bridge, not an RFG derivation",
    }


def a0_lambda_eff_audit():
    """
    Strategy 3 / M6 audit.

    Dimensional bridge:
        a0 = c H / (2π)

    This can be evaluated with H0 or with a Lambda-only de Sitter scale
    H_Lambda = c sqrt(Lambda/3), but the 2π and the MOND coupling are not
    derived from the polynomial RFG action. Therefore this remains a scaling
    coincidence until a mechanism fixes the coherence length and coupling.
    """
    c_si = 299792458.0
    mpc_to_m = 3.0856775814913673e22
    lambda_obs_m2 = 1.1056e-52
    a0_target = 1.2e-10

    scales = {
        "Planck H0 67.36": 67.36 * 1000.0 / mpc_to_m,
        "SH0ES H0 73.0": 73.0 * 1000.0 / mpc_to_m,
        "Lambda-only de Sitter": c_si * math.sqrt(lambda_obs_m2 / 3.0),
    }

    rows = []
    for label, H_value in scales.items():
        a0_value = c_si * H_value / (2.0 * math.pi)
        rows.append(
            {
                "scale": label,
                "H_s_inv": H_value,
                "a0_m_s2": a0_value,
                "percent_error_vs_Milgrom": abs(a0_value - a0_target) / a0_target * 100.0,
            }
        )

    return {
        "rows": rows,
        "lambda_obs_m2": lambda_obs_m2,
        "verdict": "NO_DERIVATION",
        "reason": (
            "Lambda_eff fixes a cosmological H-scale only after Friedmann dynamics; "
            "it does not by itself derive the MOND interpolation, the coupling to "
            "baryons, or the 2π normalization."
        ),
        "status": "dimensional bridge / numerical coincidence, not a closed RFG mechanism",
    }


def open_tasks():
    """ცარიელი ცდები, რომელიც phase23+-ში გადადის."""
    return [
        "μ(x) RFG-დან გამოყვანა (ფენომენოლოგიური Famaey-Binney არ უდრის RFG-დან გამოყვანას)",
        "a_0 = c·H_0/(2π)-ის მექანიზმი Λ_eff-დან (ნუმეროლოგიური დამთხვევაა ჯერ)",
        "SPARC ~175 გალაქტიკის χ^2 fit (phase33-ის კანდიდატი)",
        "Bullet/Abell 520/El Gordo lensing შედარება (phase20_bullet_cluster.py)",
        "χ ლაგრანჟიანის RFG-დან გამოყვანა (ჯერ ცალკე AQUAL-ის ჩასმაა)",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 11: MOND სექტორი — Bekenstein-Milgrom AQUAL")
    print("რეფერენცია: NOTATION.md, phase22, RFG_Theory.md § 4")
    print("=" * 72)

    L_chi, mu, EL_lhs, EL_rhs = aqual_lagrangian_explicit()
    print("\n1. AQUAL ლაგრანჟიანი ნორმალიზაციით")
    print(f"  L_χ = {L_chi}")
    print(f"  Euler-Lagrange: {EL_lhs} = {EL_rhs}")

    (
        mu_simple,
        F_y,
        mu_newton,
        mu_deep,
        F_deep_leading,
        F_deep_target,
        deep_mond_limit,
    ) = derive_mu_choice()
    print("\n2. μ(x) ფენომენოლოგიური არჩევანი (Famaey-Binney 2005)")
    print(f"  μ(x) = {mu_simple}")
    print(f"  ნიუტონური ლიმიტი (x → ∞): μ → {mu_newton}")
    print(f"  Deep MOND ლიმიტი (x → 0): μ ≈ {mu_deep}")
    print(f"  რეკონსტრუირებული F(y) = {F_y}")
    print(f"  F(y) deep-MOND ლიდინგ რიგი (y → 0): {F_deep_leading}")
    print(f"  სამიზნე: F(y) → {F_deep_target}")
    print(f"  შემოწმება: limit(F/F_target) = {deep_mond_limit}")

    g_N_expr, eq_full, chi_p_deep, chi_p_newton = spherical_field_equation()
    print("\n3. სფერული ფონური განტოლება")
    print(f"  g_N = {g_N_expr}")
    print(f"  სრული: {eq_full}")
    print(f"  Deep MOND χ' = {chi_p_deep}")
    print(f"  ნიუტონური χ' = {chi_p_newton}")

    g_eff, g_chi, M_b_BTFR = btfr_derivation()
    print("\n4. BTFR (Baryonic Tully-Fisher) ალგებრული ჩასმა")
    print(f"  g_χ (deep MOND) = {g_chi}")
    print(f"  g_eff = g_N + g_χ = {g_eff}")
    print(f"  M_b (deep MOND) = {M_b_BTFR}  (v^4 ∝ M_b)")

    print("\n5. a_0 კოსმოლოგიური მასშტაბი — რიცხობრივი დამთხვევა")
    results, a0_target = a0_cosmological_scale()
    print(f"  სამიზნე (Milgrom): a_0 = {a0_target:.4e} m/s²")
    for label, a0_pred, err in results:
        print(f"  {label:18s}: a_0 = {a0_pred:.4e} m/s² (ცდომილება {err:.2f}%)")

    print("\n6. Strategy 3 audit — μ(x) derivation status")
    mu_audit = mu_derivation_audit()
    for key, value in mu_audit.items():
        print(f"  {key:30s}: {value}")

    print("\n7. Strategy 3 audit — a_0 / Λ_eff bridge status")
    a0_audit = a0_lambda_eff_audit()
    for row in a0_audit["rows"]:
        print(
            f"  {row['scale']:24s}: a0={row['a0_m_s2']:.4e} m/s², "
            f"error={row['percent_error_vs_Milgrom']:.2f}%"
        )
    print(f"  verdict: {a0_audit['verdict']}")
    print(f"  status : {a0_audit['status']}")

    print("\n8. ღია ცდები (გადადის შემდეგ ფაზებში)")
    for i, task in enumerate(open_tasks(), 1):
        print(f"  {i}. {task}")

    print("\n9. სტატუსი")
    print("  - AQUAL ფარგლი დაფიქსირებულია (Bekenstein-Milgrom 1984)")
    print("  - μ(x) ფენომენოლოგიური არჩევანი, არა RFG-დან გამოყვანა")
    print("  - a_0 cH_0/(2π) რიცხობრივი დამთხვევაა (~10% ცდომილებით)")
    print("  - BTFR ალგებრული ჩასმაა, ემპირიულ ფიტს ეყრდნობა")
    print("  - SPARC შედარება, χ ლაგრანჟიანის RFG-დან გამოყვანა — ღია")
