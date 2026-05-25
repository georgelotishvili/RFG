# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 25: სრული EFT-of-Dark-Energy — no-ghost window + Bellini-Sawicki α-ები
================================================================================

სტატუსი:
Strategy 3 / X4+M3-ის შესრულება.

ამ ფაილის მიზანი:
1. phase21/phase23-ის ghost კონფლიქტის მკაცრი დახურვა:
       α_K = (-4*c_Y*X + 48*c_Y2*X^2)/(H^2*M_Pl^2)
   და α_K > 0 მოთხოვნის ცალკე ჩვენება X-სქემაში.

2. Y-სქემაში სრული background-dependent no-ghost პირობა:
       K_Phi = q00 * (c_Y + 6*c_Y2*Y0 + c_YI1*I1_bg) > 0

   FLRW normalized background:
       K_Phi = c_Y + 6*c_Y2 + 3*c_YI1/a^2 > 0
       a=1 -> c_Y + 6*c_Y2 + 3*c_YI1 > 0

3. Schwarzschild local/static ფონზე იგივე პირობის smoke-test.
   ადგილობრივი ორთონორმალური ჩარჩო უბრუნდება Minkowski/FLRW local პირობას;
   coordinate Phi=t smoke-test ცალკე იბეჭდება, რადგან კოორდინატული ნორმალიზაცია
   ფიზიკური no-ghost პირობა არ არის.

4. Solid sector აღარ იკარგება α-ებში: α_K, α_B, α_M, α_T-ში ჩნდება
   ESS/Ballesteros-Bellazzini ტიპის სიმბოლური დამატებები. სრული CMB fit მაინც
   phase32/hi_class ამოცანად რჩება.
"""

import sympy as sp


def horndeski_y_sector_alphas():
    """
    Bellini-Sawicki X scheme, Y = -2X.

    Pure Y-sector:
        G2(X) = -2*c_Y*X + 4*c_Y2*X^2
        alpha_K = (2X*G2_X + 4X^2*G2_XX)/(H^2*M_Pl^2)
                = (-4*c_Y*X + 48*c_Y2*X^2)/(H^2*M_Pl^2)

    If I1 is treated as a fixed background spurion, c_YI1*Y*I1 contributes
    by c_Y -> c_Y + c_YI1*I1_bg in G2.
    """
    X, I1_bg = sp.symbols("X I1_bg", real=True)
    H, M_Pl = sp.symbols("H M_Pl", positive=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    y_in_x = -2 * X
    G2_pure = c_Y * y_in_x + c_Y2 * y_in_x**2
    G2_with_i1 = (c_Y + c_YI1 * I1_bg) * y_in_x + c_Y2 * y_in_x**2

    def alpha_k(G2):
        G2_X = sp.diff(G2, X)
        G2_XX = sp.diff(G2, X, 2)
        return sp.simplify((2 * X * G2_X + 4 * X**2 * G2_XX) / (H**2 * M_Pl**2))

    return {
        "G2_pure": sp.expand(G2_pure),
        "alpha_T_Y_sector": sp.Integer(0),
        "alpha_M_Y_sector": sp.Integer(0),
        "alpha_B_Y_sector": sp.Integer(0),
        "alpha_K_pure": alpha_k(G2_pure),
        "alpha_K_with_I1_spurion": alpha_k(G2_with_i1),
        "ghost_rule_X_scheme": "require alpha_K_total > 0; for pure low-X branch this pushes c_Y < 0 in X convention",
    }


def y_scheme_no_ghost_conditions():
    """
    General Y-scheme quadratic coefficient for phase perturbation.

    Let Phi = Phi_bg + pi and Y = q00*(1 + pi_dot)^2 on a static time-like
    background. The pi_dot^2 coefficient is:

        K_Phi = q00 * (c_Y + 6*c_Y2*Y0 + c_YI1*I1_bg)

    q00 > 0 outside horizons, so the bracket controls the sign.
    """
    a, r, r_s, theta = sp.symbols("a r r_s theta", positive=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)
    f = 1 - r_s / r

    K_general = sp.Symbol("q00", positive=True) * (
        c_Y + 6 * c_Y2 * sp.Symbol("Y0", positive=True) + c_YI1 * sp.Symbol("I1_bg", positive=True)
    )

    K_flrw = sp.simplify(c_Y + 6 * c_Y2 + 3 * c_YI1 / a**2)
    K_flrw_today = sp.simplify(K_flrw.subs(a, 1))
    flrw_cYI1_bound = sp.solve_univariate_inequality(K_flrw_today > 0, c_YI1)
    flrw_cY2_bound = sp.solve_univariate_inequality(K_flrw_today > 0, c_Y2)

    # Schwarzschild, local orthonormal static frame: Y0=1, I1=3.
    K_schw_local = K_flrw_today

    # Coordinate smoke-test for Phi=t and solid labels (r, theta, phi), outside r>r_s.
    I1_schw_coord = sp.simplify(f + 1 / r**2 + 1 / (r**2 * sp.sin(theta) ** 2))
    K_schw_coord_bracket = sp.simplify(c_Y + 6 * c_Y2 / f + c_YI1 * I1_schw_coord)
    K_schw_coord = sp.simplify(K_schw_coord_bracket / f)
    K_schw_coord_equator = sp.simplify(K_schw_coord.subs(theta, sp.pi / 2))

    return {
        "K_general": K_general,
        "K_FLRW": K_flrw,
        "K_FLRW_today": K_flrw_today,
        "FLRW_today_c_YI1_bound": flrw_cYI1_bound,
        "FLRW_today_c_Y2_bound": flrw_cY2_bound,
        "K_Schwarzschild_local": K_schw_local,
        "I1_Schwarzschild_coordinate": I1_schw_coord,
        "K_Schwarzschild_coordinate_equator": K_schw_coord_equator,
    }


def sign_window_sweep():
    """
    Numeric smoke-test for c_Y2 > 0 and c_YI1 window.

    These are not fitted values; they are small examples that show the inequality
    catches pass/fail branches.
    """
    a_value = 1.0
    r_value = 10.0
    r_s_value = 2.0
    theta_value = sp.pi / 2

    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)
    a, r, r_s, theta = sp.symbols("a r r_s theta", positive=True)
    f = 1 - r_s / r

    K_flrw = c_Y + 6 * c_Y2 + 3 * c_YI1 / a**2
    K_schw_local = c_Y + 6 * c_Y2 + 3 * c_YI1
    I1_schw = f + 1 / r**2 + 1 / (r**2 * sp.sin(theta) ** 2)
    K_schw = (c_Y + 6 * c_Y2 / f + c_YI1 * I1_schw) / f

    cases = [
        {"name": "healthy_example", "c_Y": 1.0, "c_Y2": 0.10, "c_YI1": 0.0},
        {"name": "mixed_term_too_negative", "c_Y": 1.0, "c_Y2": 0.10, "c_YI1": -1.0},
        {"name": "c_Y2_negative_risky", "c_Y": 1.0, "c_Y2": -0.30, "c_YI1": 0.0},
    ]

    rows = []
    for case in cases:
        subs = {
            c_Y: case["c_Y"],
            c_Y2: case["c_Y2"],
            c_YI1: case["c_YI1"],
            a: a_value,
            r: r_value,
            r_s: r_s_value,
            theta: theta_value,
        }
        flrw_value = float(K_flrw.subs(subs))
        schw_local_value = float(K_schw_local.subs(subs))
        schw_value = float(K_schw.subs(subs))
        rows.append(
            {
                **case,
                "K_FLRW_today": flrw_value,
                "FLRW_status": "PASS" if flrw_value > 0 else "FAIL",
                "K_Schwarzschild_local": schw_local_value,
                "Schwarzschild_local_status": "PASS" if schw_local_value > 0 else "FAIL",
                "K_Schwarzschild_coord": schw_value,
                "Schwarzschild_coord_status": "PASS" if schw_value > 0 else "FAIL",
            }
        )
    return rows


def ess_solid_alpha_bookkeeping():
    """
    ESS/Ballesteros-Bellazzini style solid-sector bookkeeping.

    The file does not claim a full hi_class implementation. It makes the missing
    solid-sector terms explicit, so α_B/α_M/α_K are no longer silently set to zero.
    """
    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    H = sp.Function("H")(t)
    M_Pl = sp.Symbol("M_Pl", positive=True)
    X = sp.Symbol("X", real=True)
    I1_bg = sp.Symbol("I1_bg", positive=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    delta_M2 = sp.Function("delta_M2_solid")(a)
    alpha_B_solid = sp.Function("alpha_B_solid")(a)
    alpha_K_solid = sp.Function("alpha_K_solid")(a)
    alpha_T_solid = sp.Function("alpha_T_solid")(a)

    M_eff_sq = M_Pl**2 + delta_M2
    alpha_M_total = sp.simplify(sp.diff(sp.log(M_eff_sq), t) / H)
    alpha_K_y_i1 = sp.simplify(
        (-4 * X * (c_Y + c_YI1 * I1_bg) + 48 * c_Y2 * X**2) / (H**2 * M_eff_sq)
    )

    return {
        "M_eff_sq": M_eff_sq,
        "alpha_K_total": alpha_K_y_i1 + alpha_K_solid,
        "alpha_B_total": alpha_B_solid,
        "alpha_M_total": alpha_M_total,
        "alpha_T_total": alpha_T_solid,
        "GW170817_filter": "require alpha_T_solid ≈ 0 and phase34 tensor-speed filter",
        "solid_note": "delta alpha terms require ESS perturbation derivation before CLASS/hi_class fit",
    }


def observational_filters():
    return {
        "alpha_T": "|alpha_T| < O(1e-15) from GW170817/GRB170817A",
        "alpha_K": "alpha_K_total > 0 for scalar no-ghost",
        "alpha_M_alpha_B": "must be fit to CMB/LSS/BAO; not fixed by this smoke-test",
        "DESI_link": "static Lambda_eff is not enough for w(z); dynamic alpha-sector is needed",
    }


def class_camb_interface_open():
    return [
        "export alpha_K(a), alpha_B(a), alpha_M(a), alpha_T(a) arrays",
        "choose ESS closure for delta_M2_solid(a), alpha_B_solid(a), alpha_K_solid(a)",
        "run hi_class/CLASS Planck 2018 likelihood",
        "add BAO/LSS/DESI likelihoods after background H(a) is fixed",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 25: EFT-of-Dark-Energy — no-ghost + alpha sweep")
    print("=" * 72)

    print("\n1. Horndeski/Y-sector Bellini-Sawicki alphas")
    alphas = horndeski_y_sector_alphas()
    for key, value in alphas.items():
        print(f"  {key:30s}: {value}")

    print("\n2. Full Y-scheme no-ghost conditions")
    conditions = y_scheme_no_ghost_conditions()
    for key, value in conditions.items():
        print(f"  {key:34s}: {value}")

    print("\n3. c_Y2 / c_YI1 sign-window smoke-test")
    for row in sign_window_sweep():
        print(
            f"  {row['name']:24s}: "
            f"c_Y={row['c_Y']:+.2f}, c_Y2={row['c_Y2']:+.2f}, c_YI1={row['c_YI1']:+.2f} | "
            f"FLRW K={row['K_FLRW_today']:+.3f} {row['FLRW_status']} | "
            f"Schw local K={row['K_Schwarzschild_local']:+.3f} {row['Schwarzschild_local_status']} | "
            f"coord smoke K={row['K_Schwarzschild_coord']:+.3f} {row['Schwarzschild_coord_status']}"
        )
    print("  note: Schwarzschild coord smoke is not the physical ghost verdict; local K controls the sign.")

    print("\n4. ESS solid-sector alpha bookkeeping")
    ess = ess_solid_alpha_bookkeeping()
    for key, value in ess.items():
        print(f"  {key:22s}: {value}")

    print("\n5. Observational filters")
    for key, value in observational_filters().items():
        print(f"  {key:18s}: {value}")

    print("\n6. CLASS/hi_class open interface")
    for i, task in enumerate(class_camb_interface_open(), 1):
        print(f"  {i}. {task}")

    print("\n7. Status")
    print("  - Strategy 3 X4: background-dependent no-ghost window is now explicit.")
    print("  - Strategy 3 M3: alpha_K formula is corrected and solid-sector deltas are visible.")
    print("  - Full ESS perturbation derivation and Planck chi^2 fit remain phase32/hi_class work.")
