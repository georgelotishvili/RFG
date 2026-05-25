# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 24: Principal symbol, hyperbolicity, Cauchy well-posedness
================================================================================

სტატუსი:
Strategy 3 / M2-ის შესრულება. ეს ფაილი აღარ არის სქელეტი: იგი აგებს
quadratic principal Lagrangian-ს perturbation-ებისთვის და ამოწმებს
characteristic determinant-ს.

მთავარი ობიექტი:
    L2_principal = A*dPhi_dot^2 + B*pi_L_dot^2
                   + C*dPhi_z^2 + D*pi_L_z^2
                   + M*(dPhi_dot*pi_L_z + pi_L_dot*dPhi_z)

Plane wave-ით exp(i(kz - omega t)):
    det M(s) = (A*s + C)*(B*s + D) - M^2*s = 0,
    s = omega^2/k^2.

თუ mixing M=0:
    s_phi = -C/A,   s_L = -D/B.

FLRW-ზე s არის comoving speed squared; physical speed is a^2*s.
Schwarzschild-ზე principal symbol locally orthonormal frame-ში იგივეა, ხოლო
coordinate radial relation არის omega_coord^2 = f(r)^2*c_local^2*k_r^2.
"""

import sympy as sp

from phase1_action import analyze_sound_speeds


def quadratic_principal_coefficients(scale_factor=None):
    """
    Build principal coefficients around FLRW/Minkowski background.

    scale_factor=None means symbolic a. scale_factor=1 gives Minkowski.
    Perturbations propagate along comoving z. The coefficients are the closed
    FLRW scaling of the second-order expansion; at a=1 they are verified
    against phase1_action.analyze_sound_speeds().
    """
    a = sp.Symbol("a", positive=True) if scale_factor is None else sp.sympify(scale_factor)

    dPhi_dot, dPhi_z = sp.symbols("dPhi_dot dPhi_z", real=True)
    pi1_dot, pi3_dot = sp.symbols("pi1_dot pi3_dot", real=True)
    pi1_z, pi3_z = sp.symbols("pi1_z pi3_z", real=True)
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )

    A = sp.simplify(c_Y + 6 * c_Y2 + 3 * c_YI1 / a**2)
    B_long = sp.simplify(-c_I1 - c_YI1 - (6 * c_I1sq + 2 * c_I2) / a**2 - c_I3 / a**4)
    C = sp.simplify(-(c_Y + 2 * c_Y2) / a**2 - 3 * c_YI1 / a**4)
    D = sp.simplify((c_I1 + c_YI1) / a**2 + (10 * c_I1sq + 2 * c_I2) / a**4 + c_I3 / a**6)
    M_mix = sp.simplify(4 * c_YI1 / a**2)

    K_T = B_long
    C_T = sp.simplify((c_I1 + c_YI1) / a**2 + (6 * c_I1sq + c_I2) / a**4)
    cs2_T_comoving = sp.simplify(-C_T / K_T)
    L2 = sp.simplify(
        A * dPhi_dot**2
        + B_long * pi3_dot**2
        + C * dPhi_z**2
        + D * pi3_z**2
        + M_mix * (dPhi_dot * pi3_z + pi3_dot * dPhi_z)
        + K_T * pi1_dot**2
        + C_T * pi1_z**2
    )

    return {
        "a": a,
        "L2": L2,
        "A": A,
        "B_long": B_long,
        "C": C,
        "D": D,
        "M_mix": M_mix,
        "K_T": K_T,
        "C_T": C_T,
        "cs2_T_comoving": cs2_T_comoving,
        "cs2_T_physical": sp.simplify(a**2 * cs2_T_comoving),
    }


def characteristic_polynomial(coeffs, solve_symbolic=False):
    """det M(s)=0 for mixed phase-longitudinal sector."""
    s = sp.Symbol("s", real=True)
    A = coeffs["A"]
    B = coeffs["B_long"]
    C = coeffs["C"]
    D = coeffs["D"]
    M = coeffs["M_mix"]
    det = sp.factor(sp.simplify((A * s + C) * (B * s + D) - M**2 * s))
    roots = sp.solve(det, s) if solve_symbolic else None
    return s, det, roots


def minkowski_principal_symbol():
    coeffs = quadratic_principal_coefficients(scale_factor=1)
    s, det, roots = characteristic_polynomial(coeffs)
    return coeffs, s, det, roots


def flrw_principal_symbol():
    coeffs = quadratic_principal_coefficients(scale_factor=None)
    s, det, roots = characteristic_polynomial(coeffs)
    physical_det = sp.factor(sp.simplify(det.subs(s, sp.Symbol("c_phys2", real=True) / coeffs["a"]**2)))
    return coeffs, s, det, roots, physical_det


def decoupled_mixed_roots(coeffs):
    """Closed roots when c_YI1=0 removes phase-longitudinal mixing."""
    c_YI1 = sp.Symbol("c_YI1", real=True)
    A = coeffs["A"].subs(c_YI1, 0)
    B = coeffs["B_long"].subs(c_YI1, 0)
    C = coeffs["C"].subs(c_YI1, 0)
    D = coeffs["D"].subs(c_YI1, 0)
    return sp.simplify(-C / A), sp.simplify(-D / B)


def schwarzschild_local_symbol():
    """
    Local orthonormal principal symbol outside horizon.

    The local Cauchy problem is governed by the same determinant as Minkowski.
    Coordinate radial propagation is redshifted:
        omega_coord^2 = f(r)^2 * c_local^2 * k_r^2.
    """
    r, r_s, c_local2, k_r = sp.symbols("r r_s c_local2 k_r", positive=True)
    omega = sp.Symbol("omega", real=True)
    f = 1 - r_s / r
    radial_dispersion = sp.Eq(omega**2, sp.simplify(f**2 * c_local2 * k_r**2))
    return {
        "f": f,
        "local_det": minkowski_principal_symbol()[2],
        "radial_coordinate_dispersion": radial_dispersion,
        "horizon_note": "coordinate speed -> 0 as f -> 0, local characteristic remains finite",
    }


def verify_against_phase1():
    """
    Compare phase24 coefficients with phase1_action.analyze_sound_speeds().
    We compare coefficients, not phase1's printed characteristic convention.
    """
    phase24 = minkowski_principal_symbol()[0]
    cs2_T_phase1, _eq_phase1, coeffs_phase1, _roots_phase1 = analyze_sound_speeds()
    checks = {
        "A_difference": sp.simplify(phase24["A"] - coeffs_phase1["A"]),
        "B_long_difference": sp.simplify(phase24["B_long"] - coeffs_phase1["B_pi3"]),
        "C_difference": sp.simplify(phase24["C"] - coeffs_phase1["C"]),
        "D_difference": sp.simplify(phase24["D"] - coeffs_phase1["D"]),
        "M_mix_difference": sp.simplify(phase24["M_mix"] - coeffs_phase1["M_mix"]),
        "transverse_speed_difference": sp.simplify(phase24["cs2_T_comoving"] - cs2_T_phase1),
    }
    return checks


def numeric_hyperbolicity_cases():
    """
    Numeric smoke-test. PASS means:
        kinetic coefficients are positive,
        mixed roots are real and positive,
        transverse speed is positive.
    """
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )
    coeffs, s, det, _ = minkowski_principal_symbol()

    cases = [
        {
            "name": "stable_decoupled",
            c_Y: 1.0,
            c_Y2: 0.10,
            c_I1: -2.0,
            c_I1sq: 0.10,
            c_I2: -0.10,
            c_I3: 0.0,
            c_YI1: 0.0,
        },
        {
            "name": "stable_mixed_small",
            c_Y: 1.0,
            c_Y2: 0.10,
            c_I1: -2.0,
            c_I1sq: 0.10,
            c_I2: -0.10,
            c_I3: 0.0,
            c_YI1: 0.05,
        },
        {
            "name": "ghost_fail_phase",
            c_Y: -1.0,
            c_Y2: 0.05,
            c_I1: -2.0,
            c_I1sq: 0.10,
            c_I2: -0.10,
            c_I3: 0.0,
            c_YI1: 0.0,
        },
        {
            "name": "gradient_fail_solid",
            c_Y: 1.0,
            c_Y2: 0.10,
            c_I1: -0.20,
            c_I1sq: 0.00,
            c_I2: 0.50,
            c_I3: 0.0,
            c_YI1: 0.0,
        },
    ]

    rows = []
    for case in cases:
        name = case["name"]
        subs = {key: value for key, value in case.items() if key != "name"}
        kinetic_values = [
            float(sp.N(coeffs["A"].subs(subs))),
            float(sp.N(coeffs["B_long"].subs(subs))),
            float(sp.N(coeffs["K_T"].subs(subs))),
        ]
        roots = [complex(root) for root in sp.nroots(sp.Poly(det.subs(subs), s), n=20, maxsteps=100)]
        roots_real = [root.real for root in roots if abs(root.imag) < 1e-7]
        transverse = float(sp.N(coeffs["cs2_T_comoving"].subs(subs)))

        kinetic_pass = all(value > 0 for value in kinetic_values)
        roots_real_pass = len(roots_real) == 2
        roots_positive = roots_real_pass and all(value > 0 for value in roots_real)
        transverse_pass = transverse > 0
        overall = kinetic_pass and roots_positive and transverse_pass

        rows.append(
            {
                "name": name,
                "kinetic_values": kinetic_values,
                "mixed_roots": roots_real,
                "transverse_cs2": transverse,
                "kinetic": "PASS" if kinetic_pass else "FAIL",
                "mixed_roots_status": "PASS" if roots_positive else "FAIL",
                "transverse_status": "PASS" if transverse_pass else "FAIL",
                "overall": "PASS" if overall else "FAIL",
            }
        )
    return rows


def status_assessment():
    return {
        "minkowski": "det M(s)=0 computed from L2 principal coefficients",
        "flrw": "comoving det M(s)=0 computed; physical speed is a^2*s",
        "schwarzschild": "local orthonormal determinant equals Minkowski; coordinate radial redshift added",
        "remaining": "global hyperbolicity and full curved-background perturbation system remain beyond smoke-test",
    }


def compact_det_label():
    return "(A*s + C)*(B*s + D) - M_mix**2*s = 0"


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 24: Principal symbol + hyperbolicity")
    print("=" * 72)

    print("\n1. Minkowski principal symbol")
    coeffs_m, s, det_m, roots_m = minkowski_principal_symbol()
    for key in ["A", "B_long", "C", "D", "M_mix", "K_T", "C_T"]:
        print(f"  {key:10s}: {coeffs_m[key]}")
    print(f"  det M(s): {compact_det_label()}")
    print(f"  expanded degree in s: {sp.degree(det_m, s)}")
    print(f"  decoupled roots (c_YI1=0): {decoupled_mixed_roots(coeffs_m)}")
    print(f"  transverse c_s^2: {coeffs_m['cs2_T_comoving']}")

    print("\n2. FLRW principal symbol")
    coeffs_f, s_f, det_f, roots_f, det_phys = flrw_principal_symbol()
    print(f"  a: {coeffs_f['a']}")
    for key in ["A", "B_long", "C", "D", "M_mix"]:
        print(f"  {key:10s}: {coeffs_f[key]}")
    print(f"  det_comoving M(s): {compact_det_label()}")
    print(f"  expanded degree in s: {sp.degree(det_f, s_f)}")
    print("  physical substitution: s = c_phys^2/a^2")
    print(f"  transverse c_phys^2: {coeffs_f['cs2_T_physical']}")

    print("\n3. Schwarzschild local principal symbol")
    sch = schwarzschild_local_symbol()
    print(f"  f                             : {sch['f']}")
    print(f"  local_det                     : {compact_det_label()}")
    print(f"  radial_coordinate_dispersion  : {sch['radial_coordinate_dispersion']}")
    print(f"  horizon_note                  : {sch['horizon_note']}")

    print("\n4. phase1 coefficient verification")
    for key, value in verify_against_phase1().items():
        print(f"  {key:30s}: {value}")

    print("\n5. Numeric hyperbolicity smoke-test")
    for row in numeric_hyperbolicity_cases():
        roots_text = ", ".join(f"{value:.4g}" for value in row["mixed_roots"])
        kinetic_text = ", ".join(f"{value:.4g}" for value in row["kinetic_values"])
        print(f"  {row['name']:22s}: overall={row['overall']}")
        print(f"    K values: {kinetic_text} -> {row['kinetic']}")
        print(f"    mixed roots: {roots_text} -> {row['mixed_roots_status']}")
        print(f"    transverse c_s^2: {row['transverse_cs2']:.4g} -> {row['transverse_status']}")

    print("\n6. Status")
    for key, value in status_assessment().items():
        print(f"  {key:14s}: {value}")
