# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 9: gravitational waves - old ISPG predictions in RFG language.

Recovered old-theory predictions:
1. Tensor propagation is exactly luminal: alpha_T=0 -> c_g=c.
2. Two GR tensor modes remain.
3. A scalar/breathing channel is allowed by the scalar medium, but its
   source amplitude is parametrically suppressed.
4. Leading dipole radiation cancels when the compact-body sensitivity is
   universal, s=1/2; residuals require a full PN source calculation.

This file keeps the distinction sharp:
    c_g=c is closed structurally;
    breathing amplitude is a controlled working estimate;
    scalar-dipole nulling is leading-order, not a final waveform theorem.
"""

import sympy as sp


def analyze_gw_full():
    """
    Solid-sector TT check already used in earlier RFG phases.

    The RFG solid invariants do not generate h_dot^2 or h_z^2 corrections
    for a TT perturbation on FLRW.  Therefore the solid sector does not
    shift c_T.  A mass term constraint remains separate.
    """
    t, z = sp.symbols('t z', real=True)
    h = sp.Function('h')(t, z)
    h_dot = sp.diff(h, t)
    h_z = sp.diff(h, z)
    eps = sp.Symbol('eps', real=True)
    a = sp.Symbol('a', real=True, positive=True)

    g_11 = -a**2 * (1 - eps * h)
    g_22 = -a**2 * (1 + eps * h)
    g_33 = -a**2

    det_g = a**6 * (1 - eps**2 * h**2)
    sqrt_g = sp.sqrt(det_g)

    g_inv_11 = sp.series(1 / g_11, eps, 0, 3).removeO()
    g_inv_22 = sp.series(1 / g_22, eps, 0, 3).removeO()
    g_inv_33 = 1 / g_33

    b11 = -g_inv_11
    b22 = -g_inv_22
    b33 = -g_inv_33

    i1_pert = sp.simplify(b11 + b22 + b33)
    i2_pert = sp.simplify(sp.Rational(1, 2) * (i1_pert**2 - (b11**2 + b22**2 + b33**2)))
    i3_pert = sp.simplify(b11 * b22 * b33)
    y_pert = 1

    from phase1_action import get_polynomial_lagrangian

    y_s, i1_s, i2_s, i3_s = sp.symbols('Y I1 I2 I3', real=True)
    l_poly = get_polynomial_lagrangian(y_s, i1_s, i2_s, i3_s)

    l_eval = l_poly.subs({y_s: y_pert, i1_s: i1_pert, i2_s: i2_pert, i3_s: i3_pert})
    l_density = sp.series(sqrt_g * l_eval, eps, 0, 3).removeO()
    l_o2 = sp.simplify(l_density.coeff(eps, 2) / a**3)

    coeff_h_dot2 = l_o2.coeff(h_dot**2)
    coeff_h_z2 = l_o2.coeff(h_z**2)
    mass_term = sp.simplify(l_o2.coeff(h**2))

    return coeff_h_dot2, coeff_h_z2, mass_term


def analyze_horndeski_luminal_speed():
    """
    Old Appendix 10 result in RFG notation.

    Horndeski tensor-speed excess alpha_T receives contributions from
    G_{4,X} and G_5.  RFG's Einstein-Hilbert backbone has:
        G4 = const, G4_X = 0, G5 = 0.
    """
    c, c_g, alpha_T, G4_X, G5 = sp.symbols('c c_g alpha_T G4_X G5', real=True)
    alpha_t_value = sp.Integer(0)
    c_g_value = c * sp.sqrt(1 + alpha_t_value)

    return {
        "theorem": "tensor gravitational waves are exactly luminal",
        "Horndeski_conditions": "G4=const, G4_X=0, G5=0",
        "alpha_T_definition": sp.Eq(alpha_T, c_g**2 / c**2 - 1),
        "G4_X": sp.Eq(G4_X, 0),
        "G5": sp.Eq(G5, 0),
        "alpha_T": sp.Eq(alpha_T, alpha_t_value),
        "c_g": sp.Eq(c_g, c_g_value),
        "GW170817_status": "satisfied structurally, not by parameter tuning",
    }


def analyze_scalar_breathing_estimate():
    """
    Breathing-mode working estimate inherited from the old theory.

    If scalar charge per mass is universal at leading order, s=1/2, the
    monopole is stationary and the leading dipole cancels.  The first
    candidate radiative scalar channel is the trace quadrupole, estimated
    as A_b/A_t ~ v^2/c^2 ~ r_s/r.
    """
    r, r_s, v, c, s_A, s_B = sp.symbols(
        'r r_s v c s_A s_B',
        real=True,
        positive=True,
    )
    sensitivity_universal = sp.Eq(s_A, sp.Rational(1, 2))
    sensitivity_match = sp.Eq(s_A - s_B, 0)
    amplitude_ratio = sp.simplify(v**2 / c**2)
    virial_ratio = r_s / r

    return {
        "scalar_channel": "breathing polarization h_b is allowed",
        "universal_sensitivity": sensitivity_universal,
        "dipole_charge_difference": sensitivity_match,
        "dipole_status": "leading dipole cancels when s_A=s_B=1/2",
        "amplitude_ratio_working": sp.Eq(sp.Symbol('A_b/A_t'), amplitude_ratio),
        "virial_estimate": sp.Eq(sp.Symbol('A_b/A_t'), virial_ratio),
        "weak_field_example_r_10rs": sp.N(virial_ratio.subs(r, 10 * r_s), 8),
        "status": "parametric working estimate; full PN scalar quadrupole coefficient remains a waveform task",
    }


def gw_prediction_ledger():
    return [
        "Closed: alpha_T=0 -> c_g=c exactly.",
        "Closed: the solid sector adds no h_dot^2 or h_z^2 TT kinetic-gradient correction.",
        "Constraint: the TT mass term must be tuned/constrained to avoid massive graviton dispersion.",
        "Recovered old estimate: scalar breathing amplitude A_b/A_t ~ r_s/r.",
        "Leading scalar dipole is suppressed if compact-body sensitivity is universal, s=1/2.",
        "Open waveform task: exact scalar quadrupole coefficient and comparable-mass IMR templates.",
        "ISCO proxy from phase18: f_ISCO=0.931 f_ISCO_GR is a strong-field timing target, not a full waveform by itself.",
    ]


if __name__ == "__main__":
    coeff_h_dot2, coeff_h_z2, mass_term_flrw = analyze_gw_full()
    a = sp.Symbol('a', real=True, positive=True)
    mass_term_mink = sp.simplify(mass_term_flrw.subs(a, 1))

    print("--- solid-sector TT check: c_T არ იცვლება ---")
    print(f"L_solid-ის კინეტიკური წევრი (h_dot^2): {coeff_h_dot2}")
    print(f"L_solid-ის გრადიენტული წევრი (h_z^2): {coeff_h_z2}")
    print("დასკვნა: L_solid არ შეიცავს h_dot^2 ან h_z^2 წევრებს.")
    print("c_T^2 = c^2 (1 + delta), სადაც delta = 0 ზუსტად სრულდება.")
    print("\nთუმცა L_solid წარმოქმნის გრავიტონის ეფექტურ მასას (m_g^2 * h^2):")
    print(f"FLRW ფონზე მასის კოეფიციენტი: {mass_term_flrw}")
    print(f"Minkowski ფონზე (a=1): {mass_term_mink}")
    print("GW170817-ის და მასიური დისპერსიის ასარიდებლად მოითხოვება ეს mass-term constraint.")

    print("\n--- Horndeski/EFT luminal theorem ---")
    for key, value in analyze_horndeski_luminal_speed().items():
        print(f"{key:28s}: {value}")

    print("\n--- scalar breathing / dipole ledger ---")
    for key, value in analyze_scalar_breathing_estimate().items():
        print(f"{key:28s}: {value}")

    print("\n--- prediction ledger ---")
    for item in gw_prediction_ledger():
        print(f"  - {item}")
