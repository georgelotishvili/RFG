# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 41: Action-level normal form for C3 phase locking

Status:
    Algebraic support for the phase37-40 chain.

Question:
    Is the C3 phase-locking term an extra hand-added ingredient, or does it
    already exist in the RFG supersolid action P(Y, I1, I2, I3)?

Answer:
    It already exists generically through the determinant invariant I3.

Setup:
    Around the isotropic elastic background, write the 3x3 material matrix as

        B = b I + Q,       tr(Q) = 0.

    The traceless part Q is the triaxial strain. In its principal-axis frame

        Q = diag(q1, q2, q3),        q1 + q2 + q3 = 0.

    Define the complex C3 strain coordinate

        E = q1 + omega q2 + omega^2 q3,     omega = exp(2*pi*i/3).

Key identities:
    1. det(b I + Q) contains det(Q).
    2. det(Q) is the cubic triaxial invariant.
    3. E^3 + conjugate(E)^3 = 27 det(Q).

Therefore any generic I3 dependence in the action supplies the
phase-sensitive C3 anisotropy:

        V_C3 ~ -lambda_3 Re(E^3).

This gives the normal-form origin of the C3 strain-sector locking used in
phase39. Combined with the oriented-triad C3 and the h=2 framed closure from
phase40, phase38 gives theta = 2/9.
"""

import sympy as sp


def invariants_around_isotropic_background():
    """
    Expand I1, I2, I3 for B = b I + Q, tr(Q)=0.
    """
    b, q1, q2 = sp.symbols("b q1 q2", real=True)
    q3 = -q1 - q2
    eigenvalues = [b + q1, b + q2, b + q3]

    I1 = sp.simplify(sum(eigenvalues))
    I2 = sp.simplify(
        eigenvalues[0] * eigenvalues[1]
        + eigenvalues[0] * eigenvalues[2]
        + eigenvalues[1] * eigenvalues[2]
    )
    I3 = sp.simplify(eigenvalues[0] * eigenvalues[1] * eigenvalues[2])

    tr_q2 = sp.simplify(q1**2 + q2**2 + q3**2)
    tr_q3 = sp.simplify(q1**3 + q2**3 + q3**3)
    det_q = sp.simplify(q1 * q2 * q3)

    return {
        "I1": I1,
        "I2": I2,
        "I3": I3,
        "tr_Q2": tr_q2,
        "tr_Q3": tr_q3,
        "det_Q": det_q,
        "I3_expected": sp.simplify(b**3 - b * tr_q2 / 2 + det_q),
        "trQ3_minus_3detQ": sp.simplify(tr_q3 - 3 * det_q),
    }


def c3_complex_strain_identity():
    """
    Prove E^3 + Ebar^3 = 27 det(Q) for q1+q2+q3=0.
    """
    q1, q2 = sp.symbols("q1 q2", real=True)
    q3 = -q1 - q2
    omega = -sp.Rational(1, 2) + sp.I * sp.sqrt(3) / 2
    omega_bar = sp.conjugate(omega)

    E = sp.simplify(q1 + omega * q2 + omega**2 * q3)
    E_bar = sp.simplify(q1 + omega_bar * q2 + omega_bar**2 * q3)
    det_q = sp.simplify(q1 * q2 * q3)
    tr_q3 = sp.simplify(q1**3 + q2**3 + q3**3)

    identity_det = sp.simplify(sp.expand(E**3 + E_bar**3 - 27 * det_q))
    identity_tr = sp.simplify(sp.expand(E**3 + E_bar**3 - 9 * tr_q3))

    return {
        "E": E,
        "E_bar": E_bar,
        "E3_plus_Ebar3_minus_27detQ": identity_det,
        "E3_plus_Ebar3_minus_9trQ3": identity_tr,
        "det_Q": det_q,
        "tr_Q3": tr_q3,
    }


def polar_triaxial_parameterization():
    """
    Parameterize the traceless eigenvalues by a radius rho and phase beta:

        q_i = rho cos(beta + 2*pi*i/3).

    Then det(Q) = rho^3 cos(3 beta) / 4.
    """
    rho, beta = sp.symbols("rho beta", real=True)
    q1 = rho * sp.cos(beta)
    q2 = rho * sp.cos(beta + 2 * sp.pi / 3)
    q3 = rho * sp.cos(beta + 4 * sp.pi / 3)
    det_q = sp.trigsimp(sp.expand_trig(q1 * q2 * q3))
    tr_q2 = sp.trigsimp(sp.expand_trig(q1**2 + q2**2 + q3**2))
    det_expected = rho**3 * sp.cos(3 * beta) / 4
    trq2_expected = 3 * rho**2 / 2
    return {
        "q1": q1,
        "q2": q2,
        "q3": q3,
        "sum_q": sp.trigsimp(sp.expand_trig(q1 + q2 + q3)),
        "det_Q": det_q,
        "det_Q_minus_expected": sp.trigsimp(sp.expand_trig(det_q - det_expected)),
        "tr_Q2": tr_q2,
        "tr_Q2_minus_expected": sp.trigsimp(sp.expand_trig(tr_q2 - trq2_expected)),
    }


def normal_form_from_action():
    """
    Minimal effective potential inherited from the action.

    The sign and size of lambda_3 depend on the microscopic RFG coefficients,
    but the existence of the cubic C3 anisotropy is generic once I3 is present.
    """
    rho, beta = sp.symbols("rho beta", real=True)
    a2, a4, lambda3 = sp.symbols("a2 a4 lambda3", real=True)
    det_q = rho**3 * sp.cos(3 * beta) / 4
    potential = a2 * rho**2 + a4 * rho**4 - lambda3 * det_q
    d_beta = sp.diff(potential, beta)
    stationary_condition = sp.factor(d_beta)
    curvature_beta = sp.factor(sp.diff(potential, beta, 2))
    return {
        "V_eff": potential,
        "dV_d_beta": stationary_condition,
        "d2V_d_beta2": curvature_beta,
        "stationary_rule": "sin(3 beta)=0 -> three C3-locked strain sectors",
    }


def theta_chain_summary():
    return {
        "action_invariant": "I3 = det(B)",
        "triaxial_cubic": "det(Q) = Re(E^3)/27 = tr(Q^3)/3",
        "strain_lock": "V_C3 ~ -lambda3 rho^3 cos(3 beta)/4",
        "phase39": "oriented triad C3 x strain-sector C3 -> 9 closure slots",
        "phase40": "first non-trivial oriented framed branch h=2",
        "phase38": "theta = h/9 = 2/9",
        "phase37": "C3 Koide operator gives charged-lepton mass ratios",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 41: Action normal form for C3 phase locking")
    print("=" * 72)

    print("\n1. Invariants around B = b I + Q")
    inv = invariants_around_isotropic_background()
    for key, value in inv.items():
        print(f"  {key:28s}: {value}")

    print("\n2. Complex C3 strain identity")
    ident = c3_complex_strain_identity()
    for key, value in ident.items():
        print(f"  {key:34s}: {value}")

    print("\n3. Polar triaxial parameterization")
    polar = polar_triaxial_parameterization()
    for key, value in polar.items():
        print(f"  {key:28s}: {value}")

    print("\n4. Normal form inherited from I3")
    nf = normal_form_from_action()
    for key, value in nf.items():
        print(f"  {key:20s}: {value}")

    print("\n5. Theta chain")
    for key, value in theta_chain_summary().items():
        print(f"  {key:18s}: {value}")
