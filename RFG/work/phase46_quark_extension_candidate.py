# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 46: Quark-sector extension candidate from the same RFG holonomy

Status:
    Falsifiable extension sketch, not a fit.

Motivation:
    A particle-sector theory stronger than Koide-only lepton models should
    not stop at charged leptons. But quark masses are not pole-observable in
    the same clean way as charged leptons; they are scheme- and scale-dependent.
    Therefore this file defines a strict candidate rule and a falsification
    protocol rather than claiming a quark-mass derivation.

RFG rule:
    There is only one primitive reduced holonomy:

        theta_L = 2/9.

    Quark family phases, if present, must be projections of theta_L by
    topological/color/charge data. They are not new fitted angles.

Two minimal projection candidates:
    charge projection:
        theta_U = |Q_u| theta_L = 4/27,
        theta_D = |Q_d| theta_L = 2/27.

    inverse-color/empirical-Z3 projection:
        theta_U = theta_L / 3 = 2/27,
        theta_D = 2 theta_L / 3 = 4/27.

The second ordering matches a known Z3-Koide phenomenological suggestion,
but RFG should decide between them from color/framing topology, not from a
quark-mass fit.
"""

import math


THETA_LEPTON = 2.0 / 9.0


def c3_ratios(theta):
    values = sorted(
        1.0 + math.sqrt(2.0) * math.cos(theta + 2.0 * math.pi * k / 3.0)
        for k in range(3)
    )
    if values[0] <= 0:
        return None
    return [value / values[0] for value in values]


def mass_ratios_from_theta(theta):
    ratios = c3_ratios(theta)
    if ratios is None:
        return None
    return [value * value for value in ratios]


def projection_candidates():
    return [
        {
            "name": "charge_projection",
            "theta_U": 2.0 * THETA_LEPTON / 3.0,
            "theta_D": THETA_LEPTON / 3.0,
            "rule": "theta_f = |Q_f| theta_L",
            "free_angles": 0,
        },
        {
            "name": "inverse_color_empirical_z3_projection",
            "theta_U": THETA_LEPTON / 3.0,
            "theta_D": 2.0 * THETA_LEPTON / 3.0,
            "rule": "theta_U=theta_L/3, theta_D=2 theta_L/3",
            "free_angles": 0,
        },
    ]


def candidate_tables():
    rows = []
    for candidate in projection_candidates():
        theta_u = candidate["theta_U"]
        theta_d = candidate["theta_D"]
        rows.append(
            {
                **candidate,
                "freq_ratios_U": c3_ratios(theta_u),
                "mass_ratios_U": mass_ratios_from_theta(theta_u),
                "freq_ratios_D": c3_ratios(theta_d),
                "mass_ratios_D": mass_ratios_from_theta(theta_d),
            }
        )
    return rows


def falsification_protocol():
    return [
        "Choose a precise quark-mass scheme and scale, e.g. MSbar at a specified mu.",
        "Run all six quark masses to that scale with standard QCD RG equations.",
        "Normalize each sector by its lightest member: (u,c,t) and (d,s,b).",
        "Compare against the two zero-angle-fit RFG projection candidates.",
        "Reject any candidate whose required theta_U/theta_D differs from its topological projection beyond uncertainties.",
        "Do not tune theta_U or theta_D after seeing the masses; that would erase the prediction.",
    ]


def status_assessment():
    return {
        "strength": "extends RFG particle sector beyond charged leptons without new continuous angles",
        "warning": "quark masses are scheme/scale dependent, so no pole-mass claim is made here",
        "needed_theorem": "derive whether color/framing topology chooses charge_projection or inverse_color projection",
        "relation_to_leptons": "theta_L=2/9 remains the only primitive phase",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 46: Quark extension candidate")
    print("=" * 72)

    print("\n1. Projection candidates")
    for row in candidate_tables():
        print(f"  {row['name']}")
        print(f"    rule     : {row['rule']}")
        print(f"    theta_U  : {row['theta_U']:.12f}")
        print(f"    theta_D  : {row['theta_D']:.12f}")
        print(f"    U freq   : {[round(x, 6) for x in row['freq_ratios_U']]}")
        print(f"    U mass   : {[round(x, 6) for x in row['mass_ratios_U']]}")
        print(f"    D freq   : {[round(x, 6) for x in row['freq_ratios_D']]}")
        print(f"    D mass   : {[round(x, 6) for x in row['mass_ratios_D']]}")

    print("\n2. Falsification protocol")
    for item in falsification_protocol():
        print(f"  - {item}")

    print("\n3. Status")
    for key, value in status_assessment().items():
        print(f"  {key:16s}: {value}")
