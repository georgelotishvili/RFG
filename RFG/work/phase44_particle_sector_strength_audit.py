# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 44: Particle-sector strength and novelty audit

Purpose:
    Make the particle-sector claim sharper without overstating it.

Core point:
    Many models already use Koide, Z3, family symmetry, or the phase
    delta_L ~= 2/9. RFG is not novel merely because it reproduces Koide or
    notices 2/9.

RFG's stronger claim is the specific action/topology chain:

    I3 = det(B)
      -> det(Q) = Re(E^3)/27
      -> C3 strain lock from the supersolid action
      -> oriented elastic triad gives another C3
      -> C3 x C3 = 9 reduced closure slots
      -> oriented framed branch h=2
      -> theta = 2/9
      -> C3 Koide operator gives charged-lepton mass ratios.

This file also cleans up an internal tension:
    The older N=5,72,295 ladder and N=4 prediction are conditional legacy
    candidates. If phase37-43 C3 theorem-chain is adopted, the old radial
    N-ladder should no longer be presented as the main particle-sector route.
"""

from phase37_c3_koide_operator import prediction_table, koide_identity
from phase42_charged_lepton_theorem_audit import (
    audit_action_c3_lock,
    audit_mass_predictions,
    audit_z9_and_h2,
    open_theorem_items,
)


def comparison_with_existing_approaches():
    """
    High-level comparison. References are given by author/model labels rather
    than used as authority inside the code.
    """
    return [
        {
            "approach": "Koide original / scalar-potential family models",
            "has_koide": True,
            "has_delta_2_over_9": False,
            "action_origin_c3": "family scalar potential, not RFG medium",
            "topological_h2": False,
            "main_gap": "mass relation protected/engineered, but not a vacuum-medium oscillon derivation",
        },
        {
            "approach": "Z3-symmetric Koide parametrization",
            "has_koide": True,
            "has_delta_2_over_9": True,
            "action_origin_c3": "usually parametrization or family symmetry",
            "topological_h2": False,
            "main_gap": "2/9 is observed/suggested, not derived from an elastic invariant plus framed closure",
        },
        {
            "approach": "Sumino/family-gauge protection",
            "has_koide": True,
            "has_delta_2_over_9": "model-dependent",
            "action_origin_c3": "new family gauge sector",
            "topological_h2": False,
            "main_gap": "strong radiative idea, but different ontology and no supersolid oscillon route",
        },
        {
            "approach": "A4/S3/Z3 flavor models",
            "has_koide": "sometimes",
            "has_delta_2_over_9": "sometimes",
            "action_origin_c3": "discrete flavor symmetry input",
            "topological_h2": False,
            "main_gap": "symmetry is typically imposed in flavor space",
        },
        {
            "approach": "RFG phase37-43 chain",
            "has_koide": True,
            "has_delta_2_over_9": True,
            "action_origin_c3": "I3=det(B) -> Re(E^3) C3 lock",
            "topological_h2": True,
            "main_gap": "full 3D oscillon/PDE proof and radiative protection still open",
        },
    ]


def internal_consistency_cleanup():
    """
    What must be demoted or relabeled after adopting the C3 theorem-chain.
    """
    return [
        {
            "item": "phase35 radial n=14 route",
            "old_status": "candidate route for muon as a radial overtone",
            "new_status": "legacy audit / phenomenological comparison",
            "reason": "phase37 explains muon and tau together as one C3 triplet; no missing n=2..13 issue",
        },
        {
            "item": "N=5,72,295 ladder",
            "old_status": "empirical index ladder",
            "new_status": "do not present as primary derivation",
            "reason": "indices are back-solved from sqrt(m); C3 ratios are cleaner",
        },
        {
            "item": "phase36 N=4 329 keV prediction",
            "old_status": "conditional RFG-specific prediction if N-ladder is derived",
            "new_status": "conditional legacy prediction, suspended under C3 route",
            "reason": "C3 charged-lepton sector does not imply a radial N=4 companion",
        },
        {
            "item": "toy/work-2 scripts",
            "old_status": "exploratory sandbox",
            "new_status": "not part of main proof chain",
            "reason": "main particle-sector chain is now phase37-43",
        },
    ]


def rfg_strength_scorecard():
    action_checks = audit_action_c3_lock()
    z9_checks = audit_z9_and_h2()
    mass_checks = audit_mass_predictions()

    return [
        {
            "criterion": "Koide identity",
            "status": "closed algebraically",
            "evidence": f"K_C3 = {koide_identity():.12f}",
        },
        {
            "criterion": "charged-lepton ratios",
            "status": "closed as pole-mass-level postdiction",
            "evidence": all(mass_checks.values()),
        },
        {
            "criterion": "C3 action origin",
            "status": "closed in principal-axis normal form",
            "evidence": all(action_checks.values()),
        },
        {
            "criterion": "theta=2/9 route",
            "status": "candidate derivation",
            "evidence": all(z9_checks.values()),
        },
        {
            "criterion": "local h=2 stability",
            "status": "normal-form conditions derived",
            "evidence": "phase43 Hessian",
        },
        {
            "criterion": "full 3D particle proof",
            "status": "open",
            "evidence": "needs localized oscillon fluctuation operator",
        },
        {
            "criterion": "radiative stability",
            "status": "open",
            "evidence": "needs RFG analogue of Koide/Sumino protection or pole-frequency argument",
        },
        {
            "criterion": "absolute electron mass",
            "status": "open",
            "evidence": "current chain predicts ratios anchored to m_e",
        },
    ]


def strongest_defensible_claim():
    return (
        "RFG has a stronger-than-phenomenological charged-lepton candidate "
        "because its C3 structure is traced to the supersolid invariant "
        "I3=det(B) and to framed topological closure. It is not yet a final "
        "particle theory until full PDE stability, radiative protection, and "
        "absolute scale are derived."
    )


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 44: Particle-sector strength and novelty audit")
    print("=" * 72)

    print("\n1. Existing approaches vs RFG")
    for row in comparison_with_existing_approaches():
        print(f"  {row['approach']}")
        print(f"    Koide: {row['has_koide']}")
        print(f"    delta 2/9: {row['has_delta_2_over_9']}")
        print(f"    C3 origin: {row['action_origin_c3']}")
        print(f"    topological h=2: {row['topological_h2']}")
        print(f"    gap: {row['main_gap']}")

    print("\n2. Internal consistency cleanup")
    for row in internal_consistency_cleanup():
        print(f"  {row['item']}")
        print(f"    old: {row['old_status']}")
        print(f"    new: {row['new_status']}")
        print(f"    reason: {row['reason']}")

    print("\n3. RFG strength scorecard")
    for row in rfg_strength_scorecard():
        print(
            f"  {row['criterion']:24s}: {row['status']} | {row['evidence']}"
        )

    print("\n4. Mass table")
    for row in prediction_table():
        print(
            f"  {row['particle']:8s}: "
            f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
            f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
            f"rel_err={row['relative_mass_error']:.3e}"
        )

    print("\n5. Open theorem items")
    for item in open_theorem_items():
        print(f"  OPEN: {item}")

    print("\n6. Strongest defensible claim")
    print(f"  {strongest_defensible_claim()}")
