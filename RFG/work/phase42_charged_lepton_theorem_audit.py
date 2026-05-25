# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 42: Charged-lepton theorem-chain audit

Purpose:
    Collect the phase37-43 chain in one place and mark what is algebraically
    closed versus what remains a physical/theorem-level assumption.

The chain:
    phase41: I3=det(B) contains the C3 strain lock E^3.
    phase39: oriented elastic triad gives C3, strain phase gives C3 -> 9 slots.
    phase40: oriented framed closure selects h=2.
    phase43: local normal-form Hessian gives h=2 stability conditions.
    phase38: theta = h/9 = 2/9.
    phase37: C3 Koide operator with theta=2/9 gives charged-lepton ratios.

Status:
    Strong candidate derivation. The remaining open theorem is to derive the
    charged oriented-framed defect sector and its reduced holonomy directly
    as a stationary sector of the full RFG action.
"""

from phase37_c3_koide_operator import (
    LEPTON_MASSES_MEV,
    THETA_TOPOLOGICAL,
    koide_identity,
    prediction_table,
)
from phase38_z9_theta_holonomy import derivation_summary
from phase39_action_symmetry_to_z9 import z9_closure_from_action_symmetry
from phase40_projective_spinor_h2 import h2_selection_summary
from phase41_action_normal_form_theta import (
    c3_complex_strain_identity,
    invariants_around_isotropic_background,
    polar_triaxial_parameterization,
)


def pass_fail(condition):
    return "PASS" if condition else "FAIL"


def audit_action_c3_lock():
    inv = invariants_around_isotropic_background()
    ident = c3_complex_strain_identity()
    polar = polar_triaxial_parameterization()
    checks = {
        "trQ3_equals_3detQ": inv["trQ3_minus_3detQ"] == 0,
        "E3_identity_det": ident["E3_plus_Ebar3_minus_27detQ"] == 0,
        "E3_identity_trQ3": ident["E3_plus_Ebar3_minus_9trQ3"] == 0,
        "polar_det_identity": polar["det_Q_minus_expected"] == 0,
        "polar_trQ2_identity": polar["tr_Q2_minus_expected"] == 0,
    }
    return checks


def audit_z9_and_h2():
    closure = z9_closure_from_action_symmetry()
    h2 = h2_selection_summary()
    theta = h2["theta"]
    checks = {
        "axis_order_is_3": closure["axis_order"] == 3,
        "phase_order_is_3": closure["phase_order"] == 3,
        "closure_slots_is_9": closure["closure_slots"] == 9,
        "h_selected_is_2": h2["h_selected"] == 2,
        "theta_is_2_over_9": abs(theta - THETA_TOPOLOGICAL) < 1.0e-15,
    }
    return checks


def audit_mass_predictions():
    rows = prediction_table()
    checks = {}
    for row in rows:
        name = row["particle"]
        tolerance = 1.0e-4 if name != "electron" else 1.0e-15
        checks[f"{name}_mass_error_ok"] = abs(row["relative_mass_error"]) <= tolerance
    checks["koide_exact"] = abs(koide_identity() - 2.0 / 3.0) < 1.0e-15
    return checks


def open_theorem_items():
    return [
        "Derive the charged oriented-frame requirement from the RFG coupling, not as a postulate.",
        "Derive the reduced holonomy coordinate theta used by the C3 stiffness operator from the defect moduli space.",
        "Upgrade phase43 local normal-form stability to a full 3D fluctuation-operator/PDE stability proof.",
        "Show no lower-energy non-leptonic defect branch has the same charge and lower action.",
        "Extend the proof beyond principal-axis/algebraic normal form to full 3D localized oscillon fields.",
    ]


def theorem_chain_summary():
    summary = derivation_summary()
    h2 = h2_selection_summary()
    closure = z9_closure_from_action_symmetry()
    return {
        "action_lock": "I3=det(B) -> det(Q)=Re(E^3)/27 -> C3 strain locking",
        "axis_c3": closure["axis_c3_origin"],
        "phase_c3": closure["phase_c3_origin"],
        "closure_slots": closure["closure_slots"],
        "h_selected": h2["h_selected"],
        "local_stability": "phase43 Hessian conditions",
        "theta": summary["theta"],
        "koide": koide_identity(),
        "electron_anchor_MeV": LEPTON_MASSES_MEV["electron"],
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 42: Charged-lepton theorem-chain audit")
    print("=" * 72)

    print("\n1. Theorem-chain summary")
    for key, value in theorem_chain_summary().items():
        if isinstance(value, float):
            print(f"  {key:22s}: {value:.12g}")
        else:
            print(f"  {key:22s}: {value}")

    print("\n2. Algebraic C3-lock checks")
    for key, value in audit_action_c3_lock().items():
        print(f"  {key:28s}: {pass_fail(value)}")

    print("\n3. Z9 and h=2 checks")
    for key, value in audit_z9_and_h2().items():
        print(f"  {key:28s}: {pass_fail(value)}")

    print("\n4. Mass/Koide checks")
    for key, value in audit_mass_predictions().items():
        print(f"  {key:28s}: {pass_fail(value)}")

    print("\n5. Predicted charged-lepton masses")
    for row in prediction_table():
        print(
            f"  {row['particle']:8s}: "
            f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
            f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
            f"rel_err={row['relative_mass_error']:.3e}"
        )

    print("\n6. Still open theorem items")
    for item in open_theorem_items():
        print(f"  OPEN: {item}")
