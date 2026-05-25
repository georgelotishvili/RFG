# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 47: Particle-sector falsifiability map after the C3/Z9 upgrade

Purpose:
    Convert the strengthened particle sector into concrete pass/fail claims.

Important cleanup:
    The main charged-lepton route is now phase37-43 C3/Z9, not the older
    radial N-ladder. Therefore:

        - e, mu, tau are one C3 triplet.
        - missing radial modes n=2..13 are no longer expected generations.
        - N=4 329 keV is not a main prediction under the C3 route.

This makes the particle sector cleaner and easier to falsify.
"""

from phase37_c3_koide_operator import prediction_table
from phase43_h2_branch_stability import stability_conditions
from phase45_radiative_stability_requirement import rfg_radiative_options
from phase46_quark_extension_candidate import candidate_tables


def charged_lepton_pass_fail():
    return {
        "claim": "e, mu, tau are one C3/Z9 triplet with theta=2/9",
        "status": "postdiction unless derived before using the masses",
        "pass_condition": "phase37 mass ratios remain within pole-mass uncertainties after radiative protection is derived",
        "fail_condition": "full RFG action cannot support a stable h=2 C3 triplet or radiative corrections destroy the pole relation",
        "numbers": [
            {
                "particle": row["particle"],
                "predicted_MeV": row["predicted_mass_MeV"],
                "observed_MeV": row["observed_mass_MeV"],
                "relative_error": row["relative_mass_error"],
            }
            for row in prediction_table()
        ],
    }


def deprecated_legacy_predictions():
    return [
        {
            "legacy_claim": "N=5,72,295 radial ladder",
            "new_status": "not the main lepton-generation mechanism",
            "why": "C3 triplet explains all three charged leptons together",
        },
        {
            "legacy_claim": "N=4 329 keV lepton-like companion",
            "new_status": "suspended under C3 route",
            "why": "C3 charged-lepton theorem does not imply a radial N=4 state",
        },
        {
            "legacy_claim": "muon = uniquely n=14 overtone",
            "new_status": "legacy audit only",
            "why": "overlap integrals were monotonic; no unique n=14 bottleneck",
        },
    ]


def new_falsifiable_targets():
    quark_candidates = candidate_tables()
    return [
        {
            "target": "radiative protection",
            "test": "derive either dressed pole-frequency theorem or loop cancellation",
            "fail": "unprotected QED-like generation logs shift Koide by ~1e-3",
            "file": "phase45_radiative_stability_requirement.py",
        },
        {
            "target": "h=2 local stability",
            "test": "satisfy phase43 Hessian inequalities and then full fluctuation spectrum",
            "fail": "h=2 branch is a saddle or has negative fluctuation mode",
            "file": "phase43_h2_branch_stability.py",
        },
        {
            "target": "quark zero-angle extension",
            "test": "one of two projection candidates survives fixed-scheme quark masses without fitting theta_U/D",
            "fail": "both projections fail once MSbar masses are run to the chosen scale",
            "file": "phase46_quark_extension_candidate.py",
            "candidate_count": len(quark_candidates),
        },
        {
            "target": "no C3-forbidden charged h=1 branch",
            "test": "charged defects require oriented framed closure, not projective/nematic closure",
            "fail": "a stable charged h=1 branch appears in the full defect spectrum",
            "file": "phase40_projective_spinor_h2.py",
        },
    ]


def next_calculation_queue():
    return [
        "Build the full localized h=2 oscillon ansatz: Phi(t,r), Q(r,beta), framed triad U(x).",
        "Linearize the RFG action around that ansatz and compute the fluctuation operator.",
        "Check whether all non-gauge eigenvalues are non-negative.",
        "Add EM/self-field coupling and decide between phase45 options A or B.",
        "Choose one quark-mass scheme/scale and run phase46 as a real numerical falsification test.",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 47: Particle-sector falsifiability map")
    print("=" * 72)

    print("\n1. Charged-lepton C3 claim")
    cl = charged_lepton_pass_fail()
    print(f"  claim: {cl['claim']}")
    print(f"  status: {cl['status']}")
    print(f"  pass: {cl['pass_condition']}")
    print(f"  fail: {cl['fail_condition']}")
    for row in cl["numbers"]:
        print(
            f"    {row['particle']:8s}: "
            f"pred={row['predicted_MeV']:.6f} MeV, "
            f"obs={row['observed_MeV']:.6f} MeV, "
            f"rel_err={row['relative_error']:.3e}"
        )

    print("\n2. Legacy predictions demoted")
    for row in deprecated_legacy_predictions():
        print(f"  {row['legacy_claim']}: {row['new_status']} ({row['why']})")

    print("\n3. New falsifiable targets")
    for row in new_falsifiable_targets():
        extra = f", candidates={row['candidate_count']}" if "candidate_count" in row else ""
        print(f"  {row['target']} [{row['file']}{extra}]")
        print(f"    test: {row['test']}")
        print(f"    fail: {row['fail']}")

    print("\n4. h=2 stability conditions")
    for key, value in stability_conditions().items():
        print(f"  {key:26s}: {value}")

    print("\n5. Radiative protection options")
    for row in rfg_radiative_options():
        print(f"  {row['option']}: {row['claim_needed']}")

    print("\n6. Next calculation queue")
    for item in next_calculation_queue():
        print(f"  - {item}")
