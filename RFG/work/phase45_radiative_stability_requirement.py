# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 45: Radiative stability requirement for the C3 lepton sector

Status:
    Strengthening audit. This file does not solve radiative protection; it
    states the problem quantitatively and defines what RFG must prove.

Why this matters:
    Koide-type relations are usually discussed for charged-lepton pole
    masses. Ordinary QED self-energy corrections are generation-dependent
    because they contain log(m_i). Without a protection mechanism, a beautiful
    tree-level C3 relation may not survive to observed pole masses.

RFG options:
    A. Pole-frequency theorem:
       The C3 operator acts on fully dressed oscillon normal modes, so the
       observed pole masses are the eigenvalues of the already-renormalized
       medium. Then radiative corrections are not added afterward.

    B. Ward/Sumino-like cancellation:
       Phase-elastic/framed-sector loops cancel the generation-dependent QED
       logarithms, preserving the C3 ratios.

    C. Running-scale theorem:
       The C3 relation is exact at a special RFG matching scale, and its
       pole-mass accuracy is a consequence of controlled RG flow.

The strongest route for RFG is A or B. This file quantifies the required
protection as a condition on fractional frequency shifts.
"""

import math

from phase37_c3_koide_operator import (
    LEPTON_MASSES_MEV,
    c3_frequency_ratios,
    koide_ratio_from_frequencies,
)


ALPHA_EM = 1.0 / 137.035999084


def observed_frequency_vector():
    m_e = LEPTON_MASSES_MEV["electron"]
    return [
        math.sqrt(LEPTON_MASSES_MEV[name] / m_e)
        for name in ("electron", "muon", "tau")
    ]


def c3_frequency_vector():
    ratios = c3_frequency_ratios()
    return [ratios[name] for name in ("electron", "muon", "tau")]


def koide_sensitivity_to_fractional_shifts(eps):
    """
    Apply nu_i -> nu_i * (1 + eps_i) and return the Koide drift.
    eps must be a length-3 fractional-shift list.
    """
    nu = c3_frequency_vector()
    shifted = [value * (1.0 + eps_i) for value, eps_i in zip(nu, eps)]
    return koide_ratio_from_frequencies(shifted) - 2.0 / 3.0


def leading_qed_log_shifts(cutoff_mev):
    """
    Crude one-loop logarithmic mass-shift template:

        delta m_i / m_i ~ (3 alpha / 4 pi) log(cutoff^2 / m_i^2).

    Since nu ~ sqrt(m), the frequency shift is half this value.

    This is not used as a physical RFG prediction; it is a stress test.
    """
    prefactor_mass = 3.0 * ALPHA_EM / (4.0 * math.pi)
    shifts = []
    for name in ("electron", "muon", "tau"):
        mass = LEPTON_MASSES_MEV[name]
        dm_over_m = prefactor_mass * math.log((cutoff_mev / mass) ** 2)
        shifts.append(0.5 * dm_over_m)
    mean_shift = sum(shifts) / 3.0
    relative_to_mean = [shift - mean_shift for shift in shifts]
    return {
        "cutoff_MeV": cutoff_mev,
        "freq_shifts": shifts,
        "mean_shift": mean_shift,
        "generation_dependent_part": relative_to_mean,
        "koide_drift_raw": koide_sensitivity_to_fractional_shifts(shifts),
        "koide_drift_gen_dep": koide_sensitivity_to_fractional_shifts(relative_to_mean),
    }


def protection_condition_symbolic():
    """
    Linearized condition for preserving Koide.

    The common shift eps_i = const drops out. The dangerous part is the
    component of eps_i along the Koide-normal direction in frequency space.
    """
    nu = c3_frequency_vector()
    s1 = sum(nu)
    s2 = sum(value * value for value in nu)
    # dK = [2 nu_i / s1^2 - 2 s2 / s1^3] dnu_i.
    gradient = [
        2.0 * value / (s1 * s1) - 2.0 * s2 / (s1 * s1 * s1)
        for value in nu
    ]
    return {
        "linear_gradient_dK_dnu": gradient,
        "condition": "sum_i gradient_i * nu_i * eps_i = 0",
        "safe_shift": "eps_e = eps_mu = eps_tau is automatically safe",
        "dangerous_shift": "generation-dependent logarithms must be cancelled or reabsorbed into the dressed C3 operator",
    }


def rfg_radiative_options():
    return [
        {
            "option": "A. dressed pole-frequency theorem",
            "claim_needed": "C3 eigenvalues are already physical pole frequencies",
            "strength": "best fit to RFG oscillon ontology",
            "open_task": "derive dressed normal-mode equation including EM/self-field backreaction",
        },
        {
            "option": "B. Ward/Sumino-like cancellation",
            "claim_needed": "phase-elastic loops cancel QED log(m_i) terms",
            "strength": "closest to known Koide-protection logic",
            "open_task": "identify the RFG current and loop sign that cancels QED logs",
        },
        {
            "option": "C. matching-scale theorem",
            "claim_needed": "C3 relation exact at a scale and RG flow preserves pole accuracy",
            "strength": "standard EFT language",
            "open_task": "derive RG equations for the RFG lepton sector",
        },
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 45: Radiative stability requirement")
    print("=" * 72)

    print("\n1. Koide linear sensitivity")
    for key, value in protection_condition_symbolic().items():
        print(f"  {key:28s}: {value}")

    print("\n2. QED-log stress test")
    for cutoff in [1.0e3, 1.0e6, 1.0e9]:
        data = leading_qed_log_shifts(cutoff)
        print(f"  cutoff={cutoff:.3e} MeV")
        print(f"    freq shifts          : {[round(x, 8) for x in data['freq_shifts']]}")
        print(f"    generation dep. part : {[round(x, 8) for x in data['generation_dependent_part']]}")
        print(f"    Koide drift raw      : {data['koide_drift_raw']:.3e}")
        print(f"    Koide drift gen.dep. : {data['koide_drift_gen_dep']:.3e}")

    print("\n3. RFG protection options")
    for row in rfg_radiative_options():
        print(f"  {row['option']}")
        print(f"    claim needed: {row['claim_needed']}")
        print(f"    strength    : {row['strength']}")
        print(f"    open task   : {row['open_task']}")
