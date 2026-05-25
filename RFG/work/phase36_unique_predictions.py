# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 36: RFG-სპეციფიკური უნიკალური ფალსიფიცირებადი პროგნოზები
================================================================================

სტატუსი:
Strategy 3 / X3-ის შესრულება. ეს ფაილი არჩევს მხოლოდ იმ პროგნოზებს, რომლებიც
RFG-ს რეალურად გამოარჩევს GR/MOND/AQUAL/TeVeS/RMOND/AeST-გან.

მკაცრი კლასიფიკაცია:
    C1. N=4 lepton-like resonance: RFG-specific მხოლოდ მაშინ, თუ phase35
        N-ladder-ს ფესვიდან გამოიყვანს. ამ ეტაპზე პირობითია.
    C2. a_0(z): shared/phenomenological; MOND-family-სგან ბოლომდე არ გამოარჩევს.
    C3. photon ring: shared with many regular-BH models.
    C4. scalar breathing: shared with scalar-tensor theories.
    C5. framed vortex in dipolar supersolid: lab-specific RFG signature.

ამიტომ phase36-ის მკაცრი deliverable არის ორი ფანჯარა:
    1. C1: N=4 = 329 ± 2 keV, plus X17/ATOMKI consistency check.
    2. C5: Eu/Dy framed-vortex protocol with concrete labs.

განახლება phase37-44 შემდეგ:
    charged-lepton sector-ის მთავარი გზა გახდა C3 Koide operator და Z9
    framed holonomy. ამიტომ C1/N=4 აღარ არის მთავარი particle-sector
    პროგნოზი; ის რჩება მხოლოდ legacy/conditional branch-ად, თუ ძველი
    radial N-ladder ოდესმე დამოუკიდებლად დამტკიცდება.
"""

import math


M_E_MEV = 0.51099895000
N_ELECTRON = 5
N4_PREDICTED_MEV = 0.329
N4_THEORY_SIGMA_MEV = 0.002
X17_MEV = 17.0


def n_ladder_mass(n_value):
    """Existing ladder relation m_N = m_e * (N/N_e)^2."""
    return M_E_MEV * (n_value / N_ELECTRON) ** 2


def n4_fixed_prediction():
    """
    N=4 fixed number.

    q=0 ladder gives 327.04 keV. The theory text uses q-corrected 329 keV.
    Until phase35 derives the correction, ±2 keV is a theory-band, not a
    statistical experimental sigma.
    """
    q0_mev = n_ladder_mass(4)
    center_mev = N4_PREDICTED_MEV
    sigma_mev = N4_THEORY_SIGMA_MEV
    return {
        "name": "C1_N4_lepton_like_resonance",
        "n": 4,
        "q0_mass_keV": q0_mev * 1000,
        "prediction_keV": center_mev * 1000,
        "theory_sigma_keV": sigma_mev * 1000,
        "band_keV": ((center_mev - sigma_mev) * 1000, (center_mev + sigma_mev) * 1000),
        "rfg_uniqueness": "conditional: unique only if phase35 derives N-ladder",
        "phase37_update": "legacy/suspended under C3 charged-lepton route",
        "falsification": "dedicated null searches excluding 327-331 keV lepton-like resonance at required coupling",
    }


def x17_atomki_consistency():
    """
    ATOMKI/X17 is a 17 MeV anomaly, not the same object as N=4.

    The only non-fake bridge available here is a harmonic/transition audit:
    does X17 sit near an integer multiple of the 329 keV N=4 scale?
    """
    ratio = X17_MEV / N4_PREDICTED_MEV
    nearest_integer = round(ratio)
    harmonic_mev = nearest_integer * N4_PREDICTED_MEV
    fractional_error = abs(harmonic_mev - X17_MEV) / X17_MEV

    return {
        "x17_mass_MeV": X17_MEV,
        "n4_mass_MeV": N4_PREDICTED_MEV,
        "x17_over_n4": ratio,
        "nearest_harmonic": nearest_integer,
        "harmonic_mass_MeV": harmonic_mev,
        "fractional_error": fractional_error,
        "verdict": "weak/conditional bridge; do not identify X17 with N=4",
        "needed": "derive an RFG transition/coupling rule before claiming an X17 explanation",
    }


def framed_vortex_protocol():
    """
    RFG-specific laboratory protocol.

    Framed vortex prediction: vortex line has a measurable framing vector.
    RFG expects half-integer twist classes, Tw = n/2, stable under smooth
    deformations unless a reconnection event changes linking.
    """
    labs = [
        {
            "lab": "MIT / Ketterle",
            "system": "ultracold BEC/supersolid platforms",
            "role": "vortex preparation and interferometric readout candidate",
        },
        {
            "lab": "Stuttgart / Pfau",
            "system": "dipolar dysprosium supersolid",
            "role": "density-modulated supersolid and vortex protocol candidate",
        },
        {
            "lab": "Innsbruck / Ferlaino",
            "system": "Er/Dy dipolar supersolid, observed vortices in supersolid phase",
            "role": "closest existing platform for framed-vortex test",
            "strategy_label": "Strategy 3 wording: Pfau Innsbruck; treated here as the Innsbruck dipolar-supersolid platform",
        },
    ]

    return {
        "name": "C5_Eu_Dy_framed_vortex",
        "observable": "twist/framing class Tw of vortex core",
        "prediction": "Tw in half-integer classes: ..., -1, -1/2, 0, 1/2, 1, ...",
        "null_result": "continuous unquantized framing distribution -> RFG framed-vortex claim fails",
        "positive_result": "stable Tw=1/2 family -> RFG-specific lab signature",
        "labs": labs,
    }


def shared_prediction_audit():
    """Candidate-ები, რომლებიც სასარგებლოა, მაგრამ RFG-ს უნიკალურად არ არჩევს."""
    return [
        {
            "name": "C2_a0_redshift",
            "classification": "shared/phenomenological",
            "why_not_unique": "MOND-family models can also parametrize a0(z).",
            "keep_as": "cosmology consistency test, not unique signature",
        },
        {
            "name": "C3_photon_ring_substructure",
            "classification": "shared",
            "why_not_unique": "many regular-BH or quantum-corrected metrics modify rings.",
            "keep_as": "black-hole consistency test",
        },
        {
            "name": "C4_scalar_breathing",
            "classification": "shared",
            "why_not_unique": "scalar-tensor theories also predict breathing modes.",
            "keep_as": "polarization bound on scalar sector",
        },
    ]


def unique_scorecard():
    c1 = n4_fixed_prediction()
    c5 = framed_vortex_protocol()
    return [
        {
            "candidate": "C1 N=4 329 keV",
            "rfg_specific": "legacy conditional",
            "time_window": "now, if coupling/search channel specified",
            "falsifiable_number": f"{c1['prediction_keV']:.0f} ± {c1['theory_sigma_keV']:.0f} keV",
        },
        {
            "candidate": "C5 framed vortex",
            "rfg_specific": "yes",
            "time_window": "now/near-term cold atom labs",
            "falsifiable_number": "Tw = n/2 classes",
        },
    ]


def status_assessment():
    return {
        "closed_now": "5 candidates no longer presented as equally unique; C5 remains the clean RFG lab signature.",
        "c1_status": "legacy conditional after phase37-44; not a main prediction unless radial N-ladder is revived.",
        "c5_status": "concrete labs and null/positive outcomes added.",
        "dependency": "particle-sector main route is now phase37-44 C3/Z9, not N-ladder.",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 36: უნიკალური RFG ფალსიფიკატორები")
    print("=" * 72)

    print("\n1. C1 fixed N=4 prediction")
    for key, value in n4_fixed_prediction().items():
        print(f"  {key:22s}: {value}")

    print("\n2. ATOMKI/X17 consistency audit")
    for key, value in x17_atomki_consistency().items():
        if isinstance(value, float):
            print(f"  {key:22s}: {value:.6g}")
        else:
            print(f"  {key:22s}: {value}")

    print("\n3. C5 framed-vortex protocol")
    c5 = framed_vortex_protocol()
    for key, value in c5.items():
        if key == "labs":
            print("  labs:")
            for lab in value:
                print(f"    - {lab['lab']}: {lab['system']} ({lab['role']})")
                if "strategy_label" in lab:
                    print(f"      note: {lab['strategy_label']}")
        else:
            print(f"  {key:22s}: {value}")

    print("\n4. Shared predictions kept as secondary tests")
    for item in shared_prediction_audit():
        print(f"  {item['name']:25s}: {item['classification']} | {item['why_not_unique']}")

    print("\n5. Unique scorecard")
    for item in unique_scorecard():
        print(
            f"  {item['candidate']:20s} | specific={item['rfg_specific']:11s} "
            f"| number={item['falsifiable_number']} | window={item['time_window']}"
        )

    print("\n6. Status")
    for key, value in status_assessment().items():
        print(f"  {key:12s}: {value}")
