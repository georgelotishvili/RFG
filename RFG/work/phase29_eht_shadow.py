# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 29: EHT shadow — Sgr A* და M87*
================================================================================

რეფერენცია: phase18_bh_singularity.py, STRATEGY.md ეტაპი III §2

დაკვირვება:
- M87* shadow (EHT 2019): θ ≈ 42 ± 3 μas, M ≈ 6.5 × 10⁹ M_sun
- Sgr A* shadow (EHT 2022): θ ≈ 51.8 ± 2.3 μas
- Sgr A* GRAVITY prior: M ≈ 4.154 × 10⁶ M_sun, D ≈ 8178 pc
- ringdown QNM ფარგლი (GW190521)
- photon sphere r_ph = 3M/2 Schwarzschild-ში

RFG-ის ცდა:
- phase18_bh_singularity: regular center A=1+a₂r², B=B₀+b₂r²
- de Sitter ბირთვი — postulated, არა derived
- Kretschmann K=R·R არ მოწმდება
- horizon-ის ცდის სქელეტი

ცდის ფანჯარა:
- GR Schwarzschild: shadow θ = 3√3 · GM/(c²·D) · (1+correction)
- RFG-ის modification: horizon-ის ცვლა → shadow size shift
- ngEHT (2030+): θ measurement to <1 μas → potential RFG-სპეციფიკური signature
"""

import math


M_SUN = 1.98847e30  # kg
G = 6.67430e-11
C = 299792458.0
MUAS_TO_RAD = math.pi / (180 * 3600 * 1e6)
PC = 3.0856775814913673e16
KPC = 1.0e3 * PC
MPC = 1.0e6 * PC


EHT_OBSERVATIONS = {
    "M87": {
        "shadow_diameter_uas": 42.0,
        "shadow_error_uas": 3.0,
        "mass_solar": 6.5e9,
        "distance_Mpc": 16.8,
    },
    "SgrA": {
        "shadow_diameter_uas": 51.8,
        "shadow_error_uas": 2.3,
        "mass_solar": 4.154e6,
        "distance_pc": 8178,
        "mass_distance_source": "GRAVITY Collaboration priors (2019/2022)",
    },
}


def schwarzschild_shadow_prediction(M_solar, D_meters):
    """
    GR Schwarzschild shadow: θ = (3√3) · r_s / D
    r_s = 2GM/c²
    """
    M = M_solar * M_SUN
    r_s = 2 * G * M / C**2
    shadow_size = 3 * math.sqrt(3) * r_s
    theta_rad = shadow_size / D_meters
    theta_uas = theta_rad / MUAS_TO_RAD
    return r_s, shadow_size, theta_uas


def compare_with_observation(name, M_solar, D_meters, theta_obs, theta_err):
    """GR Schwarzschild შედარება EHT-ის დაკვირვებასთან."""
    r_s, size, theta_pred = schwarzschild_shadow_prediction(M_solar, D_meters)
    sigma_deviation = abs(theta_pred - theta_obs) / theta_err
    return {
        "name": name,
        "r_s_meters": r_s,
        "shadow_size_meters": size,
        "GR_prediction_uas": theta_pred,
        "EHT_observation_uas": theta_obs,
        "EHT_error_uas": theta_err,
        "deviation_sigma": sigma_deviation,
    }


def distance_to_meters(obs):
    """Convert observation distance fields to meters."""
    if "distance_Mpc" in obs:
        return obs["distance_Mpc"] * MPC
    if "distance_kpc" in obs:
        return obs["distance_kpc"] * KPC
    if "distance_pc" in obs:
        return obs["distance_pc"] * PC
    raise KeyError("distance field not found")


def rfg_modifications_open():
    """RFG-სპეციფიკური modifications — ღია ცდები."""
    return [
        "phase18_bh_singularity-ის რეგულარული A, B → horizon-ის ცვლა?",
        "RFG bi-conformal e^φ/e^(-φ) → photon path equation modification",
        "Kretschmann K(r) Schwarzschild-დან გადახრა — ცდის ღია",
        "Shadow diameter dependence on a₂, b₂ (regular center parameters)",
        "ngEHT 2030+ θ precision <1 μas → RFG vs GR comparison window",
        "Photon ring substructure (n=1, n=2) — ngEHT/BHEX next-gen test",
    ]


def predictions_summary():
    """RFG vs GR — სად შეიძლება დევიაცია."""
    return {
        "current_status": "RFG = GR shadow size (phase18 ღია)",
        "ngEHT_window": "<1 μas precision — RFG modification არის ცდის ფანჯარა",
        "BHEX_window": "photon ring substructure — RFG regular center signature?",
        "current_test": "EHT 2019/2022 within ~10% of GR — RFG consistent",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 29: EHT shadow — Sgr A* და M87*")
    print("რეფერენცია: EHT 2019/2022, phase18_bh_singularity, STRATEGY ეტაპი III")
    print("=" * 72)

    print("\n1. დაკვირვება (EHT)")
    for name, obs in EHT_OBSERVATIONS.items():
        print(f"\n  {name}")
        for key, val in obs.items():
            print(f"    {key:25s}: {val}")

    print("\n2. GR Schwarzschild prediction vs EHT")
    # M87
    M87 = EHT_OBSERVATIONS["M87"]
    D_M87 = distance_to_meters(M87)
    res_M87 = compare_with_observation(
        "M87*", M87["mass_solar"], D_M87,
        M87["shadow_diameter_uas"], M87["shadow_error_uas"]
    )
    # SgrA
    SgrA = EHT_OBSERVATIONS["SgrA"]
    D_SgrA = distance_to_meters(SgrA)
    res_SgrA = compare_with_observation(
        "Sgr A*", SgrA["mass_solar"], D_SgrA,
        SgrA["shadow_diameter_uas"], SgrA["shadow_error_uas"]
    )

    for res in [res_M87, res_SgrA]:
        print(f"\n  {res['name']}")
        print(f"    r_s = {res['r_s_meters']:.3e} m")
        print(f"    GR prediction = {res['GR_prediction_uas']:.2f} μas")
        print(f"    EHT observed  = {res['EHT_observation_uas']:.1f} ± {res['EHT_error_uas']:.1f} μas")
        print(f"    deviation     = {res['deviation_sigma']:.2f} σ")

    print("\n3. RFG modifications — ღია ცდები")
    for i, task in enumerate(rfg_modifications_open(), 1):
        print(f"  {i}. {task}")

    print("\n4. Predictions summary")
    for key, val in predictions_summary().items():
        print(f"  {key:18s}: {val}")

    print("\n5. სტატუსი")
    print("  - GR shadow consistent with EHT within ~10% (M87, SgrA)")
    print("  - RFG-სპეციფიკური modification ღია ცდაა (phase18 regularity)")
    print("  - ngEHT 2030+ <1 μas precision — RFG vs GR comparison window")
