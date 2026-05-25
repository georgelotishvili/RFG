# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 29: EHT shadow - old compact-object prediction recovered in RFG.

Phase 18 now derives the exponential exterior:

    ds^2 = -exp(-r_s/r)c^2dt^2 + exp(r_s/r)(dr^2+r^2dOmega^2).

For null circular orbits:

    d/dr [exp(-2r_s/r)/r^2] = 0 -> r_ph = r_s.

The critical impact parameter is therefore:

    b_c^RFG = e*r_s,

while Schwarzschild gives:

    b_c^GR = (3*sqrt(3)/2)*r_s.

Thus the static spherical RFG benchmark predicts a shadow diameter
larger by 2e/(3sqrt(3))-1 = 4.63%.
"""

import math


M_SUN = 1.98847e30
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


def schwarzschild_shadow_prediction(m_solar, distance_m):
    """GR Schwarzschild critical-curve diameter: theta = 3*sqrt(3)*r_s/D."""
    mass = m_solar * M_SUN
    r_s = 2 * G * mass / C**2
    diameter = 3 * math.sqrt(3) * r_s
    theta_uas = (diameter / distance_m) / MUAS_TO_RAD
    return {
        "r_s_meters": r_s,
        "critical_impact_parameter_m": 0.5 * diameter,
        "shadow_diameter_meters": diameter,
        "theta_uas": theta_uas,
    }


def rfg_shadow_prediction(m_solar, distance_m):
    """RFG exponential-exterior critical-curve diameter: theta = 2*e*r_s/D."""
    mass = m_solar * M_SUN
    r_s = 2 * G * mass / C**2
    b_c = math.e * r_s
    diameter = 2 * b_c
    theta_uas = (diameter / distance_m) / MUAS_TO_RAD
    return {
        "r_s_meters": r_s,
        "photon_sphere_r": r_s,
        "critical_impact_parameter_m": b_c,
        "shadow_diameter_meters": diameter,
        "theta_uas": theta_uas,
    }


def compare_with_observation(name, m_solar, distance_m, theta_obs, theta_err):
    """Compare GR and RFG static spherical shadow benchmarks to one observation."""
    gr = schwarzschild_shadow_prediction(m_solar, distance_m)
    rfg = rfg_shadow_prediction(m_solar, distance_m)
    ratio = rfg["theta_uas"] / gr["theta_uas"]

    return {
        "name": name,
        "r_s_meters": gr["r_s_meters"],
        "GR_prediction_uas": gr["theta_uas"],
        "RFG_prediction_uas": rfg["theta_uas"],
        "RFG_over_GR_ratio": ratio,
        "RFG_shadow_shift_percent": (ratio - 1.0) * 100.0,
        "EHT_observation_uas": theta_obs,
        "EHT_error_uas": theta_err,
        "GR_deviation_sigma": abs(gr["theta_uas"] - theta_obs) / theta_err,
        "RFG_deviation_sigma": abs(rfg["theta_uas"] - theta_obs) / theta_err,
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


def rfg_shadow_derivation_ledger():
    return [
        "exponential exterior: g_tt=-exp(-r_s/r), g_rr=exp(r_s/r)",
        "null barrier: V_null proportional to exp(-2r_s/r)/r^2",
        "photon sphere: dV_null/dr=0 -> r_ph=r_s",
        "critical impact parameter: b_c=r*exp(r_s/r) at r=r_s -> e*r_s",
        "GR reference: b_c=(3*sqrt(3)/2)*r_s",
        "static spherical prediction: RFG shadow diameter is +4.63% relative to GR",
    ]


def predictions_summary():
    """RFG vs GR shadow status."""
    ratio = 2.0 * math.e / (3.0 * math.sqrt(3.0))
    return {
        "current_status": "derived static spherical RFG benchmark, not open",
        "RFG_b_c": "e*r_s",
        "GR_b_c": "3*sqrt(3)*r_s/2",
        "RFG_over_GR": ratio,
        "shift_percent": (ratio - 1.0) * 100.0,
        "needed_for_decisive_test": "spin, accretion, mass-distance priors, and ray-traced image modelling",
        "ngEHT_BHEX_window": "few-percent shadow/ring precision can test the +4.63% benchmark",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 29: EHT shadow - RFG b_c=e*r_s benchmark")
    print("=" * 72)

    print("\n1. დაკვირვება (EHT priors used in this local script)")
    for name, obs in EHT_OBSERVATIONS.items():
        print(f"\n  {name}")
        for key, val in obs.items():
            print(f"    {key:25s}: {val}")

    print("\n2. Derivation ledger")
    for item in rfg_shadow_derivation_ledger():
        print(f"  - {item}")

    print("\n3. GR vs RFG static benchmark compared to EHT numbers")
    for name, obs in EHT_OBSERVATIONS.items():
        result = compare_with_observation(
            name,
            obs["mass_solar"],
            distance_to_meters(obs),
            obs["shadow_diameter_uas"],
            obs["shadow_error_uas"],
        )
        print(f"\n  {result['name']}")
        print(f"    r_s = {result['r_s_meters']:.3e} m")
        print(f"    GR theta  = {result['GR_prediction_uas']:.2f} microas")
        print(f"    RFG theta = {result['RFG_prediction_uas']:.2f} microas")
        print(f"    RFG/GR    = {result['RFG_over_GR_ratio']:.8f}")
        print(f"    shift     = {result['RFG_shadow_shift_percent']:.2f}%")
        print(f"    observed  = {result['EHT_observation_uas']:.1f} +/- {result['EHT_error_uas']:.1f} microas")
        print(f"    GR sigma  = {result['GR_deviation_sigma']:.2f}")
        print(f"    RFG sigma = {result['RFG_deviation_sigma']:.2f}")

    print("\n4. Predictions summary")
    for key, val in predictions_summary().items():
        print(f"  {key:26s}: {val}")

    print("\n5. სტატუსი")
    print("  - ძველი +4.6% shadow პროგნოზი ახალ RFG phase18-დან ზუსტად გამოდის.")
    print("  - მიმდინარე EHT რიცხვები არ არის საკმარისი სუფთა GR/RFG გარჩევისთვის.")
    print("  - decisive test მოითხოვს rotating RFG ray tracing-ს და ngEHT/BHEX კლასის სიზუსტეს.")
