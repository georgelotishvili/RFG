# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 27: Double Pulsar PSR J0737-3039A/B — PPK ცდები
================================================================================

რეფერენცია: STRATEGY.md ეტაპი III §2, RFG_Theory.md § 6

დაკვირვება — PSR J0737-3039A/B (Kramer et al. 2021, PRX 11:041050):
- ორმაგი პულსარი (ერთადერთი ცნობილი)
- 6 post-Keplerian (PPK) პარამეტრი 0.05% precision
- GR ცდის ერთ-ერთი მკაცრი ფარგლი

PPK პარამეტრები:
1. ω̇ (periastron advance) — RFG: PPN γ, β-დან
2. γ (Einstein delay) — gravitational redshift + time dilation
3. P_b_dot (orbital period decay) — quadrupole + dipole radiation
4. r (Shapiro range) — phase13 1PN-დან
5. s (Shapiro shape) — sin(inclination)
6. Ω̇ (geodetic precession) — phase30 Lense-Thirring related

RFG-ის პროგნოზი:
- 1PN: GR-ის იდენტური (PPN γ=β=1)
- 2.5PN orbital decay: ცდის სქელეტი (phase26)
- Scalar dipole: phase9 ღია (s_A = 0 postulated)

ცდის ფარგლი GR vs RFG:
- 1PN ემთხვევა
- 2.5PN ცდის სქელეტი
- Scalar dipole — RFG-ის ღია ცდა
"""

import math


PSR_J0737_DATA = {
    "reference": "Kramer et al. 2021, PRX 11:041050",
    "P_b_orbital_period_days": 0.10225156248,
    "eccentricity": 0.0877775,
    "M_A_total_solar": 1.338185,
    "M_B_total_solar": 1.248868,
    "periastron_advance_deg_yr": 16.899323,  # ω̇
    "Einstein_delay_ms": 0.384045,  # γ
    "P_b_dot_obs": -1.247920e-12,  # observed orbital decay
    "P_b_dot_GR_pred": -1.247843e-12,  # GR quadrupole prediction
    "Shapiro_r_M_sun_s": 6.162e-6,  # r in solar mass · seconds
    "Shapiro_s": 0.999936,  # sin(i)
    "geodetic_omega_deg_yr": 4.78,  # B's geodetic precession
}


def gr_predictions():
    """GR-ის ცხადი 1PN პრედიქცია PSR J0737-სთვის."""
    G = 6.674e-11
    c = 2.998e8
    M_sun = 1.989e30
    M_A = PSR_J0737_DATA["M_A_total_solar"] * M_sun
    M_B = PSR_J0737_DATA["M_B_total_solar"] * M_sun
    M_tot = M_A + M_B
    P_b = PSR_J0737_DATA["P_b_orbital_period_days"] * 86400
    e = PSR_J0737_DATA["eccentricity"]

    # ω̇_GR = 3 · (GM/c²a)^(5/3) · (2π/P_b)^(5/3) / (1-e²)
    # სიმარტივისთვის — სიმბოლურად
    n = 2 * math.pi / P_b
    omega_dot_GR_rad = 3 * n ** (5 / 3) * (G * M_tot / c**3) ** (2 / 3) / (1 - e**2)
    omega_dot_GR_deg_yr = omega_dot_GR_rad * (180 / math.pi) * 86400 * 365.25

    return {
        "M_total_kg": M_tot,
        "P_b_seconds": P_b,
        "omega_dot_GR_deg_per_yr": omega_dot_GR_deg_yr,
        "omega_dot_observed": PSR_J0737_DATA["periastron_advance_deg_yr"],
    }


def rfg_predictions():
    """RFG-ის პრედიქცია — phase8-დან γ=β=1 ⇒ იდენტური GR."""
    return {
        "PPN_gamma": 1.0,  # phase8 (RFG bi-conformal)
        "PPN_beta": 1.0,  # phase8 2PN
        "omega_dot_RFG": "იდენტური GR-ის (γ=β=1)",
        "Einstein_delay_RFG": "იდენტური GR-ის (Pound-Rebka + gravitational time dilation)",
        "P_b_dot_quadrupole": "იდენტური GR-ის leading order (phase9 c_T=c)",
        "P_b_dot_dipole_RFG": "phase9-დან s_A = 0 (postulated) ⇒ no dipole",
        "P_b_dot_dipole_status": "OPEN — strong-field s_A derivation აკლია (Damour-Esposito-Farèse)",
    }


def open_tests():
    """ცდის ღია ნაბიჯები."""
    return [
        "Strong-field s_A derivation Damour-Esposito-Farèse სცენარით",
        "P_b_dot_dipole_RFG რიცხობრივი ცდა |α_0| < 2e-5 (PSR J1738+0333) ფარგლი",
        "Geodetic precession Ω̇ — Lense-Thirring (phase30) გადახედვა",
        "ცდის ცხადი interface PTA კოლაბორაციით (NANOGrav, EPTA)",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 27: Double Pulsar PSR J0737-3039A/B")
    print("რეფერენცია: Kramer 2021, PRX 11:041050; phase8, phase9, RFG_Theory § 6")
    print("=" * 72)

    print("\n1. დაკვირვება (PSR J0737-3039, Kramer et al. 2021)")
    for key, val in PSR_J0737_DATA.items():
        print(f"  {key:30s}: {val}")

    print("\n2. GR-ის 1PN პრედიქცია")
    gr = gr_predictions()
    for key, val in gr.items():
        if isinstance(val, float):
            print(f"  {key:30s}: {val:.6e}")
        else:
            print(f"  {key:30s}: {val}")

    print("\n3. RFG-ის პრედიქცია (phase8 + phase9)")
    rfg = rfg_predictions()
    for key, val in rfg.items():
        print(f"  {key:30s}: {val}")

    print("\n4. ღია ცდები")
    for i, task in enumerate(open_tests(), 1):
        print(f"  {i}. {task}")

    print("\n5. სტატუსი")
    print("  - 1PN: RFG=GR ფიქსირდება (phase8 γ=β=1)")
    print("  - 2.5PN orbital decay: GR=RFG leading order, dipole ცდის სქელეტი")
    print("  - Strong-field s_A — ღია (Damour-Esposito-Farèse)")
    print("  - PSR J0737-3039 ცდის სრული χ² fit — ღია")
