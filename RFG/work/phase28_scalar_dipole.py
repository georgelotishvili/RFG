# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 28: PSR J1738+0333 — Scalar dipole radiation bound
================================================================================

რეფერენცია: phase9_gravitational_waves.py, RFG_Theory.md § 7

დაკვირვება — PSR J1738+0333 (Freire et al. 2012, Antoniadis et al. 2013):
- pulsar-white dwarf binary
- highly asymmetric system (NS vs WD)
- scalar dipole-radiation-ის ცდის მკაცრი ფარგლი
- α_0 (scalar charge coupling) bound: |α_0| < 2 × 10⁻⁵ (95% CL)

scalar-tensor theory predictions:
- Brans-Dicke: dipole ∝ (α_NS - α_WD)²
- Damour-Esposito-Farèse: strong-field "spontaneous scalarization"
- RFG (phase9): s_A = 0 postulated ⇒ no dipole

RFG-ის ცდა:
1. dipole radiation: ΔP_b/P_b ∝ (s_A - s_B)² · (orbital velocity)²
2. RFG-ში s_A = ∂ln(m_A)/∂Φ
3. Minimal coupling ⇒ s_A = 0 (postulated)
4. Strong-field correction (Damour-Esposito-Farèse): s_A ~ Ω/(mc²)
5. PSR J1738 binding energy: Ω/(mc²) ~ 0.1 (NS) vs ~10^(-6) (WD)
6. ცდის ფარგლი: |α_0| ≤ 2 × 10⁻⁵
"""

import math


PSR_J1738_DATA = {
    "type": "NS + WD binary",
    "M_pulsar_solar": 1.46,
    "M_companion_solar": 0.181,
    "orbital_period_days": 0.354790739872,
    "eccentricity": 0.34e-6,  # very circular
    "P_b_dot_obs": -25.9e-15,  # observed
    "P_b_dot_GR": -27.7e-15,  # GR quadrupole prediction
    "alpha0_bound_95CL": 2e-5,  # scalar charge coupling upper bound
    "delta_alpha0_NS_WD_bound": 6e-3,  # asymmetry bound
}


def gr_quadrupole_prediction():
    """GR quadrupole-only orbital decay (no dipole)."""
    return {
        "P_b_dot_GR": PSR_J1738_DATA["P_b_dot_GR"],
        "ratio_obs_to_GR": PSR_J1738_DATA["P_b_dot_obs"] / PSR_J1738_DATA["P_b_dot_GR"],
        "comment": "PSR J1738 P_b_dot_obs/P_b_dot_GR = 0.94 ± 0.10 — GR consistent",
    }


def scalar_dipole_prediction():
    """
    Scalar-tensor dipole correction.

    ΔP_b_dot / P_b_dot_GR = -(4/3) · (G/c²) · (αA - αB)² · M_chirp² / (P_b² · G M_tot)^(2/3)

    Damour-Esposito-Farèse strong-field:
    αA = α_0 + β_0 · (Ω/mc²) + higher order
    """
    return {
        "Brans_Dicke_form": "ΔP_b ∝ (α_A - α_B)²",
        "RFG_postulated_phase9": "s_A = ∂ln(m_A)/∂Φ = 0 (minimal coupling)",
        "Damour_strong_field": "s_A = α_0 + β_0 · (Ω_A/(m_A c²))",
        "NS_binding_energy": "Ω_NS/(m_NS c²) ~ 0.1 (R~10 km, M~1.4 M_sun)",
        "WD_binding_energy": "Ω_WD/(m_WD c²) ~ 10^(-4) (R~10^4 km)",
        "asymmetry_bound": f"|s_NS - s_WD| < {PSR_J1738_DATA['delta_alpha0_NS_WD_bound']:.1e}",
        "alpha0_bound": f"|α_0| < {PSR_J1738_DATA['alpha0_bound_95CL']:.1e}",
    }


def rfg_strong_field_open():
    """RFG-ის strong-field s_A derivation — ცდის სქელეტი."""
    return [
        "Komar-integrand argument (phase9 Appendix 16) — s_A ≈ 1/2 leading order",
        "Bi-conformal weight e^(-φ) · ρ_0 — kinematic redshift + spatial volume",
        "Structural-response correction (open task per phase9)",
        "Strong-field NS: Ω/mc² ~ 0.1 — non-perturbative regime",
        "PSR J1738 χ² fit RFG s_A-დან — ცდა ღია",
        "Future: ngVLA + SKA pulsar timing → 100× precision",
    ]


def falsification_window():
    """RFG-ის ფალსიფიკაციის ფანჯარა PSR J1738-ში."""
    return {
        "current_bound": "|α_0| < 2 × 10⁻⁵ (PSR J1738)",
        "RFG_postulated": "s_A = 0 (consistent with bound)",
        "RFG_derived_value": "OPEN — strong-field calculation needed",
        "if_s_A_nonzero_NS_O(0.01)": "FALSIFIED (RFG αA-αB ≫ 2e-5)",
        "if_s_A_zero_strictly": "RFG ემთხვევა, არ ფალსიფიცირდება",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 28: PSR J1738+0333 — Scalar dipole bound")
    print("რეფერენცია: Freire 2012, Antoniadis 2013, phase9")
    print("=" * 72)

    print("\n1. დაკვირვება (PSR J1738+0333)")
    for key, val in PSR_J1738_DATA.items():
        print(f"  {key:30s}: {val}")

    print("\n2. GR quadrupole prediction")
    gr = gr_quadrupole_prediction()
    for key, val in gr.items():
        print(f"  {key:25s}: {val}")

    print("\n3. Scalar dipole prediction")
    dipole = scalar_dipole_prediction()
    for key, val in dipole.items():
        print(f"  {key:25s}: {val}")

    print("\n4. RFG strong-field s_A — ცდის ღია ნაბიჯები")
    for i, task in enumerate(rfg_strong_field_open(), 1):
        print(f"  {i}. {task}")

    print("\n5. ფალსიფიკაციის ფანჯარა")
    fals = falsification_window()
    for key, val in fals.items():
        print(f"  {key:30s}: {val}")

    print("\n6. სტატუსი")
    print("  - PSR J1738 |α_0| < 2e-5 დაფიქსირებულია")
    print("  - RFG phase9 s_A = 0 postulated (minimal coupling) — consistent")
    print("  - Strong-field s_A derivation — ღია (Damour-Esposito-Farèse-ის ანალოგი)")
    print("  - If derived s_A > 1e-5 in NS — RFG ფალსიფიცირდება")
