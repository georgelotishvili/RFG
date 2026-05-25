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


def old_ispg_prediction_constants():
    """Numerical constants recovered from the old ISPG prediction package."""
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    isco_rfg_over_rs = phi**2
    isco_gr_iso_over_rs = (5.0 + 2.0 * math.sqrt(6.0)) / 4.0
    shadow_ratio = 2.0 * math.e / (3.0 * math.sqrt(3.0))
    shapiro_delta_b = math.pi / 4.0
    bending_ratio = 16.0 / 15.0

    # Phase18 expression: Omega_RFG^2/Omega_GR^2 at their respective ISCO radii.
    omega_rfg_sq_over_gr_sq = (
        108.0
        * math.exp(-4.0 / (math.sqrt(5.0) + 3.0))
        / (29.0 + 13.0 * math.sqrt(5.0))
    )
    frequency_ratio = math.sqrt(omega_rfg_sq_over_gr_sq)

    return {
        "golden_ratio_phi": phi,
        "r_ISCO_over_r_s": isco_rfg_over_rs,
        "GR_isotropic_ISCO_over_r_s": isco_gr_iso_over_rs,
        "ISCO_radius_ratio_RFG_over_GR_iso": isco_rfg_over_rs / isco_gr_iso_over_rs,
        "ISCO_frequency_ratio_RFG_over_GR": frequency_ratio,
        "shadow_ratio_RFG_over_GR": shadow_ratio,
        "shadow_shift_percent": (shadow_ratio - 1.0) * 100.0,
        "Shapiro_Delta_B": shapiro_delta_b,
        "bending_2PN_ratio_RFG_over_GR": bending_ratio,
        "bending_2PN_enhancement_percent": (bending_ratio - 1.0) * 100.0,
    }


def old_to_rfg_prediction_map():
    """
    Inventory of old-theory predictions and where the current RFG work
    now reproduces them or keeps them conditional.
    """
    return [
        {
            "sector": "weak field",
            "old_prediction": "gamma=beta=1; 1PN solar-system tests match GR",
            "rfg_status": "recovered",
            "current_file": "phase8_weak_field.py, phase12_perihelion.py",
            "formula": "gamma=1, beta=1; Mercury=42.98 arcsec/cy",
        },
        {
            "sector": "2PN light propagation",
            "old_prediction": "2PN Shapiro RFG-GR differential Delta_B=pi/4",
            "rfg_status": "recovered",
            "current_file": "phase14_shapiro_2pn.py",
            "formula": "Delta t=(r_g^2/(c*b))*pi/4",
        },
        {
            "sector": "2PN light propagation",
            "old_prediction": "2PN bending term enhanced by 16/15 over GR",
            "rfg_status": "recovered",
            "current_file": "phase14_shapiro_2pn.py",
            "formula": "theta_RFG=2r_s/b+pi*r_s^2/b^2",
        },
        {
            "sector": "compact objects",
            "old_prediction": "curvature invariants vanish at r->0; no finite-radius horizon",
            "rfg_status": "strengthened",
            "current_file": "phase18_bh_singularity.py",
            "formula": "R->0, K->0, Knudsen-core C2 matching completes endpoint",
        },
        {
            "sector": "compact objects",
            "old_prediction": "golden-ratio ISCO",
            "rfg_status": "recovered",
            "current_file": "phase18_bh_singularity.py",
            "formula": "r_ISCO=phi^2*r_s; f=0.931*f_GR",
        },
        {
            "sector": "compact objects",
            "old_prediction": "photon sphere r_s and shadow b_c=e*r_s",
            "rfg_status": "recovered",
            "current_file": "phase18_bh_singularity.py, phase29_eht_shadow.py",
            "formula": "RFG shadow diameter is +4.63% vs GR",
        },
        {
            "sector": "gravitational waves",
            "old_prediction": "c_g=c exactly; breathing mode suppressed",
            "rfg_status": "partly recovered",
            "current_file": "phase9_gravitational_waves.py",
            "formula": "alpha_T=0; A_b/A_t~r_s/r is a working estimate",
        },
        {
            "sector": "MOND/galaxies",
            "old_prediction": "a0=cH/(2*pi), mu=x/(1+x), BTFR, EFE",
            "rfg_status": "recovered and strengthened",
            "current_file": "phase33_sparc_rotation_curves.py",
            "formula": "g_h/g_N=a0/g -> mu=x/(1+x); v^4=G*M*a0",
        },
        {
            "sector": "MOND/formation memory",
            "old_prediction": "age/redshift-dependent a0 and BTFR residuals",
            "rfg_status": "seeded",
            "current_file": "phase33_sparc_rotation_curves.py",
            "formula": "a0(z)=cH(z)/(2*pi); vortex memory/lag needs data closure",
        },
        {
            "sector": "cluster mergers",
            "old_prediction": "Bullet Cluster lensing peaks lock to galaxies via frozen hysteresis",
            "rfg_status": "recovered and strengthened",
            "current_file": "phase20_bullet_cluster.py",
            "formula": "tau_rel=c/g_vir ~ 680 Gyr >> tau_cross ~0.33 Gyr; peaks at x~±715 kpc",
        },
        {
            "sector": "CMB/cosmology",
            "old_prediction": "linear CMB sector matches LCDM in same-matter metric limit",
            "rfg_status": "recovered and analytically closed",
            "current_file": "phase21_cmb.py",
            "formula": "Phi_0=X_0=0 -> alpha_K=alpha_B=alpha_M=alpha_T=0 -> C_l^RFG=C_l^LCDM in same-matter limit",
        },
        {
            "sector": "frame dragging",
            "old_prediction": "1.5PN Lense-Thirring matches GR; MOND rotation slot inert in Solar System",
            "rfg_status": "partly recovered; preferred-frame tightening",
            "current_file": "phase30_lense_thirring.py",
            "formula": "Omega_LT standard if g_0i sector is GR-like",
        },
        {
            "sector": "quantum tests",
            "old_prediction": "gravitational dephasing, tunneling profile effects, birefringence",
            "rfg_status": "not yet migrated",
            "current_file": "future particle/quantum phase, likely phase47 extension",
            "formula": "Gamma_phi=2m|DeltaPhi_N|/h; other channels conditional",
        },
    ]


def migrated_prediction_scorecard():
    rows = old_to_rfg_prediction_map()
    counts = {"recovered": 0, "strengthened": 0, "partly": 0, "seeded": 0, "open": 0}
    for row in rows:
        status = row["rfg_status"]
        if status.startswith("recovered"):
            counts["recovered"] += 1
        elif status.startswith("strengthened"):
            counts["strengthened"] += 1
        elif status.startswith("partly"):
            counts["partly"] += 1
        elif status.startswith("seeded"):
            counts["seeded"] += 1
        else:
            counts["open"] += 1

    return {
        "total_old_prediction_families": len(rows),
        "recovered": counts["recovered"],
        "strengthened": counts["strengthened"],
        "partly_recovered": counts["partly"],
        "seeded": counts["seeded"],
        "open_or_not_migrated": counts["open"],
        "verdict": "the main classical old-theory predictions now reappear in RFG; quantum-test channels remain the least migrated.",
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

    print("\n7. Old ISPG constants recovered in RFG")
    for key, value in old_ispg_prediction_constants().items():
        print(f"  {key:38s}: {value:.8g}")

    print("\n8. Old-to-RFG prediction migration map")
    for item in old_to_rfg_prediction_map():
        print(
            f"  {item['sector']:24s} | {item['rfg_status']:24s} "
            f"| {item['formula']} | {item['current_file']}"
        )

    print("\n9. Migration scorecard")
    for key, value in migrated_prediction_scorecard().items():
        print(f"  {key:28s}: {value}")
