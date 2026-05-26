# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.
# Active coefficient scheme: N/A; prediction ledger, no active c_Y sign claim.

"""
================================================================================
PHASE 36: RFG-სპეციფიკური უნიკალური ფალსიფიცირებადი პროგნოზები
================================================================================

სტატუსი:
Strategy 3 / X3-ის შესრულება. ეს ფაილი არჩევს მხოლოდ იმ ტესტებს, რომლებიც
შეიძლება RFG-ს გამოარჩევდეს GR/MOND/AQUAL/TeVeS/RMOND/AeST-გან.
ეს არის prediction ledger და claim gate, არა დასრულებული prediction paper.

მკაცრი კლასიფიკაცია:
    C1. N=4 lepton-like resonance: RFG-specific მხოლოდ მაშინ, თუ phase35
        N-ladder-ს ფესვიდან გამოიყვანს. ამ ეტაპზე პირობითია.
    C2. a_0(z): shared/phenomenological; MOND-family-სგან ბოლომდე არ გამოარჩევს.
    C3. photon ring: shared with many regular-BH models.
    C4. scalar breathing: shared with scalar-tensor theories.
    C5. framed vortex in dipolar supersolid: best near-term RFG candidate,
        but not proven unique until the framing observable is derived.

ამიტომ phase36-ის მკაცრი deliverable არის:
    1. C1: N=4 = 329 keV legacy conditional branch, not active main claim.
    2. C5: Dy/Er framed-vortex candidate protocol with concrete lab gates.
    3. Gate-aware migration ledger for old ISPG predictions.

განახლება phase37-44 შემდეგ:
    charged-lepton sector-ის მთავარი გზა გახდა C3 Koide operator და
    order-9 framed-holonomy candidate. p11-ის claim gate-ის მიხედვით
    cyclic lift, h=2 action derivation, m~nu^2 და radiative protection
    ჯერ ღიაა. ამიტომ C1/N=4 აღარ არის მთავარი particle-sector პროგნოზი;
    ის რჩება მხოლოდ legacy/conditional branch-ად, თუ ძველი radial
    N-ladder ოდესმე დამოუკიდებლად დამტკიცდება.
"""

import math


M_E_MEV = 0.51099895069
N_ELECTRON = 5
N4_PREDICTED_MEV = 0.329
N4_MODEL_BAND_MEV = 0.002
N4_THEORY_SIGMA_MEV = N4_MODEL_BAND_MEV
X17_MEV = 17.0


def n_ladder_mass(n_value):
    """Existing ladder relation m_N = m_e * (N/N_e)^2."""
    return M_E_MEV * (n_value / N_ELECTRON) ** 2


def n4_fixed_prediction():
    """
    Legacy N=4 fixed-number branch.

    q=0 ladder gives 327.04 keV. The theory text uses q-corrected 329 keV.
    Until phase35 derives the N-ladder and q-correction, 329 keV is a
    legacy model-band estimate, not an active particle prediction.
    """
    q0_mev = n_ladder_mass(4)
    center_mev = N4_PREDICTED_MEV
    band_mev = N4_MODEL_BAND_MEV
    return {
        "name": "C1_N4_lepton_like_resonance",
        "n": 4,
        "q0_mass_keV": q0_mev * 1000,
        "legacy_center_keV": center_mev * 1000,
        "model_band_keV": band_mev * 1000,
        "band_keV": ((center_mev - band_mev) * 1000, (center_mev + band_mev) * 1000),
        "active_prediction": False,
        "model_band_not_experimental_sigma": True,
        "rfg_uniqueness": "conditional: unique only if phase35 derives N-ladder",
        "particle_update": "legacy/suspended under C3 charged-lepton route",
        "operational_falsifier": False,
        "falsification": "blocked until coupling, decay channel, lifetime, spin/charge, and search topology are specified",
        "required_gates": [
            "derive radial N-ladder independently of lepton masses",
            "derive the q-correction from the RFG oscillon operator",
            "specify coupling, production rate, decay width, and lifetime",
            "map existing exclusion bounds for the chosen channel",
        ],
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
        "post_hoc_integer": nearest_integer,
        "operational_claim": False,
    }


def framed_vortex_protocol():
    """
    Best near-term RFG laboratory candidate.

    The strong claim is gated: Tw alone is not generally topologically
    frozen because Lk = Wr + Tw. RFG must define the physical framing
    vector, show why half-integer classes are protected, and distinguish
    the signal from conventional half-quantum/director vortices.
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
        "name": "C5_Dy_Er_framed_vortex",
        "legacy_name": "C5_Eu_Dy_framed_vortex",
        "observable": "candidate twist/framing class Tw of vortex core",
        "conditional_prediction": "Tw in half-integer classes if RFG framing vector is physically realized",
        "null_result": "validated measurement map finds continuous unprotected framing distribution",
        "positive_result": "robust half-integer class with RFG-specific discriminator survives systematics",
        "rfg_specific_status": "best candidate, not proven unique",
        "required_gates": [
            "derive physical framing vector from RFG medium variables",
            "prove half-integer class protection despite Wr/Tw exchange",
            "define image/phase reconstruction protocol for Tw",
            "quantify finite-temperature, trap, density, and reconnection systematics",
            "separate the signal from ordinary vortex quantization and half-quantum defects",
        ],
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
            "falsifiable_number": f"{c1['legacy_center_keV']:.0f} keV model band; not active prediction",
            "operational_ready": c1["operational_falsifier"],
        },
        {
            "candidate": "C5 framed vortex",
            "rfg_specific": c5["rfg_specific_status"],
            "time_window": "near-term cold atom labs after measurement protocol",
            "falsifiable_number": "conditional Tw = n/2 classes",
            "operational_ready": False,
        },
    ]


def status_assessment():
    return {
        "closed_now": "No observationally ready unique RFG prediction is closed in this file.",
        "c1_status": "legacy conditional after phase37-44; not a main prediction unless radial N-ladder is revived and couplings are specified.",
        "c5_status": "best near-term candidate, but still needs framing observable and non-RFG background rejection.",
        "dependency": "particle-sector main route is now phase37-44 C3/order-9 candidate, not N-ladder.",
    }


def external_observational_bounds():
    """
    Current external anchors that keep p12 from overclaiming.

    These entries are guardrails, not likelihood calculations.
    """
    return [
        {
            "sector": "X17",
            "anchor": "NA64 excludes part of the X17 electron-coupling region; MEG II found no significant excess in its searched channel.",
            "impact": "X17/N4 remains weak harmonic audit only.",
            "source": "NA64/CERN, MEG II/JINR reports",
        },
        {
            "sector": "dipolar supersolid vortices",
            "anchor": "Vortices have been observed in dipolar supersolids, making a vortex protocol experimentally plausible.",
            "impact": "C5 platform plausible; RFG-specific framing signal still not derived.",
            "source": "Nature 2024 Observation of vortices in a dipolar supersolid",
        },
        {
            "sector": "CMB/BAO",
            "anchor": "Planck 2018 remains a baseline CMB reference; DESI DR2 BAO data products are now available.",
            "impact": "p12 cannot claim Planck/BAO/LSS pass without a likelihood pipeline.",
            "source": "Planck 2018, DESI DR2 2025",
        },
        {
            "sector": "EHT shadow/ring",
            "anchor": "M87* crescent diameter is about 42 +/- 3 microarcsec; plasma, spin, inclination, and ray tracing are essential.",
            "impact": "static +4.63% shadow benchmark is not an EHT fit.",
            "source": "EHT M87* 2019 papers",
        },
        {
            "sector": "GW speed/polarization",
            "anchor": "GW170817/GRB170817A constrains c_g close to c at the 1e-15 level.",
            "impact": "scalar breathing and alpha_T claims require LVK posterior-level treatment.",
            "source": "GW170817 multimessenger bounds",
        },
    ]


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
    now stores them. Status words are gate-aware: this is not an empirical
    pass table.
    """
    return [
        {
            "sector": "weak field",
            "old_prediction": "gamma=beta=1; 1PN solar-system tests match GR",
            "rfg_status": "branch consistency benchmark",
            "claim_class": "branch_benchmark",
            "claim_gate": "p03 solar residual and preferred-frame gates",
            "current_file": "p03_solar.py",
            "formula": "gamma=1, beta=1; Mercury=42.98 arcsec/cy",
        },
        {
            "sector": "2PN light propagation",
            "old_prediction": "2PN Shapiro RFG-GR differential Delta_B=pi/4",
            "rfg_status": "branch consistency benchmark",
            "claim_class": "branch_benchmark",
            "claim_gate": "p03 Cassini/2PN observable gate",
            "current_file": "p03_solar.py",
            "formula": "Delta t=(r_g^2/(c*b))*pi/4",
        },
        {
            "sector": "2PN light propagation",
            "old_prediction": "2PN bending term enhanced by 16/15 over GR",
            "rfg_status": "branch consistency benchmark",
            "claim_class": "branch_benchmark",
            "claim_gate": "p03 high-precision astrometry gate",
            "current_file": "p03_solar.py",
            "formula": "theta_RFG=2r_s/b+pi*r_s^2/b^2",
        },
        {
            "sector": "compact objects",
            "old_prediction": "curvature invariants vanish at r->0; no finite-radius horizon",
            "rfg_status": "conditional compact-object construction",
            "claim_class": "conditional_ledger",
            "claim_gate": "p05 global solution, collapse, and stability gates",
            "current_file": "p05_compact.py",
            "formula": "R->0, K->0, Knudsen-core C2 matching completes endpoint",
        },
        {
            "sector": "compact objects",
            "old_prediction": "golden-ratio ISCO",
            "rfg_status": "conditional compact benchmark",
            "claim_class": "conditional_ledger",
            "claim_gate": "p05 ray-tracing and accretion-model gate",
            "current_file": "p05_compact.py",
            "formula": "r_ISCO=phi^2*r_s; f=0.931*f_GR",
        },
        {
            "sector": "compact objects",
            "old_prediction": "photon sphere r_s and shadow b_c=e*r_s",
            "rfg_status": "static benchmark only",
            "claim_class": "conditional_ledger",
            "claim_gate": "EHT plasma/spin/inclination/ray-tracing likelihood gate",
            "current_file": "p05_compact.py",
            "formula": "RFG shadow diameter is +4.63% vs GR",
        },
        {
            "sector": "gravitational waves",
            "old_prediction": "c_g=c exactly; breathing mode suppressed",
            "rfg_status": "conditional consistency with open polarization gate",
            "claim_class": "conditional_ledger",
            "claim_gate": "p04 LVK speed/polarization posterior gate",
            "current_file": "p04_gw.py",
            "formula": "alpha_T=0; A_b/A_t~r_s/r is a working estimate",
        },
        {
            "sector": "MOND/galaxies",
            "old_prediction": "a0=cH/(2*pi), mu=x/(1+x), BTFR, EFE",
            "rfg_status": "candidate data benchmark",
            "claim_class": "conditional_ledger",
            "claim_gate": "p07 SPARC/RAR/EFE likelihood gate",
            "current_file": "p07_mond.py",
            "formula": "g_h/g_N=a0/g -> mu=x/(1+x); v^4=G*M*a0",
        },
        {
            "sector": "MOND/formation memory",
            "old_prediction": "age/redshift-dependent a0 and BTFR residuals",
            "rfg_status": "seeded programme",
            "claim_class": "seeded",
            "claim_gate": "redshift galaxy sample and selection-function gate",
            "current_file": "p07_mond.py",
            "formula": "a0(z)=cH(z)/(2*pi); vortex memory/lag needs data closure",
        },
        {
            "sector": "cluster mergers",
            "old_prediction": "Bullet Cluster lensing peaks lock to galaxies via frozen hysteresis",
            "rfg_status": "candidate merger benchmark",
            "claim_class": "conditional_ledger",
            "claim_gate": "p09 lensing map, gas model, and merger simulation gate",
            "current_file": "p09_bullet.py",
            "formula": "tau_rel=c/g_vir ~ 680 Gyr >> tau_cross ~0.33 Gyr; peaks at x~±715 kpc",
        },
        {
            "sector": "CMB/cosmology",
            "old_prediction": "linear CMB sector matches LCDM in same-matter metric limit",
            "rfg_status": "same-matter compatibility branch",
            "claim_class": "branch_benchmark",
            "claim_gate": "p08 Boltzmann/Planck/BAO/BBN likelihood gate",
            "current_file": "p08_cmb.py",
            "formula": "Phi_0=X_0=0 -> alpha_K=alpha_B=alpha_M=alpha_T=0 -> C_l^RFG=C_l^LCDM in same-matter limit",
        },
        {
            "sector": "frame dragging",
            "old_prediction": "1.5PN Lense-Thirring matches GR; MOND rotation slot inert in Solar System",
            "rfg_status": "conditional preferred-frame benchmark",
            "claim_class": "conditional_ledger",
            "claim_gate": "p03 g_0i/preferred-frame gate",
            "current_file": "p03_solar.py",
            "formula": "Omega_LT standard if g_0i sector is GR-like",
        },
        {
            "sector": "quantum tests",
            "old_prediction": "gravitational dephasing, tunneling profile effects, birefringence",
            "rfg_status": "open",
            "claim_class": "open",
            "claim_gate": "p10/p11 finite-energy oscillon and EM/gauge completion gates",
            "current_file": "future particle/quantum phase, phase37 particle-sector extension",
            "formula": "Gamma_phi=2m|DeltaPhi_N|/h; other channels conditional",
        },
    ]


def migrated_prediction_scorecard():
    rows = old_to_rfg_prediction_map()
    counts = {
        "active_unique_ready": 0,
        "branch_benchmark": 0,
        "conditional_ledger": 0,
        "seeded": 0,
        "open": 0,
    }
    for row in rows:
        counts[row.get("claim_class", "open")] += 1

    return {
        "total_old_prediction_families": len(rows),
        "active_unique_ready": counts["active_unique_ready"],
        "branch_benchmark": counts["branch_benchmark"],
        "conditional_ledger": counts["conditional_ledger"],
        "seeded": counts["seeded"],
        "open_or_not_migrated": counts["open"],
        "verdict": "no old prediction family is counted as observationally closed here; p12 is a gate-aware ledger.",
    }


# =============================================================================
# STAGE D5: OLD quantum/lab prediction gate
# =============================================================================

def stage_d5_old_quantum_prediction_ledger():
    """
    Drain of the testable but lower-priority quantum predictions from
    OLD/ISPG_Quantum.tex into the prediction file.
    """
    return {
        "old_source": "OLD/ISPG_Quantum.tex",
        "new_file": "p12_predictions.py",
        "migrated": True,
        "quantum_dephasing": {
            "formula": "Gamma_phi = 2*m*|Delta Phi_N|/h",
            "best_regime": "massive mesoscopic interferometry, not optical photons",
            "status": "candidate lab prediction; requires noise/environment scan",
        },
        "fine_structure_variation": {
            "old_range": "Delta alpha/alpha ~ 1e-7--1e-6 at z~1--3",
            "status": "conditional; must be tied to a real EM/gauge-completion mechanism before publication",
        },
        "gravitational_birefringence": {
            "status": "very small, secondary prediction",
            "best_targets": "strong gradients around compact objects with high-frequency photons",
            "warning": "do not use as a main claim until the EM sector is derived",
        },
        "framed_vortex_lab": {
            "observable": "half-integer vortex framing/twist classes Tw=n/2 in dipolar supersolids",
            "status": "best near-term candidate; RFG-specificity and measurement map still open",
            "candidate_platforms": ["Dy/Er dipolar supersolids", "BEC vortex interferometry"],
        },
        "particle_legacy": {
            "N4_329keV": "legacy conditional branch, suspended under C3 charged-lepton route",
            "X17": "only a weak harmonic audit; no identification claimed",
        },
    }


def stage_d5_prediction_priority_order():
    """Publication/prioritization filter for old predictions."""
    return [
        {
            "priority": 1,
            "target": "gravity article",
            "use": "1PN match, 2PN Shapiro/bending, compact-object ISCO/shadow, GW/CMB consistency",
        },
        {
            "priority": 2,
            "target": "MOND/cluster articles",
            "use": "a0(z), SPARC, EFE, Bullet, residual cluster binding",
        },
        {
            "priority": 3,
            "target": "particle/quantum programme",
            "use": "C3/order-9 lepton-ratio candidate, EM/charge/spin/QCD/SM embedding, lab quantum tests",
        },
        {
            "priority": 4,
            "target": "low-priority speculative tests",
            "use": "alpha variation, birefringence, X17/N4 legacy branch",
        },
    ]


def p12_do_not_claim():
    """Prediction-sector overclaim blacklist."""
    return [
        "p12 contains an observationally confirmed unique RFG prediction.",
        "N=4 329 keV is an active main particle prediction.",
        "The 2 keV N4 band is an experimental sigma.",
        "X17 is explained or identified with the N4 branch.",
        "C5 framed vortex is already proven RFG-specific.",
        "Old ISPG predictions are recovered as observational passes.",
        "Planck/BAO/LSS/CMB likelihoods have been passed.",
        "EHT supports the static +4.63% shadow benchmark.",
        "GW scalar breathing is observationally allowed or excluded by this file.",
        "Bullet/SPARC/RAR are solved without the corresponding data likelihoods.",
        "Particle/quantum predictions are ready before p10/p11 gates close.",
    ]


def p12_prediction_claim_gate():
    """Single importable gate for the prediction file."""
    c1 = n4_fixed_prediction()
    c5 = framed_vortex_protocol()
    score = migrated_prediction_scorecard()
    return {
        "overall_status": "PREDICTION_LEDGER_NOT_PUBLICATION_READY",
        "active_unique_predictions_ready": False,
        "best_near_term_candidate": c5["name"],
        "c1_active_prediction": c1["active_prediction"],
        "c1_operational_falsifier": c1["operational_falsifier"],
        "c5_rfg_specific_status": c5["rfg_specific_status"],
        "migration_scorecard": score,
        "external_bounds": external_observational_bounds(),
        "hard_gates": [
            "C5: derive measurable RFG framing vector and distinguish from conventional half-quantum/director vortices",
            "C1: derive N-ladder, q-correction, coupling, width, lifetime, and exclusion map",
            "X17: derive transition/coupling rule or keep as weak harmonic audit",
            "Classical predictions: import each p03-p09 claim gate, not free recovered labels",
            "Cosmology: build Planck/BAO/BBN/LSS likelihood pipeline before empirical pass claims",
            "EHT/GW: use ray-tracing and LVK posterior-level comparisons before observational claims",
            "Particle/quantum: wait for p10/p11 finite-energy, m~nu^2, radiative, and gauge gates",
        ],
        "do_not_claim": p12_do_not_claim(),
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
            f"| ready={item['operational_ready']} "
            f"| number={item['falsifiable_number']} | window={item['time_window']}"
        )

    print("\n6. Status")
    for key, value in status_assessment().items():
        print(f"  {key:12s}: {value}")

    print("\n7. External observational guardrails")
    for item in external_observational_bounds():
        print(f"  {item['sector']:28s}: {item['impact']} ({item['source']})")

    print("\n8. Old ISPG constants ledger")
    for key, value in old_ispg_prediction_constants().items():
        print(f"  {key:38s}: {value:.8g}")

    print("\n9. Old-to-RFG prediction migration map")
    for item in old_to_rfg_prediction_map():
        print(
            f"  {item['sector']:24s} | {item['rfg_status']:24s} "
            f"| gate={item['claim_gate']} | {item['current_file']}"
        )

    print("\n10. Migration scorecard")
    for key, value in migrated_prediction_scorecard().items():
        print(f"  {key:28s}: {value}")

    print("\n11. Final p12 claim gate")
    gate = p12_prediction_claim_gate()
    print(f"  overall_status: {gate['overall_status']}")
    print(f"  active_unique_predictions_ready: {gate['active_unique_predictions_ready']}")
    print(f"  best_near_term_candidate: {gate['best_near_term_candidate']}")
    print("  hard gates:")
    for item in gate["hard_gates"]:
        print(f"    - {item}")
    print("  do-not-claim:")
    for item in gate["do_not_claim"]:
        print(f"    - {item}")
