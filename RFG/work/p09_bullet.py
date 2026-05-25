# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 20: Bullet Cluster - chi ველის მეხსიერების მექანიკა
================================================================================

სტატუსი:
ეს ფაილი აღარ არის მხოლოდ კონცეპტუალური ესკიზი. ძველი MOND/ISPG
Bullet-Cluster გამოთვლიდან გადმოტანილია frozen-hysteresis benchmark:
    freeze theorem  -> tau_rel=c/g_vir >> tau_cross
    localization    -> Helmholtz chi response peaks at galaxy centroids
    peak locking    -> lensing peaks stay with galaxies, not shocked gas
    threshold       -> O(0.3 M_sub) transported hysteretic mass is enough

ეს მაინც არ არის სრული N-body+gas+chi time-dependent simulation და არ აცხადებს
pixel-level Clowe et al. shear-map likelihood-ს. მაგრამ უკვე არის რაოდენობრივი
snapshot benchmark ძველი ცალკე Bullet calculation-ის ძირითადი რიცხვებით.

ცენტრალური ფაქტი (1E 0657-558, 2006):
    ორი გალაქტიკის გროვა შეჯახდა ~150 მილიონი წლის წინ.
    - გაზი (X-ray, ~90% ბარიონული მასა): შენელდა ram pressure-ით
    - გალაქტიკები (~10%): გაიარეს თავისუფლად (collisionless)
    - გრავიტაციული ლენზირება: პიკი გალაქტიკებთან,
      *არა* გაზთან, სადაც ბარიონული მასის უმეტესობაა

ეს არის *მთავარი წინააღმდეგობა MOND-ისთვის*:
    MOND-ის დამატებითი გრავიტაცია ეფუძნება ლოკალურ ბარიონულ მასას.
    გაზი (უმეტესობა) -> MOND მოელის ლენზირების პიკს გაზთან.
    დაკვირვება: პიკი გალაქტიკებთან -> MOND-ის სტანდარტული ფორმა ვერ ხსნის.

RFG-ის მექანიზმი:
    chi ველი ვორტექსულ-მეხსიერებითი ფაზაა (PLAN.md ფაზა 4).
    მისი თავისებურება: *მეხსიერების შენარჩუნება დიდი დროით*
    (tau_chi >> t_collision).

    chi-ის რეაგირების განტოლება:
        Box chi + (1/tau_chi) * (chi - chi_eq) = J

    სადაც chi_eq წყარო არის კინემატიკური/ვორტიკალური დინება (omega dot v).
    ეს არის სამუშაო ჰიპოთეზა: გაზი შოკის დროს თერმალიზდება და კარგავს ვორტიციტს,
    ხოლო გალაქტიკები ინარჩუნებენ კინემატიკურ სტრუქტურას.
    ძველი თეორიის გაძლიერებული ვერსია tau_chi-ს თავისუფალ პარამეტრად აღარ
    ტოვებს merger-freeze არგუმენტში. გამოიყენება უნივერსალური მაკრო-რელაქსაცია
        tau_rel = c / g_vir,
    რომელიც გალაქტიკების MOND/transport სექტორთან იგივე სკალირების ნაწილია.

დინამიკა:
    tau_chi ~ Gyr >> t_collision ~ 100 Myr -> chi ინარჩუნებს თავის
    პოზიციას (გალაქტიკებთან), არ მიჰყვება გაზს.

შესაბამისად:
    Phi_lens(x) = Phi_baryon(x) + Phi_chi(x)
    პიკი:     gas-ში (~90%) + გალაქტიკები-შიგნით (~10%) + chi (გალაქტიკებთან)
    შედეგი:   გალაქტიკები + chi >> გაზი -> პიკი გალაქტიკებთან

References:
    - Clowe et al. 2006 (ApJL 648:L109) - Bullet Cluster discovery
    - Paraficz et al. 2016 - 250 kpc aperture mass budget benchmark
    - OLD/ISPG_MOND.tex Sec. "The Bullet Cluster in the frozen-hysteresis limit"
    - Intuitive_Theory.md §6.2 - vortex memory mechanism
    - PLAN.md ფაზა 4 - chi field role
"""

import math

import sympy as sp
from sympy import Heaviside, Symbol, exp, oo, simplify, symbols


C_LIGHT = 299_792_458.0
G_NEWTON = 6.67430e-11
M_SUN = 1.98847e30
KPC = 3.0856775814913673e19
MPC = 1.0e3 * KPC
YEAR = 365.25 * 24.0 * 3600.0
GYR = 1.0e9 * YEAR
A0_OBS = 1.2e-10


# ==============================================================================
# Setup
# ==============================================================================

def setup_collision():
    """შეჯახების ბაზური პარამეტრები."""
    t, x = symbols("t x", real=True)
    v_coll = Symbol("v_coll", positive=True)  # შეჯახების სიჩქარე
    L_cluster = Symbol("L_cluster", positive=True)
    tau_chi = Symbol("tau_chi", positive=True)  # local response placeholder; merger freeze uses tau_rel=c/g
    return t, x, v_coll, L_cluster, tau_chi


# ==============================================================================
# ნაბიჯი 1: კომპონენტთა დინამიკა შეჯახებაში
# ==============================================================================

def step1_collision_dynamics():
    """
    შეჯახების სამი კომპონენტი:

    1. გაზი (hydrodynamic):
       v_gas(t) = v_0 * exp(-t/tau_drag)  <- ram pressure-ით შენელება
       tau_drag ~ 10 Myr (გროვის გადაკვეთის დროზე მცირე)
       -> გაზი ცენტრში გროვდება

    2. გალაქტიკები (collisionless):
       v_gal(t) = v_0  <- მუდმივი (ბალისტიკური მოძრაობა)
       -> გალაქტიკები გაივლიან თავიდან ბოლომდე

    3. chi ველი (memory):
       v_chi მიბმულია გალაქტიკების კინემატიკურ სტრუქტურაზე,
       რადგან tau_chi მნიშვნელოვნად დიდია.
       -> chi ჩამორჩება გაზს და მიჰყვება კინემატიკურ სტრუქტურას.
    """
    t, v_0, tau_drag = symbols("t v_0 tau_drag", positive=True)

    # გაზის სიჩქარე (exponential damping)
    v_gas = v_0 * exp(-t / tau_drag)
    x_gas = simplify(sp.integrate(v_gas, (t, 0, t)))

    # გალაქტიკების ბალისტიკური მოძრაობა
    v_gal = v_0
    x_gal = v_0 * t

    # გალაქტიკებსა და გაზს შორის დაშორება
    separation = simplify(x_gal - x_gas)

    return v_gas, x_gas, v_gal, x_gal, separation


# ==============================================================================
# ნაბიჯი 2: chi ველის რეაგირების განტოლება
# ==============================================================================

def step2_chi_field_equation():
    """
    chi ველის რეაგირების განტოლება (toy model):
        d chi / dt = -(1/tau_chi) * (chi - chi_eq(x,t))

    სადაც:
        chi_eq - წონასწორული მნიშვნელობა, რომელიც მოდის ვორტიკალური დინებიდან
                 (სამუშაო ჰიპოთეზა)
        tau_chi - რელაქსაციის/მეხსიერების დრო

    ინტეგრალური ფორმა Green-ის ფუნქციით:
        G(t,t') = exp(-(t - t')/tau_chi) * Heaviside(t - t') / tau_chi
        chi(t) = chi(0) * exp(-t/tau_chi) + integral_0^t G(t,t') * chi_eq(t') dt'

    აგენტთა საბჭოს დაზუსტება:
    vorticity source (omega dot v) არის ჰიპოთეზა, კოდში ჯერ არ არის გამოყვანილი.
    ასევე, ცვლადი წყაროს ინტეგრალური ფორმა აქ არ მოწმდება; კოდი იყენებს
    მხოლოდ constant-source toy solution-ს.
    """
    t, tau_chi = symbols("t tau_chi", positive=True)
    chi_0, omega_v, k_chi = symbols("chi_0 omega_v k_chi", real=True)
    t_p = Symbol("t_prime", positive=True)

    # Causal Green kernel
    green_kernel = exp(-(t - t_p) / tau_chi) * Heaviside(t - t_p) / tau_chi

    # ვორტიკალური წყარო კოდში მონიშნულია როგორც სიმბოლური მუდმივა.
    chi_eq_const = k_chi * omega_v

    # ინტეგრალური ამოხსნა მუდმივი წყაროსთვის (test case)
    chi_t = chi_0 * exp(-t / tau_chi) + chi_eq_const * (1 - exp(-t / tau_chi))
    chi_t = simplify(chi_t)

    # მეხსიერების რეჟიმი: t << tau_chi
    chi_short = simplify(chi_t.series(t, 0, 2).removeO())

    # რელაქსირებული რეჟიმი: t >> tau_chi
    chi_long = sp.limit(chi_t, t, oo)

    return chi_t, chi_short, chi_long, green_kernel


# ==============================================================================
# ნაბიჯი 3: Bullet Cluster-ის რიცხვობრივი მასშტაბები
# ==============================================================================

def step3_bullet_cluster_timescales():
    """
    ძველი თეორიის frozen-hysteresis დროითი ბირთვი.

    merger-freeze არგუმენტში რელაქსაციის დრო აღარ არის თავისუფალი tau_chi.
    ის იგივე transport-side relaxation scale-ია, რაც MOND სექტორში:

        tau_rel = c / g_vir,
        g_vir   = G M_sub / R_vir^2,
        tau_cross = R_vir / v_coll.

    fiducial values:
        M_sub ~ 1e14 M_sun, R_vir ~ 1 Mpc, v_coll ~ 3000 km/s
        g_vir ~ 1.4e-11 m/s^2
        tau_rel ~ 6.8e2 Gyr
        tau_cross ~ 3.3e-1 Gyr
        tau_rel/tau_cross ~ 2e3 >> 1

    Interpretation:
        chi cannot re-equilibrate onto shocked gas during the collision.
        It remains frozen into the collisionless galaxy component.
    """
    m_sub_msun = 1.0e14
    r_vir_mpc = 1.0
    v_coll_km_s = 3000.0
    d_collision_kpc = 300.0

    m_sub = m_sub_msun * M_SUN
    r_vir = r_vir_mpc * MPC
    v_coll = v_coll_km_s * 1000.0
    d_collision = d_collision_kpc * KPC

    g_vir = G_NEWTON * m_sub / r_vir**2
    tau_rel_gyr = (C_LIGHT / g_vir) / GYR
    tau_cross_gyr = (r_vir / v_coll) / GYR
    freeze_ratio = tau_rel_gyr / tau_cross_gyr
    a_collision = v_coll**2 / d_collision

    return {
        "M_sub_Msun": m_sub_msun,
        "R_vir_Mpc": r_vir_mpc,
        "v_coll_km_s": v_coll_km_s,
        "d_collision_kpc": d_collision_kpc,
        "g_vir_m_s2": g_vir,
        "tau_rel_Gyr": tau_rel_gyr,
        "tau_cross_Gyr": tau_cross_gyr,
        "tau_rel_over_tau_cross": freeze_ratio,
        "a_collision_m_s2": a_collision,
        "a_collision_over_a0": a_collision / A0_OBS,
        "verdict": "FROZEN: tau_rel/tau_cross >> 1, chi remains with collisionless galaxies",
    }


def freeze_theorem_symbolic():
    """
    Dimensionally closed freeze condition.

        tau_rel/tau_cross = (c/g_vir)/(R_vir/v_coll)
                            = c v_coll R_vir/(G M_sub).

    If this is much larger than 1, redistribution during the merger is
    mathematically suppressed without introducing a cluster-specific tau_chi.
    """
    c, v_coll, R_vir, G, M_sub = symbols("c v_coll R_vir G M_sub", positive=True)
    g_vir = G * M_sub / R_vir**2
    tau_rel = c / g_vir
    tau_cross = R_vir / v_coll
    ratio = simplify(tau_rel / tau_cross)

    return {
        "g_vir": g_vir,
        "tau_rel": tau_rel,
        "tau_cross": tau_cross,
        "freeze_ratio": ratio,
        "condition": "c*v_coll*R_vir/(G*M_sub) >> 1",
    }


# ==============================================================================
# ნაბიჯი 4: ლენზირების პოტენციალი - frozen-hysteresis benchmark
# ==============================================================================

def step4_frozen_hysteresis_lensing_benchmark():
    """
    ძველი ცალკე Bullet calculation-ის snapshot benchmark.

    Boundary conditions:
        galaxy centroids: ±720 kpc
        gas lags: main ~450 kpc, bullet/sub ~320 kpc
        aperture masses: M_main~2.5e14 Msun, M_sub~2.0e14 Msun
        baryonic masses: M_main,bary~2.03e14 Msun, M_sub,bary~1.51e14 Msun

    Dynamics:
        frozen positive Helmholtz chi response centered on collisionless galaxies
        L_chi ~ 100 kpc, f_hyst ~ 1.3

    Outputs copied from the old frozen-hysteresis benchmark:
        peaks at x ~ ±715 kpc, within ~5 kpc of galaxy centroids
        kappa_peak ~ 1.39 and 1.04, kappa_mid ~ 0.17
        main aperture fractions: gas~0.091, galaxy~0.121, chi~0.788
        sub aperture fractions:  gas~0.156, galaxy~0.112, chi~0.731
    """
    return {
        "boundary_conditions": {
            "galaxy_centroids_kpc": (-720.0, 720.0),
            "gas_lag_main_kpc": 450.0,
            "gas_lag_sub_kpc": 320.0,
            "M_main_total_250kpc_Msun": 2.5e14,
            "M_sub_total_250kpc_Msun": 2.0e14,
            "M_main_bary_Msun": 2.03e14,
            "M_sub_bary_Msun": 1.51e14,
        },
        "dynamics": {
            "kernel": "positive screened Helmholtz response of collisionless galaxy sources",
            "L_chi_kpc": 100.0,
            "f_hyst": 1.3,
            "M_chi_rule": "M_chi ≈ f_hyst * M_bary, order-unity transported hysteretic component",
        },
        "benchmark_outputs": {
            "kappa_peak_positions_kpc": (-715.0, 715.0),
            "peak_locking_error_kpc": 5.0,
            "kappa_peak_main": 1.39,
            "kappa_peak_sub": 1.04,
            "kappa_midplane": 0.17,
            "main_to_mid_contrast": 1.39 / 0.17,
            "sub_to_mid_contrast": 1.04 / 0.17,
        },
        "mass_budget": {
            "main": {"f_gas": 0.091, "f_gal": 0.121, "f_chi": 0.788},
            "sub": {"f_gas": 0.156, "f_gal": 0.112, "f_chi": 0.731},
            "observed_reference": "Paraficz et al. 2016: f_gas≈0.09±0.03, f_gal≈0.11±0.05, residual≈0.80±0.06",
        },
        "status": "quantitative frozen-hysteresis snapshot; not full pixel-level shear likelihood",
    }

def helmholtz_peak_locking_theorem():
    """
    Why the frozen chi component peaks at the galaxy position.

    In the frozen limit the chi map is a positive convolution of the
    collisionless-galaxy source with a screened Helmholtz kernel. A positive
    radially decreasing kernel centered at x_g has its maximum at x_g. The gas
    foreground can reduce contrast but cannot move the chi-only peak away from
    its collisionless source.
    """
    x, x_g, L_chi = symbols("x x_g L_chi", positive=True)
    kernel_1d_right = exp(-(x - x_g) / L_chi) / (2 * L_chi)
    derivative_right = simplify(sp.diff(kernel_1d_right, x))

    return {
        "kernel_right_branch": kernel_1d_right,
        "d_kernel_dx_right_of_source": derivative_right,
        "maximum": "at x=x_g for a positive screened kernel",
        "consequence": "chi-only convergence peaks are locked to galaxy centroids",
    }


def bullet_threshold_and_robustness():
    """
    Old benchmark robustness grid, compressed into data.

    Cells are (main-cluster peak shift in kpc, main-cluster f_chi).
    """
    robustness_grid = {
        0.8: {60: (5, 0.73), 80: (5, 0.72), 100: (5, 0.70), 120: (5, 0.67), 150: (5, 0.64)},
        1.0: {60: (0, 0.78), 80: (5, 0.76), 100: (5, 0.74), 120: (5, 0.72), 150: (5, 0.69)},
        1.3: {60: (0, 0.82), 80: (0, 0.80), 100: (5, 0.79), 120: (5, 0.77), 150: (5, 0.74)},
        1.6: {60: (0, 0.85), 80: (0, 0.84), 100: (5, 0.82), 120: (5, 0.81), 150: (5, 0.78)},
        2.0: {60: (0, 0.87), 80: (0, 0.86), 100: (0, 0.85), 120: (5, 0.84), 150: (5, 0.82)},
    }
    return {
        "existence_threshold": "M_chi/M_sub ≳ 0.22--0.41 for L_chi=80--150 kpc",
        "order_of_magnitude": "O(0.3 M_sub) transported hysteretic mass moves the dominant peak to the galaxy side",
        "fiducial": {"f_hyst": 1.3, "L_chi_kpc": 100, "main_peak_shift_kpc": 5, "main_f_chi": 0.79},
        "robustness_grid": robustness_grid,
        "verdict": "main peak locking is stable throughout the phenomenologically relevant grid",
    }


# ==============================================================================
# ნაბიჯი 5: MOND-ის კონტრასტი
# ==============================================================================

def step5_mond_contrast():
    """
    MOND (Milgrom) vs RFG-chi:
        MOND:  a_extra ∝ sqrt(g_N_baryonic) -> ბარიონების უმეტესობა გაზშია,
               პიკი უნდა იყოს გაზთან.
        RFG:   a_extra ∝ sqrt(g_chi_memory) -> პიკი რჩება გალაქტიკებთან
               chi ვორტექსული მეხსიერების გამო.
    """
    return {
        "MOND_failure": "a_extra ∝ sqrt(g_N_baryonic) -> peak at gas (ეწინააღმდეგება დაკვირვებას)",
        "RFG_success": "frozen chi/hysteresis source -> peak at galaxies in the snapshot benchmark",
        "mechanism": "tau_rel=c/g_vir >> tau_cross freezes the transported component onto collisionless galaxies",
        "no_new_particle": "chi is an effective collective memory field, not a sterile-neutrino/particle-DM species",
        "remaining_test": "full time-dependent N-body + gas hydro + chi evolution and pixel-level shear comparison",
    }


# ==============================================================================
# ნაბიჯი 6: ფალსიფიცირებადი პროგნოზები
# ==============================================================================

def step6_falsifiable_predictions():
    """
    RFG-ის ფალსიფიცირებადი პროგნოზები (ტექსტური ესკიზი):
    """
    return {
        "Bullet Cluster": "two convergence peaks locked to galaxy centroids; gas center remains subdominant",
        "Abell 520 (Train Wreck)": "peak between gas and galaxies (მეხსიერების რელაქსაციის ტესტი)",
        "MACS J0717.5 / El Gordo": "დამატებითი ტესტები tau_chi-ის შესაზღუდად",
        "Relaxation Law": "tau_rel=c/g_vir must work across cluster mergers without per-cluster tuning",
        "Screening Length": "L_chi should derive from intracluster thermodynamics, not from Bullet-only fitting",
        "Falsification": "pixel-level kappa maps fail if galaxy-locked chi peaks require cluster-by-cluster parameters",
    }


def status_upgrade_audit():
    return {
        "old_status": "conceptual sketch with free tau_chi",
        "new_status": "quantitative frozen-hysteresis snapshot benchmark",
        "closed_now": [
            "tau_rel=c/g_vir gives tau_rel/tau_cross ~ 2e3",
            "positive Helmholtz chi response locks peaks to galaxy centroids",
            "old benchmark peak positions, contrasts, mass budget, and threshold imported",
            "RFG avoids the standard MOND gas-peak failure in this frozen limit",
        ],
        "still_open": [
            "full time-dependent N-body+gas+chi simulation",
            "pixel-level weak-lensing likelihood against Clowe et al.",
            "first-principles derivation of L_chi for intracluster plasma conditions",
        ],
    }


def stage_b5_cluster_resonant_tail_binding_benchmark():
    """
    Drain of OLD/ISPG_MOND.tex cluster-binding section.

    Standard MOND/AQUAL gives much of the cluster enhancement but leaves a
    residual rich-cluster binding gap.  The old theory's useful idea is that
    dense intracluster emitters can retain an enhanced resonant-tail background,
    adding a local Bernoulli pressure deficit on top of the vortex/MOND channel.
    """
    m_bary_msun = 1.0e14
    r_cluster_mpc = 1.0
    m_bary = m_bary_msun * M_SUN
    r_cluster = r_cluster_mpc * MPC
    tau_esc_gyr = (r_cluster / C_LIGHT) / GYR
    tau_ret_50_gyr = 0.63
    v_vir = math.sqrt(G_NEWTON * m_bary / r_cluster)
    t_dyn_gyr = (r_cluster / v_vir) / GYR
    ell_mfp_50_kpc = (
        3.0 * r_cluster**2 / (C_LIGHT * tau_ret_50_gyr * GYR)
    ) / KPC

    retention_table = [
        {"tau_ret_over_tau_esc": 1, "tau_ret_Gyr": 0.003, "rho_extra_over_rho_cl": "0.003"},
        {"tau_ret_over_tau_esc": 100, "tau_ret_Gyr": 0.33, "rho_extra_over_rho_cl": "0.25--0.27"},
        {"tau_ret_over_tau_esc": 190, "tau_ret_Gyr": 0.63, "rho_extra_over_rho_cl": "0.50 threshold"},
        {"tau_ret_over_tau_esc": 300, "tau_ret_Gyr": 0.98, "rho_extra_over_rho_cl": "0.71--0.88"},
    ]

    return {
        "source": "OLD/ISPG_MOND.tex cluster-binding section",
        "problem": "simple MOND often leaves about a factor-of-two residual in rich clusters",
        "RFG_mechanism": (
            "overlapping irreversible resonant tails in dense ICM/galaxy environments "
            "raise the local background vibration and Bernoulli pressure deficit"
        ),
        "fiducial_cluster": {
            "M_bary_Msun": m_bary_msun,
            "R_Mpc": r_cluster_mpc,
            "delta_cl_0": 3.79e3,
            "delta_w0_high_amplitude_benchmark": 0.24,
            "tau_escape_Gyr": tau_esc_gyr,
        },
        "retention_table": retention_table,
        "tau_ret_50_Gyr": tau_ret_50_gyr,
        "t_dyn_Gyr": t_dyn_gyr,
        "tau_ret_50_over_t_dyn": tau_ret_50_gyr / t_dyn_gyr,
        "ell_mfp_50_kpc": ell_mfp_50_kpc,
        "interpretation": (
            "50% residual binding is reachable if the tail background is retained "
            "for about 0.4 dynamical times; ballistic escape would not close the gap."
        ),
    }


def stage_b5_cluster_retention_ode_symbolic():
    """Symbolic form of the old cluster retention benchmark."""
    H0, rho_de0, delta_w0, delta_cl, z, eps, tau_ret = symbols(
        "H0 rho_DE0 delta_w0 delta_cl z epsilon_cl tau_ret",
        positive=True,
    )
    phi_z, phi0 = symbols("phi_z phi0", real=True)
    source = (
        3 * H0 * rho_de0 * delta_w0 * delta_cl
        * (1 + z) ** 3
        * exp(2 * phi_z) / exp(2 * phi0)
    )
    return {
        "retention_ODE": sp.Eq(Symbol("epsilon_dot_cl"), source - eps / tau_ret),
        "normalization": "epsilon_dot(0)=3*H0*rho_DE0*delta_w0",
        "free_physical_variable": "tau_ret, the cluster residence time of the resonant-tail background",
        "not_a_new_particle": "the retained component is a collective resonant-tail background, not sterile-neutrino DM",
    }


def stage_b5_cluster_predictions_and_open_tasks():
    """Falsifiable outputs and the exact calculations still missing."""
    return {
        "predictions": [
            "residual binding should correlate with cold-front/sloshing activity at fixed baryonic mass",
            "extra binding should correlate with total baryonic density, gas plus galaxies",
            "rho_extra(r) should look like baryons convolved with a scattering/retention kernel",
            "clusters with weak ICM scattering should show a larger unresolved residual",
        ],
        "required_calculations": [
            "3D scalar solve on Chandra-derived cluster profiles",
            "wave-scattering calculation for ell_mfp from cold-front thickness and ICM power spectra",
            "weak-lensing comparison of the radial rho_extra/rho_bary profile",
            "joint Bullet/Abell520/ElGordo time-dependent N-body+gas+chi simulations",
        ],
        "status": (
            "Stage B5 imports the cluster-residual benchmark; it is a strong "
            "future paper target, not a completed cluster proof."
        ),
    }


# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 20: Bullet Cluster - chi ველის მეხსიერების მექანიკა")
    print("=" * 72)

    # ნაბიჯი 1
    print("\n--- ნაბიჯი 1: შეჯახების სამი კომპონენტი ---")
    v_gas, x_gas, v_gal, x_gal, sep = step1_collision_dynamics()
    print(f"  გაზი: v(t) = {v_gas}")
    print(f"  გალაქტიკები: v = {v_gal}, x = {x_gal}")
    print("  -> შედეგი: გაზი ნელდება, გალაქტიკები წინ მიდიან")

    # ნაბიჯი 2
    print("\n--- ნაბიჯი 2: chi ველის რეაგირების განტოლება ---")
    chi_t, chi_short, chi_long, green_k = step2_chi_field_equation()
    print(f"  Green Kernel G(t,t') = {green_k}")
    print(f"  chi(t) = {chi_t}")
    print(f"  t << tau_chi-ში: chi(t) ≈ {chi_short}")
    print("  -> მეხსიერების რეჟიმში chi ინარჩუნებს თავდაპირველ კინემატიკურ კავშირს")

    # ნაბიჯი 3
    print("\n--- ნაბიჯი 3: frozen-hysteresis დროითი იერარქია ---")
    times = step3_bullet_cluster_timescales()
    for k, v in times.items():
        if isinstance(v, float):
            print(f"  {k:28s}: {v:.6g}")
        else:
            print(f"  {k:28s}: {v}")

    print("\n--- ნაბიჯი 3b: freeze theorem სიმბოლურად ---")
    freeze = freeze_theorem_symbolic()
    for k, v in freeze.items():
        print(f"  {k:18s}: {v}")

    # ნაბიჯი 4
    print("\n--- ნაბიჯი 4: frozen-hysteresis lensing snapshot benchmark ---")
    benchmark = step4_frozen_hysteresis_lensing_benchmark()
    for section, values in benchmark.items():
        print(f"  {section}: {values}")

    print("\n--- ნაბიჯი 4b: Helmholtz peak-locking theorem ---")
    locking = helmholtz_peak_locking_theorem()
    for k, v in locking.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 4c: threshold და robustness grid ---")
    robust = bullet_threshold_and_robustness()
    for k, v in robust.items():
        print(f"  {k:24s}: {v}")

    # ნაბიჯი 5
    print("\n--- ნაბიჯი 5: MOND-ის კლასიკური მოლოდინი vs RFG-ის chi-მოდელი ---")
    contrast = step5_mond_contrast()
    for k, v in contrast.items():
        print(f"  {k:18s}: {v}")

    # ნაბიჯი 6
    print("\n--- ნაბიჯი 6: ფალსიფიცირებადი პროგნოზები ---")
    pred = step6_falsifiable_predictions()
    for k, v in pred.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 7: სტატუსის upgrade audit ---")
    audit = status_upgrade_audit()
    for k, v in audit.items():
        print(f"  {k:18s}: {v}")

    print("\n--- ნაბიჯი 8: STAGE B5 cluster residual-binding benchmark ---")
    cluster = stage_b5_cluster_resonant_tail_binding_benchmark()
    for k, v in cluster.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 8b: cluster retention ODE ---")
    retention = stage_b5_cluster_retention_ode_symbolic()
    for k, v in retention.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 8c: cluster predictions/open tasks ---")
    cluster_tasks = stage_b5_cluster_predictions_and_open_tasks()
    for k, v in cluster_tasks.items():
        print(f"  {k:30s}: {v}")

    print("\n" + "=" * 72)
    print("სტატუსი და აგენტთა საბჭოს შენიშვნების დადასტურება:")
    print("=" * 72)
    print("""
    1. სტატუსი: quantitative frozen-hysteresis snapshot benchmark.
       სრული validation კვლავ საჭიროებს κ(x,y), shear map, gas profile,
       galaxy distribution, projected mass და relativistic potentials შედარებას.
    2. tau_chi თავისუფალი merger პარამეტრი აღარ არის მთავარი საყრდენი:
       freeze hierarchy მოდის tau_rel=c/g_vir კანონიდან.
    3. Green-ის ფუნქციის ნიშანი: e^(-(t-t')/tau_chi) კოდში დაფიქსირდა,
       თუმცა ცვლადი წყაროს ინტეგრალური ფორმა აქ არ მოწმდება; გამოიყენება
       constant-source toy solution.
    4. chi ველის წყარო: კოდში დაემატა k_chi * omega_v (ვორტიკალური დინება)
       როგორც სიმბოლური ჰიპოთეზა/toy-წყარო.
    5. ძველი ცალკე Bullet calculation-ის ძირითადი snapshot outputs გადმოტანილია:
       peak positions, kappa contrast, mass-budget fractions, threshold და robustness grid.
    """)
    print("=" * 72)
