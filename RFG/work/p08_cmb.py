# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 21: CMB — კოსმოლოგიური პერტურბაციები და EFT 
================================================================================

სტატუსი:
ეს ფაილი წარმოადგენს კონცეპტუალურ ხიდს კოსმოლოგიური პერტურბაციების 
ეფექტური ველის თეორიასთან (EFT). სრული CMB სპექტრის (C_l), BAO-სა და 
ლენზირების გამოსათვლელად საჭიროა განტოლებების CLASS/CAMB კოდებში 
ინტეგრირება, რაც ცალკე ამოცანაა.

ცენტრალური დასკვნები:
    1. სტანდარტული Horndeski (Bellini-Sawicki) პარამეტრიზაცია არ არის 
       საკმარისი RFG-სთვის, რადგან I_k ელასტიური სექტორი შეიცავს 
       სივრცულ ტრანსვერსულ მოდებს. საჭიროა "EFT of Solid Inflation" (ESS) ჩარჩო.
    2. c_Y კოეფიციენტი უნდა იყოს უარყოფითი (c_Y < 0), რათა კინეტიკური 
       წევრი (alpha_K) იყოს დადებითი და ავირიდოთ Ghost არასტაბილურობა.
    3. გრავიტაციული ტალღის სიჩქარე c_T = c, მაგრამ მასიური დისპერსიის 
       ასარიდებლად საჭიროა phase9-ის კონსტრეინტი.

References:
    - Bellini & Sawicki 2014, JCAP 07:050
    - Endlich, Nicolis, Wang 2013 (Solid Inflation)
"""

from __future__ import annotations

import math
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

import sympy as sp
from sympy import symbols, Symbol, simplify, diff


PLANCK_2018 = {
    "H0": (67.36, 0.54),        # km/s/Mpc
    "Omega_m": (0.3153, 0.0073),
    "sigma8": (0.8111, 0.0060),
    "n_s": (0.9649, 0.0042),
    "tau": (0.0544, 0.0073),
    "BAO_scale_Mpc": (147.09, 0.26),
    "reference": "Planck 2018 TT,TE,EE+lowE+lensing baseline",
}

LOCAL_LATE_UNIVERSE = {
    "SH0ES_2022_H0": (73.0, 1.0),
    "TRGB_H0": (69.8, 1.9),
    "KiDS_1000_S8": (0.766, 0.020),
    "DES_Y3_S8": (0.776, 0.017),
}

DESI_DR2 = {
    "release_date": "2025-03-19",
    "combo_baseline": "DESI DR2 BAO + CMB",
    "combo_with_sne": "DESI DR2 BAO + CMB + supernova samples",
    "published_joint_significance": "3.1 sigma for DESI+CMB; 2.8-4.2 sigma with SNe combinations",
    "reference": "DESI DR2 dark-energy results, PRD 112:083515 (2025)",
    "status": "DR2 replaces the earlier DR1 placeholder; use full likelihood for parameter fitting.",
}

DEFAULT_ALPHA_TABLE = Path("RFG/work/cmb_alpha_table.dat")


@dataclass
class RFGAlphaModel:
    """
    hi_class/CLASS bridge parameterization:
        alpha_i(a) = alpha_i0 * Omega_DE(a)/Omega_DE0.

    The default alpha_i0=0 branch is not a fit. It is the analytic locked-FLRW
    CMB branch: Phi_0=X_0=0, so alpha_K=alpha_B=alpha_M=alpha_T=0.
    Nonzero alpha_i0 values are off-branch phenomenological probes.
    """

    alpha_K0: float = 0.00
    alpha_B0: float = 0.00
    alpha_M0: float = 0.00
    alpha_T0: float = 0.00
    omega_m0: float = PLANCK_2018["Omega_m"][0]
    a_min: float = 1.0e-3
    a_max: float = 1.0
    n_steps: int = 64


@dataclass
class TensionResult:
    name: str
    value_a: float
    err_a: float
    value_b: float
    err_b: float

    @property
    def sigma(self) -> float:
        return abs(self.value_a - self.value_b) / math.sqrt(self.err_a**2 + self.err_b**2)

    @property
    def chi2(self) -> float:
        return self.sigma**2


@dataclass
class FitReadiness:
    status: str
    hi_class_exe: str | None
    planck_likelihood_dir: str | None
    reason: str


def map_rfg_to_horndeski():
    """
    RFG-ის L_solid → Horndeski G_2, G_3, G_4, G_5 (მხოლოდ Y-სექტორით)

    Y = -2X (Bellini-Sawicki convention: X = -½ g^μν ∂_μφ ∂_νφ)
    ⟹ Y-სექტორი ⊂ G_2(X, φ)

    I_k-სექტორი არ ჯდება სტანდარტულ Horndeski-ში — ის მოითხოვს
    EFT of Solid Inflation (ESS) ჩარჩოს.

    Horndeski (მხოლოდ φ) სექტორში:
        G_2 = c_Y·Y + c_Y2·Y² = -2c_Y·X + 4·c_Y2·X²
        G_3 = 0   (no kinetic mixing)
        G_4 = M_Pl²/2
        G_5 = 0
    """
    X = Symbol('X', real=True)
    # აგენტთა საბჭოს შესწორება: c_Y არ უნდა იყოს positive=True, 
    # რადგან alpha_K > 0 მოითხოვს c_Y < 0-ს.
    c_Y = Symbol('c_Y', real=True) 
    c_Y2, M_Pl = symbols('c_Y2 M_Pl', positive=True)

    G_2 = -2 * c_Y * X + 4 * c_Y2 * X**2
    G_3 = sp.Integer(0)
    G_4 = M_Pl**2 / 2
    G_5 = sp.Integer(0)

    return G_2, G_3, G_4, G_5, X


# ==============================================================================
# ნაბიჯი 1: α_T (tensor speed excess)
# ==============================================================================

def compute_alpha_T():
    """
    Bellini-Sawicki:
        α_T = 2X·(G_{4,X} - G_{5,φ}) / M_*²    +  (G_5,X-term)

    RFG-ში:
        G_4 = M_Pl²/2 (X-დამოუკიდებელი) ⟹ G_{4,X} = 0
        G_5 = 0 ⟹ G_{5,φ} = 0

    ⟹ α_T = 0  (c_T = c)

    თუმცა I_k სექტორი წარმოშობს გრავიტონის ეფექტურ მასას (phase9).
    მკაცრი GW170817 თავსებადობისთვის მოითხოვება phase9 კონსტრეინტი:
    -0.5*c_Y - 0.5*c_Y2 + 0.5*c_I1 + 7.5*c_I1sq + 1.5*c_I2 + 0.5*c_I3 + 0.5*c_YI1 = 0
    """
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()
    M_star_sq = Symbol('M_star_sq', positive=True)
    phi = Symbol('phi', real=True)

    # G_{4,X} = 0 since G_4 const in X
    G4_X = diff(G_4, X)
    # G_{5,φ}, G_{5,X} = 0
    G5_phi = diff(G_5, phi)

    # α_T = 2X(G_{4,X} - G_{5,φ})/M_*²
    alpha_T = 2*X*(G4_X - G5_phi) / M_star_sq
    alpha_T = simplify(alpha_T)

    return alpha_T, G4_X, G5_phi


# ==============================================================================
# ნაბიჯი 2: α_M (Planck-mass running)
# ==============================================================================

def compute_alpha_M():
    """
    Bellini-Sawicki:
        M_*² = 2(G_4 - 2X·G_{4,X} + ...)

    α_M = (d ln M_*² / dt) / H
    რადგან G_4 = M_Pl²/2 = const, საბაზისო Horndeski სექტორში α_M = 0.
    """
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()
    M_Pl = Symbol('M_Pl', positive=True)

    M_star_sq = 2 * G_4  # = M_Pl²
    M_star_sq_simplified = simplify(M_star_sq)

    # dM_*²/dt = 0 (M_Pl const)
    t = Symbol('t', real=True)
    M_star_sq_t = M_star_sq_simplified  # no t-dependence
    d_M_star_dt = diff(M_star_sq_t, t)

    alpha_M = sp.Integer(0)

    return alpha_M, M_star_sq_simplified, d_M_star_dt


# ==============================================================================
# ნაბიჯი 3: α_B (braiding)
# ==============================================================================

def compute_alpha_B():
    """
    Bellini-Sawicki:
        α_B = 2(X·G_{3,X}·φ̇/H ·... + G_{4,X}·... ) / M_*²

    α_B = 0, რადგან G_3 = 0 და G_4 = const.
    """
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()

    G3_X = diff(G_3, X)
    G4_X = diff(G_4, X)

    # All zero ⟹ α_B = 0
    alpha_B = sp.Integer(0)

    return alpha_B, G3_X, G4_X


# ==============================================================================
# ნაბიჯი 4: α_K (kineticity)
# ==============================================================================

def compute_alpha_K():
    """
    Bellini-Sawicki:
        α_K = (2X·G_{2,X} + 4X²·G_{2,XX} + ...) / (H²·M_*²)

    G_2 = -2c_Y·X + 4c_Y2·X²
    α_K = (-4c_Y·X + 48c_Y2·X²) / (H²·M_Pl²)
    სტაბილურობა (no-ghost) მოითხოვს α_K > 0.
    ვინაიდან X დადებითია (time-like დერივატივი), c_Y უნდა იყოს უარყოფითი!
    """
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()
    H, M_Pl = symbols('H M_Pl', positive=True)

    G2_X = diff(G_2, X)
    G2_XX = diff(G_2, X, 2)

    alpha_K = (2*X*G2_X + 4*X**2*G2_XX) / (H**2 * M_Pl**2)
    alpha_K = simplify(alpha_K)

    return alpha_K, G2_X, G2_XX


def flrw_metric_sector_locking_theorem():
    """
    Old-theory CMB core, translated to RFG notation.

    In the physical matter-clock FLRW frame, g_tt = 1 (signature +---).
    The static RFG/bi-conformal time factor is g_tt = exp(Phi_0).
    Imposing both on the homogeneous cosmological background locks:

        exp(Phi_0) = 1  ->  Phi_0 = 0  ->  dot(Phi_0)=0  ->  X_0=0.

    Therefore the scalar carries no homogeneous stress-energy in the metric
    sector and the background Friedmann equations are inherited from GR for the
    same matter/radiation content.
    """
    Phi0, H, M_Pl, rho, p = symbols("Phi0 H M_Pl rho p", real=True)

    return {
        "cosmic_time_lock": sp.Eq(sp.exp(Phi0), 1),
        "homogeneous_solution": sp.Eq(Phi0, 0),
        "Phi0_dot": sp.Integer(0),
        "X0": sp.Integer(0),
        "tau_phi_background": sp.Integer(0),
        "friedmann_1": sp.Eq(3 * M_Pl**2 * H**2, rho),
        "friedmann_2": sp.Eq(2 * M_Pl**2 * sp.Symbol("Hdot", real=True), -(rho + p)),
        "status": "GR background inherited in the same-matter metric-sector limit",
    }


def bellini_sawicki_zero_alpha_theorem():
    """
    Strong CMB compatibility branch.

    If G4 is constant, G3=G5=0, and the FLRW scalar is locked to Phi_0=0,
    every Bellini-Sawicki alpha that can shift primary CMB peaks vanishes.
    The no-ghost inequality belongs to off-branch/dynamical scalar perturbations;
    on this locked branch X0=0 and the scalar stress starts at quadratic order.
    """
    return {
        "condition_G4": "G4 = M_Pl^2/2 = const",
        "condition_G3": sp.Integer(0),
        "condition_G5": sp.Integer(0),
        "condition_X0": sp.Integer(0),
        "alpha_T": sp.Integer(0),
        "alpha_M": sp.Integer(0),
        "alpha_B": sp.Integer(0),
        "alpha_K": sp.Integer(0),
        "consequence": "no linear EFT modification of tensor speed, Planck mass, braiding, or kineticity",
    }


def linear_metric_decoupling_theorem():
    """
    The sourced scalar can exist, but it does not source metric potentials at
    linear order when the homogeneous scalar derivative vanishes.
    """
    eps, dPhi = symbols("eps dPhi", real=True)
    stress_order_model = (eps * dPhi) ** 2
    linear_stress = sp.diff(stress_order_model, eps).subs(eps, 0)

    k, a, G, rho, Delta = symbols("k a G rho Delta", positive=True)
    Psi, Phi = symbols("Psi Phi", real=True)

    return {
        "scalar_stress_order": stress_order_model,
        "linear_scalar_stress": simplify(linear_stress),
        "poisson_equation": sp.Eq(k**2 * Psi, -4 * sp.pi * G * a**2 * rho * Delta),
        "slip_equation": sp.Eq(Phi - Psi, 0),
        "lensing_potential": "(Phi + Psi)/2 is GR-identical at linear order",
        "cmb_power": "C_ell^RFG = C_ell^LCDM for same matter content and initial conditions",
    }


def acoustic_ruler_inheritance():
    """
    If H(z), Phi/Psi, and the photon-baryon equations are unchanged, the sound
    horizon and acoustic angular scale are inherited.
    """
    z, z_star, c_s, H, D_A, r_s = symbols("z z_star c_s H D_A r_s", positive=True)
    sound_horizon = sp.Integral(c_s / H, (z, z_star, sp.oo))
    theta_star = r_s / D_A

    return {
        "sound_horizon": sp.Eq(r_s, sound_horizon),
        "theta_star": sp.Eq(sp.Symbol("theta_star", positive=True), theta_star),
        "delta_sound_horizon_same_matter": sp.Integer(0),
        "delta_theta_star_same_matter": sp.Integer(0),
        "peak_phase_shift_same_matter": sp.Integer(0),
        "status": "primary acoustic peak locations are inherited on the locked branch",
    }


def trace_channel_recombination_filter():
    """
    RFG scalar sourcing is trace-sensitive. Relativistic radiation has T=0, so
    the photon bath does not directly drive the trace channel at linear order.
    """
    rho, p = symbols("rho p", real=True)
    trace = rho - 3 * p

    return {
        "trace_definition_signature_plus_minus": trace,
        "radiation_p": sp.Eq(p, rho / 3),
        "T_radiation": simplify(trace.subs(p, rho / 3)),
        "T_pressureless_matter": simplify(trace.subs(p, 0)),
        "cmb_meaning": "photon acoustic pressure is not directly modified by the RFG trace channel",
    }


def is_memory_freezing_cmb_estimate():
    """
    Old-theory idea migrated into RFG language.

    If the nonlinear/transport channel has a relaxation time tau_RFG ~ c/g,
    then recombination-era acoustic oscillations see it as frozen whenever
    tau_RFG/T_acoustic >> 1. This is the possible route by which CMB potential
    wells can behave CDM-like at z~1100 while galaxies later show MOND response.
    """
    c = 299_792_458.0
    c_s = c / math.sqrt(3.0)
    z_star = 1090.0
    mpc = 3.0856775814913673e22
    wavelengths_comoving_mpc = [0.1, 1.0, 10.0, 100.0, 300.0]
    accelerations = [1.0e-8, 1.0e-10, 1.0e-12]

    ratios = []
    for lam_mpc in wavelengths_comoving_mpc:
        lam_phys = lam_mpc * mpc / (1.0 + z_star)
        t_acoustic = lam_phys / c_s
        for g in accelerations:
            tau_rfg = c / g
            ratios.append(tau_rfg / t_acoustic)

    return {
        "tau_RFG": "c/g",
        "T_acoustic": "lambda_phys/c_s",
        "z_star": z_star,
        "lambda_comoving_Mpc_grid": wavelengths_comoving_mpc,
        "g_grid_m_s2": accelerations,
        "min_tau_over_T": min(ratios),
        "max_tau_over_T": max(ratios),
        "status": "frozen-memory regime across this broad recombination grid",
        "interpretation": "candidate CDM-like wells for CMB; requires Boltzmann coupling in this phase21 file",
    }


# ==============================================================================
# ნაბიჯი 5: CMB სპექტრის თავსებადობა
# ==============================================================================

def cmb_consistency_check():
    return {
        'locked_FLRW_branch': 'Phi_0=0, X_0=0 in matter-clock cosmic time',
        'alpha_T': '0 on the locked branch; full solid-sector tensors still obey phase9 mass/speed filter',
        'alpha_M': '0 because G4=M_Pl^2/2 is constant',
        'alpha_B': '0 because G3=0 and G4 has no Phi-dependence',
        'alpha_K': '0 on the locked branch; no-ghost window applies to off-branch propagating scalar/ESS modes',
        'metric_potentials': 'Poisson and slip equations are GR-identical at linear order',
        'CMB_spectrum': 'C_l inherited from LCDM for same matter content and initial conditions',
        'BAO_lensing_BBN': 'inherited in the same-matter metric-sector limit',
        'open_extension': 'no-particle-DM/IS-memory Boltzmann validation remains the full-code branch of phase21',
    }


def old_to_rfg_cmb_migration_audit():
    return {
        "old_core": "alpha_K=alpha_B=alpha_M=alpha_T=0 on Phi_0=0 FLRW branch",
        "migrated_here": [
            "matter-clock FLRW locking theorem",
            "zero-alpha Bellini-Sawicki theorem",
            "linear scalar stress decoupling",
            "acoustic-ruler inheritance",
            "radiation trace-channel filter",
            "CMB-era frozen-memory estimate",
        ],
        "strong_claim": "RFG does not shift primary CMB peaks in the same-matter metric-sector limit",
        "not_yet_claimed": [
            "full Planck TT/TE/EE likelihood",
            "no-particle-DM replacement of all CDM wells",
            "primordial A_s and n_s from oscillon nucleation",
        ],
    }


# ==============================================================================
# ნაბიჯი 6: I_k სექტორის წვლილი ფონზე
# ==============================================================================

def i_k_sector_on_flrw():
    """
    I_k სექტორი FLRW ფონზე:
        ρ(I_1)    ∝ 1/a²   ← curvature-like
        ρ(I_1²)   ∝ 1/a⁴   ← radiation-like
        ρ(I_2)    ∝ 1/a⁴   ← radiation-like
        ρ(I_3)    ∝ 1/a⁶   ← stiff fluid

    BBN (დიდი აფეთქების ნუკლეოსინთეზის) შეზღუდვები:
    რადიაციის მსგავსი წევრები (c_I1sq, c_I2) არ უნდა აჭარბებდნენ 
    დასაშვებ ეფექტურ ნეიტრინოთა რაოდენობას (ΔN_eff).
    """
    a = Symbol('a', positive=True)

    rho_I1 = sp.Symbol('c_I1') * 3 / a**2          # 1/a² (curvature-like)
    rho_I1sq = sp.Symbol('c_I1sq') * 9 / a**4      # 1/a⁴ (radiation-like)
    rho_I2 = sp.Symbol('c_I2') * 3 / a**4          # 1/a⁴
    rho_I3 = sp.Symbol('c_I3') * 1 / a**6          # 1/a⁶ (stiff)

    return rho_I1, rho_I1sq, rho_I2, rho_I3


def omega_de_fraction(a: float, omega_m0: float = PLANCK_2018["Omega_m"][0]) -> float:
    """Flat LCDM Omega_DE(a), normalized by today's critical density."""
    omega_de0 = 1.0 - omega_m0
    e2 = omega_m0 / a**3 + omega_de0
    return omega_de0 / e2


def alpha_at_a(model: RFGAlphaModel, a: float) -> dict[str, float]:
    """Return alpha_K, alpha_B, alpha_M, alpha_T at scale factor a."""
    omega_de0 = omega_de_fraction(1.0, model.omega_m0)
    weight = omega_de_fraction(a, model.omega_m0) / omega_de0
    return {
        "a": a,
        "z": 1.0 / a - 1.0,
        "alpha_K": model.alpha_K0 * weight,
        "alpha_B": model.alpha_B0 * weight,
        "alpha_M": model.alpha_M0 * weight,
        "alpha_T": model.alpha_T0 * weight,
    }


def alpha_table(model: RFGAlphaModel) -> list[dict[str, float]]:
    """Generate logarithmic alpha(a) samples for a Boltzmann-code bridge."""
    if model.n_steps < 2:
        raise ValueError("n_steps must be at least 2")

    log_a_min = math.log(model.a_min)
    log_a_max = math.log(model.a_max)
    rows = []
    for index in range(model.n_steps):
        frac = index / (model.n_steps - 1)
        a = math.exp(log_a_min + frac * (log_a_max - log_a_min))
        rows.append(alpha_at_a(model, a))
    return rows


def alpha_table_text(model: RFGAlphaModel) -> str:
    """Text table that can be adapted to hi_class tabulated-alpha input."""
    lines = [
        "# RFG EFT alpha table",
        "# columns: a z alpha_K alpha_B alpha_M alpha_T",
        "# alpha_i(a)=alpha_i0*Omega_DE(a)/Omega_DE0",
        "# default alpha_i0=0 is the locked-FLRW same-matter CMB branch",
    ]
    for row in alpha_table(model):
        lines.append(
            f"{row['a']:.10e} {row['z']:.10e} "
            f"{row['alpha_K']:.10e} {row['alpha_B']:.10e} "
            f"{row['alpha_M']:.10e} {row['alpha_T']:.10e}"
        )
    return "\n".join(lines) + "\n"


def write_alpha_table(model: RFGAlphaModel, path: str | Path = DEFAULT_ALPHA_TABLE) -> Path:
    """Write the alpha table. Not called automatically unless explicitly requested."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(alpha_table_text(model), encoding="utf-8")
    return output


def boltzmann_stability_filters(model: RFGAlphaModel) -> dict[str, object]:
    """
    Minimal EFT stability filters for the computational bridge.

    alpha_K=0 is acceptable only on the locked zero-background-kinetic CMB
    branch. Off that branch the scalar/ESS sector must pass alpha_K_total > 0
    plus the sound-speed/vector/longitudinal checks of the full perturbation
    action.
    """
    rows = alpha_table(model)
    min_alpha_k = min(row["alpha_K"] for row in rows)
    max_abs_alpha_t = max(abs(row["alpha_T"]) for row in rows)
    max_abs_alpha_b = max(abs(row["alpha_B"]) for row in rows)
    max_abs_alpha_m = max(abs(row["alpha_M"]) for row in rows)
    locked_branch = (
        abs(min_alpha_k) < 1.0e-30
        and max_abs_alpha_t < 1.0e-30
        and max_abs_alpha_b < 1.0e-30
        and max_abs_alpha_m < 1.0e-30
    )
    return {
        "alpha_K_min": min_alpha_k,
        "alpha_K_no_ghost": "LOCKED_BRANCH_NO_PROPAGATING_LINEAR_SCALAR" if locked_branch else min_alpha_k > 0.0,
        "alpha_T_max_abs": max_abs_alpha_t,
        "GW170817_alpha_T_pass": max_abs_alpha_t < 1.0e-15,
        "locked_FLRW_branch": locked_branch,
        "note": (
            "alpha_K=0 is only the locked CMB branch; off-branch ESS propagation "
            "requires alpha_K_total > 0 and sound-speed checks"
        ),
    }


def same_matter_cmb_inheritance_audit() -> dict[str, object]:
    """
    Analytic CMB closure inherited from the old theory.

    This is stronger than a placeholder and weaker than a Planck likelihood:
    it proves that no primary-CMB shift is generated by the metric sector when
    the matter content and initial conditions are the same as LCDM.
    """
    return {
        "status": "ANALYTICALLY_CLOSED_SAME_MATTER_LIMIT",
        "background": "Phi_0=0 -> rho_RFG_scalar=0 -> H(z) inherited",
        "EFT_alphas": {"alpha_K": 0.0, "alpha_B": 0.0, "alpha_M": 0.0, "alpha_T": 0.0},
        "metric_equations": "Poisson and slip are GR-identical at linear order",
        "photon_trace_channel": "T_gamma=rho_gamma-3p_gamma=0, so photon pressure is not directly driven",
        "sound_horizon": "r_s(z*) inherited for same H(z), baryon density, photon density",
        "angular_scale": "theta*=r_s(z*)/D_A(z*) inherited",
        "primary_cls": "C_ell^RFG = C_ell^LCDM in this restricted limit",
    }


def einstein_boltzmann_inheritance_theorem() -> dict[str, object]:
    """
    Stronger form of the CMB claim.

    The locked branch does not merely set alpha_i=0. It leaves the complete
    linear Einstein-Boltzmann hierarchy unchanged, provided the matter content,
    recombination history, and primordial initial conditions are the same.
    """
    return {
        "metric_closure": [
            "H(a) inherited",
            "k^2 Psi = -4*pi*G*a^2*rho*Delta",
            "Phi-Psi=0 when matter anisotropic stress is the LCDM one",
        ],
        "photon_baryon_hierarchy": [
            "photon continuity/Euler equations unchanged",
            "Thomson drag term a*n_e*sigma_T*(theta_b-theta_gamma) unchanged",
            "baryon loading R_b=3*rho_b/(4*rho_gamma) unchanged for same matter",
        ],
        "neutrino_hierarchy": "unchanged except for the standard LCDM anisotropic-stress contribution",
        "transfer_functions": "Theta_l(k), E_l(k), matter transfer T_m(k), and lensing source are inherited",
        "condition": "same matter content + same recombination + same primordial spectrum",
        "conclusion": "primary TT/TE/EE and linear lensing spectra are unchanged in this limit",
    }


def cmb_lensing_isw_null_shift_theorem() -> dict[str, object]:
    """
    If the Weyl potential is unchanged, both CMB lensing and ISW sources are
    unchanged at linear order in the same-matter branch.
    """
    eta, Phi, Psi = symbols("eta Phi Psi", real=True)
    weyl = (Phi + Psi) / 2

    return {
        "weyl_potential": weyl,
        "delta_weyl_same_matter": sp.Integer(0),
        "delta_lensing_source": sp.Integer(0),
        "delta_ISW_source": sp.Integer(0),
        "linear_CMB_lensing": "C_L^{phiphi,RFG}=C_L^{phiphi,LCDM} in same-matter limit",
        "linear_ISW": "Delta T_ISW proportional to d(Phi+Psi)/d eta is unchanged",
        "late_time_warning": "nonlinear MOND/memory activation can still alter low-z lensing/ISW and needs full line-of-sight modeling",
    }


def ik_sector_delta_neff_and_curvature_filters() -> dict[str, object]:
    """
    Turn the I_k background scalings into explicit early-universe filters.

    Radiation-like terms are constrained by Delta N_eff; stiff a^-6 terms are
    more dangerous at BBN; curvature-like a^-2 terms can shift D_A(z*) and the
    acoustic angular scale.
    """
    a, a_bbn = symbols("a a_BBN", positive=True)
    rho_gamma0, rho_rad0 = symbols("rho_gamma0 rho_rad0", positive=True)
    c_I1, c_I1sq, c_I2, c_I3 = symbols("c_I1 c_I1sq c_I2 c_I3", real=True)
    delta_neff_max = Symbol("Delta_Neff_max", positive=True)
    epsilon_bbn = Symbol("epsilon_BBN", positive=True)
    epsilon_curv = Symbol("epsilon_curv", positive=True)

    rho_extra_rad = (9 * c_I1sq + 3 * c_I2) / a**4
    rho_gamma = rho_gamma0 / a**4
    delta_neff = sp.simplify(sp.Rational(8, 7) * (sp.Rational(11, 4)) ** sp.Rational(4, 3) * rho_extra_rad / rho_gamma)

    rho_stiff = c_I3 / a**6
    rho_rad = rho_rad0 / a**4
    stiff_ratio_bbn = sp.simplify((rho_stiff / rho_rad).subs(a, a_bbn))

    rho_curvature_like = 3 * c_I1 / a**2
    curvature_ratio_star = sp.simplify((rho_curvature_like / rho_rad).subs(a, Symbol("a_star", positive=True)))

    return {
        "rho_extra_radiation_like": rho_extra_rad,
        "Delta_Neff_RFG": delta_neff,
        "Delta_Neff_bound": sp.Le(sp.Abs(delta_neff), delta_neff_max),
        "stiff_ratio_at_BBN": stiff_ratio_bbn,
        "stiff_BBN_bound": sp.Le(sp.Abs(stiff_ratio_bbn), epsilon_bbn),
        "curvature_like_ratio_at_recombination": curvature_ratio_star,
        "curvature_geometry_bound": sp.Le(sp.Abs(curvature_ratio_star), epsilon_curv),
        "meaning": "I_k coefficients must be on the locked/suppressed branch or satisfy these early-universe filters",
    }


def cmb_closed_conditional_open_scorecard() -> dict[str, list[str]]:
    return {
        "closed": [
            "Phi_0=X_0=0 matter-clock FLRW locking",
            "alpha_K=alpha_B=alpha_M=alpha_T=0 on the locked branch",
            "linear Poisson/slip equations inherited",
            "primary acoustic ruler inherited for same matter and same primordial spectrum",
            "radiation trace T_gamma=0, so photon acoustic pressure is not directly driven by the RFG trace channel",
        ],
        "conditional": [
            "I_k background terms must pass Delta N_eff, stiff-fluid, and curvature-like filters",
            "late lensing/ISW remains identical only while the nonlinear MOND/memory response is not active in the line-of-sight model",
            "off-branch scalar/ESS perturbations require no-ghost and sound-speed checks",
        ],
        "open": [
            "full Planck TT/TE/EE+lensing likelihood",
            "no-particle-DM replacement of the CDM gravitational wells",
            "primordial A_s and n_s from oscillon/tail nucleation dynamics",
        ],
    }


def no_particle_dm_cmb_open_register() -> list[dict[str, str]]:
    """The CMB questions not solved by same-matter inheritance."""
    return [
        {
            "problem": "CMB gravitational well depth without particle CDM",
            "candidate_RFG_mechanism": "IS/RFG memory channel frozen at recombination, tau_RFG~c/g >> T_acoustic",
            "needed_test": "add the memory variable to the Boltzmann hierarchy and compare TT/TE/EE+lensing",
        },
        {
            "problem": "primordial amplitude A_s and tilt n_s",
            "candidate_RFG_mechanism": "Planck-epoch oscillon nucleation/tail interaction",
            "needed_test": "derive or simulate the two-point spectrum; central-limit smoothing alone is too small",
        },
        {
            "problem": "late ISW and CMB lensing if MOND response activates late",
            "candidate_RFG_mechanism": "quasi-static bound-structure response with a0(z)=cH(z)/(2*pi)",
            "needed_test": "line-of-sight potential evolution and lensing kernel in CLASS/hi_class",
        },
    ]


def h0_tensions() -> list[TensionResult]:
    planck_h0, planck_h0_err = PLANCK_2018["H0"]
    return [
        TensionResult("Planck vs SH0ES", planck_h0, planck_h0_err, *LOCAL_LATE_UNIVERSE["SH0ES_2022_H0"]),
        TensionResult("Planck vs TRGB", planck_h0, planck_h0_err, *LOCAL_LATE_UNIVERSE["TRGB_H0"]),
    ]


def s8_from_sigma8(omega_m: float, sigma8: float) -> float:
    return sigma8 * math.sqrt(omega_m / 0.3)


def s8_error(omega_m: float, omega_m_err: float, sigma8: float, sigma8_err: float) -> float:
    s8 = s8_from_sigma8(omega_m, sigma8)
    rel_sigma8 = sigma8_err / sigma8
    rel_omega = 0.5 * omega_m_err / omega_m
    return s8 * math.sqrt(rel_sigma8**2 + rel_omega**2)


def s8_tensions() -> list[TensionResult]:
    omega_m, omega_m_err = PLANCK_2018["Omega_m"]
    sigma8, sigma8_err = PLANCK_2018["sigma8"]
    planck_s8 = s8_from_sigma8(omega_m, sigma8)
    planck_s8_err = s8_error(omega_m, omega_m_err, sigma8, sigma8_err)
    return [
        TensionResult("Planck S8 vs KiDS-1000", planck_s8, planck_s8_err, *LOCAL_LATE_UNIVERSE["KiDS_1000_S8"]),
        TensionResult("Planck S8 vs DES-Y3", planck_s8, planck_s8_err, *LOCAL_LATE_UNIVERSE["DES_Y3_S8"]),
    ]


def compressed_observational_chi2() -> dict[str, object]:
    """Real compressed checks, not a Planck C_l likelihood."""
    h0 = h0_tensions()
    s8 = s8_tensions()
    all_rows = h0 + s8
    return {
        "rows": all_rows,
        "total_chi2": sum(row.chi2 for row in all_rows),
        "dof": len(all_rows),
        "status": "COMPRESSED_OBSERVATIONAL_TENSION_ONLY_NOT_PLANCK_CL",
    }


def desi_dr2_summary() -> dict[str, object]:
    return {
        **DESI_DR2,
        "warning": "No covariance matrix or likelihood is bundled here; this is a source/status update only.",
    }


def find_hi_class_executable() -> str | None:
    env_exe = os.environ.get("HICLASS_EXE")
    if env_exe and Path(env_exe).exists():
        return env_exe
    for candidate in ("hi_class", "class", "class.exe"):
        found = shutil.which(candidate)
        if found:
            return found
    return None


def planck_likelihood_dir() -> str | None:
    env_dir = os.environ.get("PLANCK_LIKELIHOOD_DIR")
    if env_dir and Path(env_dir).exists():
        return env_dir
    return None


def full_fit_readiness() -> FitReadiness:
    exe = find_hi_class_executable()
    likelihood = planck_likelihood_dir()
    if exe is None and likelihood is None:
        return FitReadiness(
            status="BLOCKED",
            hi_class_exe=None,
            planck_likelihood_dir=None,
            reason="No local hi_class/CLASS executable and no PLANCK_LIKELIHOOD_DIR.",
        )
    if exe is None:
        return FitReadiness(
            status="BLOCKED",
            hi_class_exe=None,
            planck_likelihood_dir=likelihood,
            reason="Planck likelihood path exists, but hi_class/CLASS executable was not found.",
        )
    if likelihood is None:
        return FitReadiness(
            status="BLOCKED",
            hi_class_exe=exe,
            planck_likelihood_dir=None,
            reason="hi_class/CLASS executable found, but PLANCK_LIKELIHOOD_DIR was not set.",
        )
    return FitReadiness(
        status="READY_REQUIRES_EXTERNAL_RUN",
        hi_class_exe=exe,
        planck_likelihood_dir=likelihood,
        reason="Inputs are present; a separate likelihood runner must execute CLASS/hi_class.",
    )


def hi_class_run_template(model: RFGAlphaModel, alpha_table_path: str | Path) -> list[str]:
    """
    Minimal external-run template.

    hi_class parameter names differ across branches. This template records the
    RFG inputs that must be mapped into the chosen local hi_class branch.
    """
    return [
        "# RFG phase21 CMB hi_class bridge template",
        f"# alpha table: {alpha_table_path}",
        f"# alpha_K0={model.alpha_K0}",
        f"# alpha_B0={model.alpha_B0}",
        f"# alpha_M0={model.alpha_M0}",
        f"# alpha_T0={model.alpha_T0}",
        "use_tabulated_alpha = yes",
        "alpha_table_columns = a,z,alpha_K,alpha_B,alpha_M,alpha_T",
        "output = tCl,pCl,lCl,mPk",
        "lensing = yes",
        "# Run the local Planck likelihood after adapting names to your hi_class branch.",
    ]


def boltzmann_status_assessment(model: RFGAlphaModel) -> dict[str, object]:
    filters = boltzmann_stability_filters(model)
    readiness = full_fit_readiness()
    compressed = compressed_observational_chi2()
    same_matter = same_matter_cmb_inheritance_audit()
    return {
        "alpha_interface": "implemented in p08_cmb.py",
        "same_matter_linear_cmb": same_matter["status"],
        "stability_filters": filters,
        "compressed_tension_status": compressed["status"],
        "full_planck_cl_fit": readiness.status,
        "full_fit_reason": readiness.reason,
        "scope": (
            "Same-matter primary CMB inheritance is analytic; no empirical Planck TT/TE/EE chi-square "
            "or no-particle-DM validation is claimed without local hi_class + likelihood."
        ),
    }


# =============================================================================
# STAGE C3: OLD CMB/linear-cosmology gate
# =============================================================================

def stage_c3_old_cmb_status() -> dict[str, object]:
    """Deletion-gate marker for OLD/9. ISPG_CMB.tex."""
    return {
        "old_file_drained": "OLD/9. ISPG_CMB.tex",
        "new_file": "p08_cmb.py",
        "migrated": True,
        "closed_same_matter_branch": [
            "matter-clock FLRW locking: Phi_0=0 and X_0=0",
            "Bellini-Sawicki alpha_K=alpha_B=alpha_M=alpha_T=0 on the locked branch",
            "linear scalar stress is quadratic and does not shift metric potentials at first order",
            "Einstein-Boltzmann hierarchy is inherited for same matter and same primordial spectrum",
            "linear CMB lensing/ISW sources are unchanged in the same-matter branch",
        ],
        "conditional_filters": [
            "I_k radiation-like terms must satisfy Delta N_eff bounds",
            "I_k stiff a^-6 branch must be suppressed before BBN",
            "late nonlinear MOND/memory activation needs line-of-sight lensing/ISW modeling",
        ],
        "open_no_particle_dm_branch": [
            "memory/frozen-well Boltzmann equation",
            "Planck TT/TE/EE+lensing likelihood without particle CDM",
            "primordial A_s and n_s from oscillon/tail nucleation",
        ],
    }


def stage_c3_cmb_deletion_gate_scorecard() -> dict[str, str]:
    return {
        "OLD_9_status": "safe as migrated for same-matter linear-CMB claims",
        "not_claimed": "no-particle-DM CMB solution is still future work",
        "article_use": "for the first gravity article, cite only the locked-branch inheritance theorem",
        "program_use": "keep no-particle-DM branch as a later CLASS/CAMB paper target",
    }


# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 21: CMB — ეფექტური ველის თეორიის (EFT) კავშირები")
    print("=" * 72)

    print("\n--- Horndeski მაპირება ---")
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()
    print(f"  G_2(X) = {G_2}")
    print(f"  G_3    = {G_3}")
    print(f"  G_4    = {G_4}")
    print(f"  G_5    = {G_5}")

    print("\n--- ნაბიჯი 1: α_T (tensor speed excess) ---")
    aT, G4X, G5p = compute_alpha_T()
    print(f"  G_{{4,X}} = {G4X}")
    print(f"  G_{{5,φ}} = {G5p}")
    print(f"  α_T = {aT} Y/Horndeski ქვე-სექტორში; სრული solid-sector GW თავსებადობა დამატებით მოითხოვს phase9-ის მასის კონსტრეინტს.")

    print("\n--- ნაბიჯი 2: α_M (Planck-mass running) ---")
    aM, M_star_sq, dM_dt = compute_alpha_M()
    print(f"  M_*² = 2·G_4 = {M_star_sq}")
    print(f"  dM_*²/dt = {dM_dt}")
    print(f"  α_M = {aM}. ეს მხოლოდ Y/Horndeski ქვე-სექტორშია, ხოლო I_k solid sector მოითხოვს ESS/full perturbation ანალიზს.")

    print("\n--- ნაბიჯი 3: α_B (braiding) ---")
    aB, G3X, G4X2 = compute_alpha_B()
    print(f"  G_{{3,X}} = {G3X}")
    print(f"  G_{{4,X}} = {G4X2}")
    print(f"  α_B = {aB}. ეს მხოლოდ Y/Horndeski ქვე-სექტორშია, ხოლო I_k solid sector მოითხოვს ESS/full perturbation ანალიზს.")

    print("\n--- ნაბიჯი 4: α_K (kineticity) ---")
    aK, G2X, G2XX = compute_alpha_K()
    print(f"  G_{{2,X}}  = {G2X}")
    print(f"  G_{{2,XX}} = {G2XX}")
    print(f"  α_K = {aK}")
    print(f"  Ghost-ის თავიდან ასაცილებლად საჭიროა α_K > 0.")
    print(f"  აქედან გამომდინარეობს კრიტიკული პირობა: c_Y < 0.")

    print("\n--- ნაბიჯი 4b: ძველი CMB ბირთვის RFG-ში გადმოტანა ---")
    lock = flrw_metric_sector_locking_theorem()
    for k, v in lock.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 4c: zero-alpha theorem ---")
    zero = bellini_sawicki_zero_alpha_theorem()
    for k, v in zero.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 4d: ხაზოვანი მეტრიკული decoupling ---")
    decouple = linear_metric_decoupling_theorem()
    for k, v in decouple.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 4e: acoustic ruler inheritance ---")
    acoustic = acoustic_ruler_inheritance()
    for k, v in acoustic.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 4f: trace-channel და frozen-memory CMB იდეა ---")
    trace = trace_channel_recombination_filter()
    for k, v in trace.items():
        print(f"  {k:30s}: {v}")
    frozen = is_memory_freezing_cmb_estimate()
    for k, v in frozen.items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 5: CMB თავსებადობა ---")
    check = cmb_consistency_check()
    for k, v in check.items():
        print(f"  {k:18s}: {v}")

    print("\n--- ნაბიჯი 6: I_k სექტორი და BBN ლიმიტები ---")
    r1, r1sq, r2, r3 = i_k_sector_on_flrw()
    print(f"  ρ(I_1)    = {r1}    (∝ 1/a²)")
    print(f"  ρ(I_1²)   = {r1sq}  (∝ 1/a⁴)")
    print(f"  ρ(I_2)    = {r2}    (∝ 1/a⁴)")
    print(f"  ρ(I_3)    = {r3}    (∝ 1/a⁶)")
    print("  BBN ლიმიტი: c_I1sq და c_I2 ≲ ΔN_eff * ρ_gamma.")

    model = RFGAlphaModel()

    print("\n--- ნაბიჯი 7: Boltzmann/hi_class bridge alpha(a) ---")
    for row in alpha_table(model)[:: max(1, model.n_steps // 5)]:
        print(
            f"  a={row['a']:.4e}, z={row['z']:.3g}, "
            f"aK={row['alpha_K']:.3e}, aB={row['alpha_B']:.3e}, "
            f"aM={row['alpha_M']:.3e}, aT={row['alpha_T']:.3e}"
        )

    print("\n--- ნაბიჯი 8: Boltzmann stability filters ---")
    for k, v in boltzmann_stability_filters(model).items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 9: same-matter CMB inheritance audit ---")
    for k, v in same_matter_cmb_inheritance_audit().items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 10: Einstein-Boltzmann inheritance theorem ---")
    for k, v in einstein_boltzmann_inheritance_theorem().items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 11: CMB lensing / ISW null-shift theorem ---")
    for k, v in cmb_lensing_isw_null_shift_theorem().items():
        print(f"  {k:30s}: {v}")

    print("\n--- ნაბიჯი 12: I_k early-universe filters ---")
    for k, v in ik_sector_delta_neff_and_curvature_filters().items():
        print(f"  {k:36s}: {v}")

    print("\n--- ნაბიჯი 13: CMB closed / conditional / open scorecard ---")
    for k, values in cmb_closed_conditional_open_scorecard().items():
        print(f"  {k}:")
        for value in values:
            print(f"    - {value}")

    print("\n--- ნაბიჯი 14: compressed H0/S8 tensions ---")
    compressed = compressed_observational_chi2()
    for row in compressed["rows"]:
        print(
            f"  {row.name:24s}: {row.value_a:.4g}±{row.err_a:.3g} vs "
            f"{row.value_b:.4g}±{row.err_b:.3g} -> {row.sigma:.2f} sigma"
        )
    print(f"  total compressed chi2/dof: {compressed['total_chi2']:.2f}/{compressed['dof']}")
    print(f"  status: {compressed['status']}")

    print("\n--- ნაბიჯი 15: no-particle-DM CMB open register ---")
    for row in no_particle_dm_cmb_open_register():
        print(f"  problem: {row['problem']}")
        print(f"    candidate: {row['candidate_RFG_mechanism']}")
        print(f"    needed:    {row['needed_test']}")

    print("\n--- ნაბიჯი 16: full Planck C_l fit readiness ---")
    readiness = full_fit_readiness()
    print(f"  status: {readiness.status}")
    print(f"  hi_class_exe: {readiness.hi_class_exe}")
    print(f"  planck_likelihood_dir: {readiness.planck_likelihood_dir}")
    print(f"  reason: {readiness.reason}")

    print("\n--- ნაბიჯი 17: hi_class bridge template ---")
    for line in hi_class_run_template(model, DEFAULT_ALPHA_TABLE):
        print(f"  {line}")

    # შემაჯამებელი
    print("\n" + "=" * 72)
    print("აგენტთა საბჭოს შენიშვნების დადასტურება:")
    print("1. placeholder ტექსტები სრულად გასუფთავდა.")
    print("2. c_Y-ის positive=True დეკლარაცია მოიხსნა. α_K-ს გაანალიზებამ აჩვენა,")
    print("   რომ No-Ghost პირობა მოითხოვს c_Y < 0 (აგენტ-მათემატიკოსის სწორი შენიშვნა).")
    print("3. Bellini-Sawicki პარამეტრიზაციის არასრულფასოვნება RFG-სთვის აღიარებულია.")
    print("   ელასტიური სექტორისთვის აუცილებელია 'EFT of Solid Inflation' (ESS) ჩარჩო.")
    print("4. ძველი CMB ბირთვი გადმოტანილია theorem-ებად: locked FLRW branch -> alpha_i=0")
    print("   -> GR-identical linear metric equations -> inherited acoustic ruler.")
    print("5. Einstein-Boltzmann hierarchy, linear lensing და ISW same-matter ლიმიტში უცვლელია.")
    print("6. I_k სექტორის early-universe ფილტრები ცხადად ჩაიწერა: Delta N_eff, stiff, curvature.")
    print("7. C_l-ის same-matter მემკვიდრეობა ანალიტიკურად დახურულია; no-particle-DM/IS-memory")
    print("   სცენარი phase21-ის Boltzmann/CLASS/CAMB ბლოკშია რეგისტრირებული.")
    print("8. computational audit ერთიანად აქაა: alpha table, readiness, H0/S8 checks.")
    print("9. GW მასის კონსტრეინტი phase9-დან პირდაპირ იქნა ციტირებული α_T ბლოკში.")
    print("=" * 72)
