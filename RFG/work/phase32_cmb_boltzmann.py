# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 32: CMB Boltzmann — hi_class ინტერფეისი და კომპრესირებული შემოწმებები
================================================================================

რეფერენცია: phase21_cmb.py, phase23_horndeski_dhost_embedding.py,
            phase25_full_eft_de.py, STRATEGY.md S3/E7

სტატუსი:
ეს ფაილი არ აცხადებს სრულ Planck TT/TE/EE likelihood გაშვებას. რეალური Planck 
chi-square მოითხოვს ლოკალურ CLASS/hi_class გარემოს და Planck-ის 
მონაცემებს. მიმდინარე რეპოზიტორიუმში არცერთი მათგანი არ ინახება.

რა არის იმპლემენტირებული:
- RFG EFT alpha(a) ცხრილის გენერატორი hi_class/CLASS-თან ინტეგრაციისთვის.
- GW170817 alpha_T ფილტრი და alpha_K no-ghost ფილტრი.
- კომპრესირებული ობსერვაციული შემოწმებები H0 და S8/sigma8 tension-ებისთვის.
- Blocked/full-fit სტატუსის დეტექტორი, რათა კოდმა ჩუმად არ გამოაცხადოს, 
  თითქოს Planck C_l likelihood მართლა გაეშვა.
"""

from __future__ import annotations

import math
import os
import shutil
from dataclasses import dataclass
from pathlib import Path


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

DEFAULT_ALPHA_TABLE = Path("RFG/work/phase32_rfg_alpha_table.dat")


@dataclass
class RFGAlphaModel:
    """
    Common EFT parameterization: alpha_i(a) = alpha_i0 * Omega_DE(a)/Omega_DE0.

    This is a hi_class-ready phenomenological bridge, not a derivation of the
    ESS solid-sector perturbations. phase25 keeps those deltas symbolic.
    """

    alpha_K0: float = 0.10
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
        "# alpha_i(a)=alpha_i0*Omega_DE(a)/Omega_DE0; phenomenological bridge",
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


def stability_filters(model: RFGAlphaModel) -> dict[str, object]:
    """
    Minimal EFT stability filters.

    Full scalar/sound-speed stability requires the actual ESS perturbation
    action. Here we enforce only filters already fixed in phase25/phase34.
    """
    rows = alpha_table(model)
    min_alpha_k = min(row["alpha_K"] for row in rows)
    max_abs_alpha_t = max(abs(row["alpha_T"]) for row in rows)
    return {
        "alpha_K_min": min_alpha_k,
        "alpha_K_no_ghost": min_alpha_k > 0.0,
        "alpha_T_max_abs": max_abs_alpha_t,
        "GW170817_alpha_T_pass": max_abs_alpha_t < 1.0e-15,
        "note": "sound speed and ESS vector/longitudinal stability require a full perturbation derivation",
    }


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
    """
    Real compressed checks, not a Planck C_l likelihood.

    The output quantifies the existing H0/S8 tensions that phase32 must address
    once a real hi_class likelihood is available.
    """
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
        "# RFG phase32 hi_class bridge template",
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


def status_assessment(model: RFGAlphaModel) -> dict[str, object]:
    filters = stability_filters(model)
    readiness = full_fit_readiness()
    compressed = compressed_observational_chi2()
    return {
        "alpha_interface": "implemented",
        "stability_filters": filters,
        "compressed_tension_status": compressed["status"],
        "full_planck_cl_fit": readiness.status,
        "full_fit_reason": readiness.reason,
        "scope": "No empirical Planck TT/TE/EE chi-square is claimed without local hi_class + likelihood.",
    }


def main() -> None:
    print("=" * 72)
    print("PHASE 32: CMB Boltzmann — hi_class interface and compressed checks")
    print("=" * 72)

    model = RFGAlphaModel()

    print("\n1. Planck 2018 compressed inputs")
    for key, value in PLANCK_2018.items():
        print(f"  {key:18s}: {value}")

    print("\n2. RFG alpha(a) sample")
    for row in alpha_table(model)[:: max(1, model.n_steps // 5)]:
        print(
            f"  a={row['a']:.4e}, z={row['z']:.3g}, "
            f"aK={row['alpha_K']:.3e}, aB={row['alpha_B']:.3e}, "
            f"aM={row['alpha_M']:.3e}, aT={row['alpha_T']:.3e}"
        )

    print("\n3. Stability filters")
    for key, value in stability_filters(model).items():
        print(f"  {key:28s}: {value}")

    print("\n4. Compressed H0/S8 tensions")
    compressed = compressed_observational_chi2()
    for row in compressed["rows"]:
        print(
            f"  {row.name:24s}: {row.value_a:.4g}±{row.err_a:.3g} vs "
            f"{row.value_b:.4g}±{row.err_b:.3g} -> {row.sigma:.2f} sigma"
        )
    print(f"  total compressed chi2/dof: {compressed['total_chi2']:.2f}/{compressed['dof']}")
    print(f"  status: {compressed['status']}")

    print("\n5. DESI DR2 dark-energy summary")
    for key, value in desi_dr2_summary().items():
        print(f"  {key:28s}: {value}")

    print("\n6. Full Planck C_l fit readiness")
    readiness = full_fit_readiness()
    print(f"  status: {readiness.status}")
    print(f"  hi_class_exe: {readiness.hi_class_exe}")
    print(f"  planck_likelihood_dir: {readiness.planck_likelihood_dir}")
    print(f"  reason: {readiness.reason}")

    print("\n7. hi_class bridge template")
    for line in hi_class_run_template(model, DEFAULT_ALPHA_TABLE):
        print(f"  {line}")

    print("\n8. Status")
    for key, value in status_assessment(model).items():
        print(f"  {key:24s}: {value}")


if __name__ == "__main__":
    main()
