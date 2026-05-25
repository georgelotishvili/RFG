# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 33: SPARC ბრუნვის მრუდები — RFG/AQUAL chi-square ფაიფლაინი
================================================================================

რეფერენცია: phase11_mond_metric.py, STRATEGY.md S3/E5

სტატუსი:
ეს ფაილი წარმოადგენს SPARC ფიტინგის გაშვებად ფაიფლაინს და არა უკვე მიღებულ 
შედეგს. მიმდინარე რეპოზიტორიუმში არ არის ატვირთული 175 SPARC ბრუნვის მრუდის ფაილი, 
ამიტომ სკრიპტმა არ უნდა გამოაცხადოს ემპირიული SPARC chi-square შედეგი, სანამ ეს 
ფაილები ლოკალურად არ იარსებებს.

რა არის იმპლემენტირებული:
- SPARC rotmod-სტილის მონაცემების ჩამტვირთავი.
- AQUAL/MOND სიჩქარის მოდელი mu(x) = x/(1+x)-ით, რაც ეკვივალენტურია
  nu(y) = 0.5 * (1 + sqrt(1 + 4/y))-ს.
- გლობალური ბადისებრი ფიტი (grid fit) უნივერსალური a0-ისა და chi_coupling-ისთვის.
- RFG ანიზოტროპიის (Delta p) პროფილის პირდაპირი ექსტრაქცია დაკვირვებული სიჩქარეებიდან.
- ინდივიდუალური გალაქტიკების nuisance პარამეტრების ფიტი (მასა-ნათობის ფარდობა).
- სინთეტიკური smoke-test, რომელიც ამოწმებს, შეუძლია თუ არა ფაიფლაინს ხელოვნურად 
  ჩასმული a0 და chi_coupling მნიშვნელობების აღდგენა. ეს არ არის ემპირიული შედეგი.

მოსალოდნელი ლოკალური მონაცემების სტრუქტურა:
    RFG/data/SPARC/*_rotmod.dat
    RFG/work/data/SPARC/*_rotmod.dat
    data/SPARC/*_rotmod.dat
ან გამოიყენეთ SPARC_DATA_DIR გარემოს ცვლადი.
"""

from __future__ import annotations

import math
import os
import sympy as sp
from dataclasses import dataclass
from pathlib import Path


KPC_M = 3.0856775814913673e19
A0_CANONICAL = 1.2e-10
RHO_SOLID_COSMOLOGICAL = 1.0e-26  # kg/m^3 (დაახლოებით კოსმოლოგიური მუდმივას/ვაკუუმის სიმკვრივე)

SPARC_SUMMARY = {
    "n_galaxies": 175,
    "reference": "Lelli, McGaugh, Schombert 2016 (AJ 152:157)",
    "model_used_here": "AQUAL/Famaey-Binney mu(x)=x/(1+x)",
    "rfg_status": "extracts required supersolid anisotropy (Delta p); empirical verdict pending local SPARC data",
}

DEFAULT_DATA_DIRS = (
    Path("RFG/data/SPARC"),
    Path("RFG/work/data/SPARC"),
    Path("data/SPARC"),
    Path("SPARC"),
)


@dataclass
class RotationCurve:
    name: str
    radius_kpc: list[float]
    v_obs_km_s: list[float]
    v_err_km_s: list[float]
    v_gas_km_s: list[float]
    v_disk_km_s: list[float]
    v_bulge_km_s: list[float]

    @property
    def n_points(self) -> int:
        return len(self.radius_kpc)

    @property
    def has_bulge(self) -> bool:
        return any(abs(value) > 1e-9 for value in self.v_bulge_km_s)


@dataclass
class GalaxyFit:
    name: str
    chi2: float
    dof_local: int
    ml_disk: float
    ml_bulge: float
    nuisance_count: int
    n_points: int
    rms_residual_km_s: float

    @property
    def chi2_dof_local(self) -> float:
        return self.chi2 / self.dof_local if self.dof_local > 0 else math.inf


@dataclass
class GlobalFit:
    a0_m_s2: float
    chi_coupling: float
    total_chi2: float
    dof: int
    galaxy_fits: list[GalaxyFit]
    n_galaxies: int
    n_points: int
    status: str

    @property
    def chi2_dof(self) -> float:
        return self.total_chi2 / self.dof if self.dof > 0 else math.inf


def parse_numeric_row(line: str) -> list[float] | None:
    """Return a numeric table row; ignore comments and metadata."""
    stripped = line.strip()
    if not stripped or stripped[0] in "#!;":
        return None

    values = []
    for part in stripped.replace(",", " ").split():
        try:
            values.append(float(part))
        except ValueError:
            return None
    return values if len(values) >= 5 else None


def clean_galaxy_name(path: Path) -> str:
    name = path.stem
    lower = name.lower()
    for suffix in ("_rotmod", ".rotmod", "-rotmod"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def load_rotmod_file(path: Path) -> RotationCurve:
    """
    Load a SPARC rotmod-style file.

    Expected first columns:
        R[kpc], Vobs[km/s], e_Vobs[km/s], Vgas, Vdisk, Vbulge

    The gas/disk/bulge columns are velocity contributions; the baryonic
    Newtonian term is built as signed velocity squared:
        Vbar^2 = Vgas|Vgas| + Upsilon_d Vdisk|Vdisk| + Upsilon_b Vbul|Vbul|
    """
    radius: list[float] = []
    v_obs: list[float] = []
    v_err: list[float] = []
    v_gas: list[float] = []
    v_disk: list[float] = []
    v_bulge: list[float] = []

    text = path.read_text(encoding="utf-8", errors="ignore")
    for line in text.splitlines():
        row = parse_numeric_row(line)
        if row is None:
            continue

        r_kpc, obs, err, gas, disk = row[:5]
        bulge = row[5] if len(row) >= 6 else 0.0
        if r_kpc <= 0 or obs <= 0:
            continue

        radius.append(r_kpc)
        v_obs.append(obs)
        v_err.append(abs(err) if err != 0 else 1.0)
        v_gas.append(gas)
        v_disk.append(disk)
        v_bulge.append(bulge)

    if not radius:
        raise ValueError(f"no usable SPARC numeric rows found in {path}")

    return RotationCurve(
        name=clean_galaxy_name(path),
        radius_kpc=radius,
        v_obs_km_s=v_obs,
        v_err_km_s=v_err,
        v_gas_km_s=v_gas,
        v_disk_km_s=v_disk,
        v_bulge_km_s=v_bulge,
    )


def discover_sparc_files(data_dir: str | Path | None = None) -> list[Path]:
    """Find local SPARC rotmod files without downloading anything."""
    candidate_dirs: list[Path] = []
    if data_dir is not None:
        candidate_dirs.append(Path(data_dir))

    env_dir = os.environ.get("SPARC_DATA_DIR")
    if env_dir:
        candidate_dirs.append(Path(env_dir))

    candidate_dirs.extend(DEFAULT_DATA_DIRS)

    files: list[Path] = []
    seen: set[Path] = set()
    patterns = ("*rotmod*.dat", "*rotmod*.txt", "*.rotmod", "*.dat", "*.txt")

    for directory in candidate_dirs:
        if directory.is_file():
            resolved = directory.resolve()
            if resolved not in seen:
                files.append(directory)
                seen.add(resolved)
            continue
        if not directory.exists():
            continue
        for pattern in patterns:
            for path in sorted(directory.glob(pattern)):
                resolved = path.resolve()
                if resolved not in seen:
                    files.append(path)
                    seen.add(resolved)

    return sorted(files)


def load_sparc_dataset(
    data_dir: str | Path | None = None,
    max_galaxies: int | None = None,
) -> tuple[list[RotationCurve], list[str]]:
    """Load all usable local SPARC files; collect parse errors as warnings."""
    curves: list[RotationCurve] = []
    warnings: list[str] = []

    for path in discover_sparc_files(data_dir):
        try:
            curves.append(load_rotmod_file(path))
        except ValueError as exc:
            warnings.append(str(exc))
        if max_galaxies is not None and len(curves) >= max_galaxies:
            break

    return curves, warnings


def signed_square(value: float) -> float:
    """SPARC convention for gas components that can carry a negative sign."""
    return value * abs(value)


def baryonic_velocity_sq_km2_s2(
    curve: RotationCurve,
    ml_disk: float,
    ml_bulge: float,
) -> list[float]:
    return [
        signed_square(gas)
        + ml_disk * signed_square(disk)
        + ml_bulge * signed_square(bulge)
        for gas, disk, bulge in zip(
            curve.v_gas_km_s,
            curve.v_disk_km_s,
            curve.v_bulge_km_s,
        )
    ]


def nu_from_mu_simple(y: float) -> float:
    """
    Inverse interpolation for mu(x)=x/(1+x).

    With x = g/a0 and y = g_N/a0:
        y = x^2/(1+x)
        g = nu(y) g_N
        nu(y) = 0.5 * (1 + sqrt(1 + 4/y))
    """
    if y <= 0:
        return 0.0
    return 0.5 * (1.0 + math.sqrt(1.0 + 4.0 / y))


def predict_velocity_km_s(
    curve: RotationCurve,
    a0_m_s2: float,
    chi_coupling: float,
    ml_disk: float,
    ml_bulge: float,
) -> list[float]:
    """
    RFG/AQUAL rotation curve.

    chi_coupling = 0 gives the Newtonian baryonic curve.
    chi_coupling = 1 gives the AQUAL/MOND curve used in phase11.
    Values different from 1 are a phenomenological stress-test parameter.
    """
    v2_bar = baryonic_velocity_sq_km2_s2(curve, ml_disk, ml_bulge)
    predicted: list[float] = []

    for radius_kpc, v2_km2_s2 in zip(curve.radius_kpc, v2_bar):
        if radius_kpc <= 0 or v2_km2_s2 <= 0:
            predicted.append(0.0)
            continue

        radius_m = radius_kpc * KPC_M
        g_newton = v2_km2_s2 * 1.0e6 / radius_m
        y = g_newton / a0_m_s2
        g_aqual = g_newton * nu_from_mu_simple(y)
        g_model = g_newton + chi_coupling * (g_aqual - g_newton)

        predicted.append(math.sqrt(max(g_model * radius_m, 0.0)) / 1000.0)

    return predicted


def rfg_required_anisotropy_delta_p(
    curve: RotationCurve,
    ml_disk: float,
    ml_bulge: float,
) -> list[float]:
    """
    RFG-ის პირდაპირი მიკროფიზიკური გათვლა.
    ჩვენ ვიყენებთ RFG-ის ანიზოტროპიულ განტოლებას (phase1-დან):
        a_obs = v_obs^2 / r = g_N + 2 * Delta_p / (r * rho_solid)
    აქედან ვითვლით, რა ანიზოტროპიული სტრესი (Delta p) არის საჭირო ფუძე-მედიუმში, 
    რომ გალაქტიკამ იარსებოს ბნელი მატერიის გარეშე.
    """
    v2_bar = baryonic_velocity_sq_km2_s2(curve, ml_disk, ml_bulge)
    delta_p_list = []

    for radius_kpc, v_obs, v2_bar_val in zip(curve.radius_kpc, curve.v_obs_km_s, v2_bar):
        if radius_kpc <= 0:
            delta_p_list.append(0.0)
            continue

        radius_m = radius_kpc * KPC_M
        g_newton = v2_bar_val * 1.0e6 / radius_m
        a_obs = (v_obs * 1000.0)**2 / radius_m
        
        g_extra = a_obs - g_newton
        delta_p = g_extra * radius_m * RHO_SOLID_COSMOLOGICAL / 2.0
        delta_p_list.append(delta_p)

    return delta_p_list


def rfg_theoretical_anisotropy_scaling():
    """
    მათემატიკურად ამოწმებს ლაგრანჟიანიდან მიღებული სტატიკური ანიზოტროპიის
    თავსებადობას MOND-ის ბრტყელ მრუდებთან.
    """
    r, v_flat, rho_solid = sp.symbols('r v_flat rho_solid', real=True, positive=True)
    c_I1 = sp.Symbol('c_I1', real=True)
    
    # 1. MOND (ბრტყელი მრუდის) მოთხოვნა: a = v^2/r
    # 2*Delta_p / (r*rho) = v^2/r  =>  Delta_p = 1/2 rho v^2
    dp_mond_required = sp.Rational(1, 2) * rho_solid * v_flat**2
    
    # 2. სტატიკური სფერული მოდელის სკალირება phase1_tensor-დან (f(r)=r)
    # O(1) წევრები p_tan-სა და p_rad-ში იბათილება; რჩება გრადიენტები 1/r^2-დან.
    dp_static_scaling = c_I1 / r**2
    
    return {
        "required_for_flat_curves": dp_mond_required,
        "lagrangian_static_prediction": dp_static_scaling,
        "extra_acceleration_from_static": sp.simplify(2 * dp_static_scaling / (r * rho_solid)),
        "velocity_scaling_from_static": "v^2 = a*r ∝ 1/r^2, so v ∝ 1/r",
        "mathematical_conflict": "dp_static (1/r^2) ძალიან სწრაფად ქრება და ბრტყელ მრუდს ვერ იძლევა.",
        "physical_resolution": "ბრტყელი მრუდისთვის (v=const) ანიზოტროპია არ უნდა ქრებოდეს 1/r^2-ით. "
                               "დაკვირვებულ რადიუსებზე საჭიროა ვორტექსული plateau-რეჟიმი, "
                               "ხოლო საბოლოოდ finite-radius cutoff, რომ ენერგია უსასრულოდ არ გაიზარდოს."
    }


def rfg_a0_vortex_emergence():
    """
    აკავშირებს MOND-ის a_0 პარამეტრს RFG-ის ანიზოტროპიასთან (Delta p).
    აჩვენებს, რა ფორმის უნდა იყოს ვორტექსის სტრესი შორ მანძილებზე.
    """
    r, M, G, a0, rho_solid = sp.symbols('r M G a0 rho_solid', real=True, positive=True)
    Delta_p = sp.Symbol('Delta_p', real=True)
    
    # 1. MOND-ის ემპირიული აჩქარება Deep-MOND (შორ) ლიმიტში: a = sqrt(a0 * g_N)
    g_mond = sp.sqrt(a0 * G * M / r**2)
    
    # 2. RFG-ის ანიზოტროპიული აჩქარება: a = 2*Delta_p / (r*rho_solid)
    g_rfg = 2 * Delta_p / (r * rho_solid)
    
    # 3. ვუტოლებთ და ვხსნით Delta_p-სთვის
    dp_eq = sp.Eq(g_rfg, g_mond)
    dp_sol = sp.solve(dp_eq, Delta_p)[0]
    
    return {
        "MOND_deep_acceleration": g_mond,
        "RFG_anisotropic_acceleration": g_rfg,
        "required_vortex_stress_Delta_p": sp.simplify(dp_sol),
        "physical_conclusion": "MOND-ის რეჟიმი ითხოვს, რომ ვორტექსის ანიზოტროპია (Delta_p) შორ მანძილებზე "
                               "იყოს მუდმივი (r-ზე დამოუკიდებელი) და პროპორციული sqrt(M)-ის."
    }

def rfg_vortex_stress_generation():
    """
    იკვლევს, რა ტიპის დინამიკური პროფილი ქმნის
    მუდმივ Delta p ანიზოტროპიას შორ მანძილებზე.
    """
    r = sp.Symbol('r', real=True, positive=True)
    v_phi = sp.Function('v_phi')(r)
    eta_eff = sp.Symbol('eta_eff', real=True, positive=True) # მედიუმის ეფექტური ელასტიური სიხისტე
    C_dp = sp.Symbol('Delta_p_const', real=True)
    
    # დიფერენციალური ბრუნვის (Shear) სტრესი
    shear_stress = eta_eff * r * sp.diff(v_phi / r, r)
    eq = sp.Eq(shear_stress, C_dp)
    
    v_phi_sol = sp.dsolve(eq, v_phi).rhs
    
    return {
        "shear_stress_eq": eq,
        "v_phi_solution": v_phi_sol,
        "meaning": "მუდმივი ანიზოტროპია მოითხოვს მედიუმის დიფერენციალურ ბრუნვას (v_phi ∝ r ln r). "
                   "ეს მიუთითებს აქტიურ ვორტექსულ რეჟიმზე, მაგრამ გლობალურად საჭიროებს "
                   "core/cutoff საზღვრებს, რადგან მუდმივი სტრესი უსასრულო არეალში "
                   "ფიზიკურად დაუშვებელია."
    }


def rfg_finite_vortex_plateau_model():
    """
    finite-radius ვორტექსის მინიმალური phenomenological პროფილი.

    მიზანი არ არის ახალი fit-ის ჩასმა, არამედ იმ ფიზიკური პრობლემის გამოსწორება,
    რომ მუდმივი Delta p უსასრულობამდე ვერ გაგრძელდება. Plateau შიგნით იძლევა
    ბრტყელ მრუდს, გარეთ კი სტრესი ქრება.
    """
    r, R_v, rho_solid, Delta_p0 = sp.symbols(
        'r R_v rho_solid Delta_p0',
        real=True,
        positive=True,
    )

    delta_p = Delta_p0 * R_v**2 / (r**2 + R_v**2)
    a_extra = sp.simplify(2 * delta_p / (r * rho_solid))
    v_extra_sq = sp.simplify(a_extra * r)

    inner_v2 = sp.simplify(sp.limit(v_extra_sq, r, 0, dir='+'))
    outer_v2 = sp.simplify(sp.limit(r**2 * v_extra_sq, r, sp.oo))

    return {
        "Delta_p_profile": delta_p,
        "a_extra": a_extra,
        "v_extra_squared": v_extra_sq,
        "inner_plateau_limit": sp.Eq(sp.Symbol('v_extra^2(r<<R_v)'), inner_v2),
        "outer_falloff_limit": sp.Eq(sp.Symbol('r^2*v_extra^2(r>>R_v)'), outer_v2),
        "interpretation": "r << R_v-ზე ბრტყელი plateau; r >> R_v-ზე finite-energy cutoff და v ∝ 1/r.",
        "status": "minimal finite-vortex regularization; full profile must come from RFG vortex equations.",
    }


def rfg_vortex_cosmological_evolution():
    """
    ამოწმებს მომხმარებლის ინტუიციას: როგორც სინათლე ატარებს წარსულის 
    (უფრო სწრაფი დროის/მკვრივი მედიუმის) კვალს წითელი წანაცვლების სახით, 
    ასევე გალაქტიკური ვორტექსიც მაღალ რედშიფტზე (z) უნდა ავლენდეს იმ ეპოქის 
    მედიუმის მახასიათებლებს.
    """
    z, H_z, c_sym = sp.symbols('z H_z c', real=True, positive=True)
    
    # a_0(z) დამოკიდებულება კოსმოლოგიურ ეპოქაზე
    a0_z = c_sym * H_z / (2 * sp.pi)
    
    return {
        "hypothesis": "გალაქტიკური ვორტექსი ატარებს წარსული ეპოქის მედიუმის დინამიკის კვალს (როგორც სინათლის წითელი წანაცვლება).",
        "a0_evolution_eq": sp.Eq(sp.Symbol('a_0(z)'), a0_z),
        "physical_meaning": "ვინაიდან წარსულში ფუძე-მედიუმი უფრო აქტიური იყო და H(z) დიდი იყო, ვორტექსის დაჭიმულობაც მეტი უნდა ყოფილიყო.",
        "BTFR_shift_prediction": "BTFR: v^4 = G * M_b * a_0(z). ვინაიდან a_0(z) იზრდება, მაღალ რედშიფტზე (z>1) იგივე მასის გალაქტიკები უფრო სწრაფად უნდა ბრუნავდნენ.",
        "observational_test": "JWST-ის ახალი მონაცემები z=1-3 გალაქტიკებისთვის წარმოადგენს ამ ვარაუდის პირდაპირ ტესტს."
    }

def chi_square_curve(
    curve: RotationCurve,
    a0_m_s2: float,
    chi_coupling: float,
    ml_disk: float,
    ml_bulge: float,
    error_floor_km_s: float = 3.0,
) -> tuple[float, float]:
    model = predict_velocity_km_s(curve, a0_m_s2, chi_coupling, ml_disk, ml_bulge)
    chi2 = 0.0
    residual2_sum = 0.0

    for obs, err, pred in zip(curve.v_obs_km_s, curve.v_err_km_s, model):
        sigma = max(abs(err), error_floor_km_s)
        residual = obs - pred
        chi2 += (residual / sigma) ** 2
        residual2_sum += residual * residual

    rms = math.sqrt(residual2_sum / curve.n_points)
    return chi2, rms


def fit_curve_nuisance(
    curve: RotationCurve,
    a0_m_s2: float,
    chi_coupling: float,
    ml_disk_grid: list[float] | None = None,
    ml_bulge_grid: list[float] | None = None,
) -> GalaxyFit:
    """Fit disk/bulge mass-to-light ratios for one galaxy."""
    if ml_disk_grid is None:
        ml_disk_grid = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
    if ml_bulge_grid is None:
        ml_bulge_grid = [0.0, 0.3, 0.5, 0.7, 1.0]

    bulge_grid = ml_bulge_grid if curve.has_bulge else [0.0]
    best: GalaxyFit | None = None
    nuisance_count = 1 + (1 if curve.has_bulge else 0)
    dof_local = max(curve.n_points - nuisance_count, 1)

    for ml_disk in ml_disk_grid:
        for ml_bulge in bulge_grid:
            chi2, rms = chi_square_curve(
                curve,
                a0_m_s2=a0_m_s2,
                chi_coupling=chi_coupling,
                ml_disk=ml_disk,
                ml_bulge=ml_bulge,
            )
            candidate = GalaxyFit(
                name=curve.name,
                chi2=chi2,
                dof_local=dof_local,
                ml_disk=ml_disk,
                ml_bulge=ml_bulge,
                nuisance_count=nuisance_count,
                n_points=curve.n_points,
                rms_residual_km_s=rms,
            )
            if best is None or candidate.chi2 < best.chi2:
                best = candidate

    if best is None:
        raise RuntimeError(f"no fit candidate for {curve.name}")
    return best


def global_grid_fit(
    curves: list[RotationCurve],
    a0_grid: list[float] | None = None,
    chi_grid: list[float] | None = None,
    ml_disk_grid: list[float] | None = None,
    ml_bulge_grid: list[float] | None = None,
    status: str = "EMPIRICAL_LOCAL_DATA",
) -> GlobalFit:
    """Fit universal a0 and chi_coupling, with per-galaxy M/L nuisance grids."""
    if not curves:
        raise ValueError("global_grid_fit requires at least one rotation curve")

    if a0_grid is None:
        a0_grid = [0.8e-10, 1.0e-10, 1.2e-10, 1.4e-10, 1.6e-10]
    if chi_grid is None:
        chi_grid = [0.0, 0.5, 1.0, 1.5]

    best_global: GlobalFit | None = None

    for a0_m_s2 in a0_grid:
        for chi_coupling in chi_grid:
            galaxy_fits: list[GalaxyFit] = []
            total_chi2 = 0.0
            total_points = 0
            nuisance_params = 0

            for curve in curves:
                fit = fit_curve_nuisance(
                    curve,
                    a0_m_s2=a0_m_s2,
                    chi_coupling=chi_coupling,
                    ml_disk_grid=ml_disk_grid,
                    ml_bulge_grid=ml_bulge_grid,
                )
                galaxy_fits.append(fit)
                total_chi2 += fit.chi2
                total_points += fit.n_points
                nuisance_params += fit.nuisance_count

            dof = max(total_points - nuisance_params - 2, 1)
            candidate = GlobalFit(
                a0_m_s2=a0_m_s2,
                chi_coupling=chi_coupling,
                total_chi2=total_chi2,
                dof=dof,
                galaxy_fits=galaxy_fits,
                n_galaxies=len(curves),
                n_points=total_points,
                status=status,
            )
            if best_global is None or candidate.total_chi2 < best_global.total_chi2:
                best_global = candidate

    if best_global is None:
        raise RuntimeError("global fit did not evaluate any grid point")
    return best_global


def make_synthetic_curve(
    name: str,
    disk_scale: float,
    gas_scale: float,
    bulge_scale: float,
    ml_disk_true: float,
    ml_bulge_true: float,
    a0_true: float = A0_CANONICAL,
    chi_true: float = 1.0,
) -> RotationCurve:
    """Create deterministic mock data for code validation only."""
    radii = [0.7 + 0.55 * i for i in range(1, 31)]
    v_gas = [gas_scale * (1.0 - math.exp(-r / 7.0)) for r in radii]
    v_disk = [
        disk_scale * (r / 2.5) / (1.0 + (r / 2.5) ** 2) ** 0.65
        for r in radii
    ]
    v_bulge = [bulge_scale * math.exp(-r / 2.0) for r in radii]

    template = RotationCurve(
        name=name,
        radius_kpc=radii,
        v_obs_km_s=[1.0 for _ in radii],
        v_err_km_s=[4.0 for _ in radii],
        v_gas_km_s=v_gas,
        v_disk_km_s=v_disk,
        v_bulge_km_s=v_bulge,
    )
    v_obs = predict_velocity_km_s(
        template,
        a0_m_s2=a0_true,
        chi_coupling=chi_true,
        ml_disk=ml_disk_true,
        ml_bulge=ml_bulge_true,
    )
    return RotationCurve(
        name=name,
        radius_kpc=radii,
        v_obs_km_s=v_obs,
        v_err_km_s=[4.0 for _ in radii],
        v_gas_km_s=v_gas,
        v_disk_km_s=v_disk,
        v_bulge_km_s=v_bulge,
    )


def synthetic_smoke_test() -> GlobalFit:
    """Verify that the fitter recovers a known injected AQUAL model."""
    curves = [
        make_synthetic_curve("mock_LSB", 80.0, 35.0, 0.0, 0.5, 0.0),
        make_synthetic_curve("mock_HSB", 150.0, 45.0, 70.0, 0.6, 0.7),
        make_synthetic_curve("mock_dwarf", 45.0, 28.0, 0.0, 0.4, 0.0),
    ]
    return global_grid_fit(
        curves,
        a0_grid=[0.8e-10, 1.0e-10, 1.2e-10, 1.4e-10],
        chi_grid=[0.0, 0.5, 1.0, 1.5],
        ml_disk_grid=[0.3, 0.4, 0.5, 0.6, 0.8],
        ml_bulge_grid=[0.0, 0.5, 0.7, 1.0],
        status="SYNTHETIC_SMOKE_TEST_NOT_SPARC",
    )


def btfr_deep_mond_mass_solar(v_flat_km_s: float, a0_m_s2: float = A0_CANONICAL) -> float:
    """Deep-MOND BTFR: M_b = v^4/(G a0), returned in solar masses."""
    g_si = 6.67430e-11
    m_sun = 1.98847e30
    v_m_s = v_flat_km_s * 1000.0
    return v_m_s**4 / (g_si * a0_m_s2) / m_sun


def format_global_fit(fit: GlobalFit, max_galaxies: int = 8) -> list[str]:
    lines = [
        f"status: {fit.status}",
        f"galaxies: {fit.n_galaxies}, points: {fit.n_points}, dof: {fit.dof}",
        f"best a0: {fit.a0_m_s2:.3e} m/s^2",
        f"best chi_coupling: {fit.chi_coupling:.3g}",
        f"total chi2/dof: {fit.total_chi2:.3f}/{fit.dof} = {fit.chi2_dof:.3f}",
    ]

    ranked = sorted(fit.galaxy_fits, key=lambda item: item.chi2_dof_local, reverse=True)
    for galaxy in ranked[:max_galaxies]:
        lines.append(
            f"{galaxy.name}: chi2/dof={galaxy.chi2_dof_local:.3f}, "
            f"M/Ld={galaxy.ml_disk:.2f}, M/Lb={galaxy.ml_bulge:.2f}, "
            f"rms={galaxy.rms_residual_km_s:.2f} km/s"
        )
    return lines


def model_scope_notes() -> list[str]:
    return [
        "RFG extracts the supersolid anisotropy (Delta p) required by a given rotation curve.",
        "AQUAL is kept as an empirical baseline to compare RFG's Delta p profile against.",
        "Static 1/r^2 anisotropy is not enough; observed flat regions require a finite vortex plateau.",
        "A real SPARC verdict requires the local 175-galaxy rotmod dataset and distance/inclination priors.",
        "Without local SPARC files, only the synthetic smoke-test is executed.",
    ]


def main() -> None:
    print("=" * 72)
    print("PHASE 33: SPARC rotation curves — RFG/AQUAL chi-square pipeline")
    print("=" * 72)

    print("\n1. SPARC scope")
    for key, value in SPARC_SUMMARY.items():
        print(f"  {key:18s}: {value}")

    print("\n2. Local data discovery")
    curves, warnings = load_sparc_dataset()
    print(f"  loaded curves: {len(curves)}")
    if warnings:
        print(f"  parse warnings: {len(warnings)}")
        for warning in warnings[:5]:
            print(f"    {warning}")

    if curves:
        print("\n3. Empirical local-data fit")
        fit = global_grid_fit(curves)
        for line in format_global_fit(fit):
            print(f"  {line}")
    else:
        print("  local SPARC rotmod files were not found.")
        print("  searched:")
        for directory in DEFAULT_DATA_DIRS:
            print(f"    {directory}")

        print("\n3a. RFG მიკროფიზიკა — ანიზოტროპიის ექსტრაქცია სინთეტიკური მონაცემებიდან")
        mock_curve = make_synthetic_curve("mock_HSB", 150.0, 45.0, 70.0, 0.6, 0.7)
        dp_profile = rfg_required_anisotropy_delta_p(mock_curve, 0.6, 0.7)
        print(f"  გალაქტიკა {mock_curve.name}-ისთვის ვაკუუმის სტრესის (Delta p) პროფილი:")
        for r, dp in zip(mock_curve.radius_kpc[::5], dp_profile[::5]):
            print(f"    r = {r:4.1f} kpc -> Delta p = {dp:.3e} Pa")
        print("  დასკვნა: RFG ახლა უკვე ზუსტად ითვლის, რა ფიზიკური დაჭიმულობა სჭირდება")
        print("  მედიუმს ბნელი მატერიის ჩასანაცვლებლად. ეს აღარ არის უბრალოდ MOND ფიტი.")

        print("\n3b. RFG ანიზოტროპიის თეორიული სკალირების აუდიტი")
        scaling_audit = rfg_theoretical_anisotropy_scaling()
        print(f"  MOND-ის მოთხოვნა (v=const): Delta_p ∝ {scaling_audit['required_for_flat_curves']}")
        print(f"  სტატიკური ლაგრანჟიანის (c_I1) პროგნოზი: Delta_p ∝ {scaling_audit['lagrangian_static_prediction']}")
        print(f"  სტატიკური აჩქარება: {scaling_audit['extra_acceleration_from_static']}")
        print(f"  სიჩქარის სკალირება: {scaling_audit['velocity_scaling_from_static']}")
        print(f"  კონფლიქტი: {scaling_audit['mathematical_conflict']}")
        print(f"  გადაჭრა: {scaling_audit['physical_resolution']}")

        print("\n3c. a_0 პარამეტრის და ვორტექსის სტრესის კავშირი")
        vortex = rfg_a0_vortex_emergence()
        for k, v in vortex.items():
            print(f"  {k:32s}: {v}")

        print("\n3d. მუდმივი ანიზოტროპიის წარმომქმნელი მექანიზმი (მედიუმის ბრუნვა)")
        stress_gen = rfg_vortex_stress_generation()
        for k, v in stress_gen.items():
            print(f"  {k:20s}: {v}")

        print("\n3e. finite-radius vortex plateau regularization")
        finite_vortex = rfg_finite_vortex_plateau_model()
        for k, v in finite_vortex.items():
            print(f"  {k:24s}: {v}")

        print("\n3f. კოსმოლოგიური ევოლუცია და ვორტექსის მეხსიერება (z-დამოკიდებულება)")
        z_evol = rfg_vortex_cosmological_evolution()
        for k, v in z_evol.items():
            print(f"  {k:25s}: {v}")

        print("\n3g. Synthetic smoke-test (pipeline only, not SPARC)")
        fit = synthetic_smoke_test()
        for line in format_global_fit(fit):
            print(f"  {line}")
        recovered = (
            abs(fit.a0_m_s2 - A0_CANONICAL) < 1e-20
            and abs(fit.chi_coupling - 1.0) < 1e-12
        )
        print(f"  recovery status: {'OK' if recovered else 'CHECK'}")

    print("\n4. BTFR deep-MOND scale check")
    for v_flat in (50.0, 100.0, 200.0):
        mass = btfr_deep_mond_mass_solar(v_flat)
        print(f"  v_flat={v_flat:5.1f} km/s -> M_b={mass:.3e} M_sun")

    print("\n5. Scope notes")
    for note in model_scope_notes():
        print(f"  - {note}")


if __name__ == "__main__":
    main()
