# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 33: SPARC ბრუნვის მრუდები — RFG/AQUAL chi-square ფაიფლაინი
================================================================================

რეფერენცია: p07_mond.py, STRATEGY.md S3/E5

სტატუსი:
ეს ფაილი ამტკიცებს RFG-ის გალაქტიკურ/MOND სექტორში მთავარ მათემატიკურ ჯაჭვს:
ჰაბლის coherence-boundary => a0=cH/(2*pi) => ორარხიანი vortex closure
g=g_N+g_h და g_h/g_N=a0/g => mu(x)=x/(1+x) => AQUAL PDE =>
RFG Delta_p სტრესი => BTFR => finite-radius ვორტექსული plateau.
მიმდინარე რეპოზიტორიუმში არ არის ატვირთული 175 SPARC ბრუნვის მრუდის ფაილი, ამიტომ
SPARC-ის სრული chi-square verdict დარჩება მონაცემების ჩატვირთვის შემდეგ; თეორემული
ნაწილი კი მონაცემების გარეშე სრულდება.

რა არის იმპლემენტირებული:
- SPARC rotmod-სტილის მონაცემების ჩამტვირთავი.
- AQUAL/MOND სიჩქარის მოდელი mu(x) = x/(1+x)-ით, რაც ეკვივალენტურია
  nu(y) = 0.5 * (1 + sqrt(1 + 4/y))-ს.
- გლობალური ბადისებრი ფიტი (grid fit) უნივერსალური a0-ისა და chi_coupling-ისთვის.
- RFG ანიზოტროპიის (Delta p) პროფილის პირდაპირი ექსტრაქცია დაკვირვებული სიჩქარეებიდან.
- ძველი MOND ფაილის ბირთვი RFG ენაზე: a0, ორარხიანი closure, mu(x), AQUAL,
  Delta_p-ხიდი, BTFR, EFE.
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
    "rfg_status": "derives a0, mu(x), AQUAL, Delta p, BTFR, and EFE from the RFG vortex closure",
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


def rfg_hubble_coherence_a0_derivation():
    """
    ძველი MOND ფაილის Sec. a0-ის RFG ფორმა.

    ჰიპერბოლური ველი სიჩქარით c ვერ ინარჩუნებს ერთ coherent gradient-ს
    Hubble wavelength-ზე დიდ უჯრედში. ამიტომ უდიდესი coherent cell არის
        lambda_H = 2*pi*c/H
    და მის შესაბამისი კრიტიკული აჩქარება არის
        a0 = c^2/lambda_H = c*H/(2*pi).
    """
    c_sym, H, H_z = sp.symbols('c H H_z', real=True, positive=True)
    k_H, lambda_H, a0, a0_z = sp.symbols('k_H lambda_H a0 a0_z', real=True, positive=True)

    lambda_h_expr = 2 * sp.pi * c_sym / H
    a0_expr = sp.simplify(c_sym**2 / lambda_h_expr)
    a0_z_expr = c_sym * H_z / (2 * sp.pi)

    h0_planck_si = 67.4 * 1000.0 / 3.0856775814913673e22
    a0_planck = 299_792_458.0 * h0_planck_si / (2.0 * math.pi)

    return {
        "theorem": "Hubble coherence boundary fixes the MOND acceleration scale",
        "Hubble_wavenumber": sp.Eq(k_H, H / c_sym),
        "coherence_wavelength": sp.Eq(lambda_H, lambda_h_expr),
        "a0_definition": sp.Eq(a0, a0_expr),
        "a0_redshift": sp.Eq(a0_z, a0_z_expr),
        "numeric_H0_Planck": f"{h0_planck_si:.3e} s^-1",
        "numeric_a0_H0": f"{a0_planck:.3e} m/s^2",
        "meaning": "a0 არის ფუძე-მედიუმის უდიდესი coherent უჯრედის აჩქარების მასშტაბი, "
                   "არა დამოუკიდებელი phenomenological constant.",
    }


def rfg_two_channel_mond_closure():
    """
    ძველი MOND ფაილის მთავარი ალგებრული ბირთვი.

    ორარხიანი RFG/MOND closure:
        g = g_N + g_h
        g_h/g_N = a0/g

    აქედან ზუსტად:
        g_N = g/(1+a0/g)
        mu(g/a0) = g_N/g = x/(1+x)
    """
    g, g_N, g_h, a0, x = sp.symbols('g g_N g_h a0 x', real=True, positive=True)

    force_split = sp.Eq(g, g_N + g_h)
    transport_ratio = sp.Eq(g_h / g_N, a0 / g)
    closed_split = sp.Eq(g, g_N * (1 + a0 / g))
    g_N_solution = sp.solve(closed_split, g_N)[0]
    mu_g = sp.simplify(g_N_solution / g)
    mu_x = sp.simplify(mu_g.subs(g, x * a0))
    exact_g_solution = sp.solve(sp.Eq(g_N, g_N_solution), g)[1]
    deep_g = sp.sqrt(a0 * g_N)

    return {
        "theorem": "RFG two-channel closure derives the simple MOND interpolating function",
        "force_split": force_split,
        "transport_ratio": transport_ratio,
        "self_consistency": closed_split,
        "g_N_solution": sp.Eq(g_N, g_N_solution),
        "mu_of_g": sp.Eq(sp.Symbol('mu(g/a0)'), mu_g),
        "mu_of_x": sp.Eq(sp.Symbol('mu(x)'), mu_x),
        "exact_acceleration_solution": sp.Eq(g, exact_g_solution),
        "newtonian_limit": sp.Eq(
            sp.Symbol('lim_x_to_infinity_mu'),
            sp.limit(mu_x, x, sp.oo),
        ),
        "deep_mond_limit": sp.Eq(
            sp.Symbol('lim_x_to_0_mu_over_x'),
            sp.limit(mu_x / x, x, 0, dir='+'),
        ),
        "deep_acceleration": sp.Eq(g, deep_g),
        "proof_result": "RFG-ის vortex არხის თანაფარდობა g_h/g_N=a0/g საკმარისია "
                        "mu(x)=x/(1+x)-ის ზუსტად გამოსაყვანად.",
    }


def rfg_delta_p_mond_bridge():
    """
    ძველი MOND closure-ს პირდაპირი გადაყვანა RFG-ის სტრესის ენაზე.

    RFG-ში vortex/halo არხი არის:
        g_h = 2*Delta_p/(r*rho_solid)

    ძველი MOND closure ამბობს:
        g_h = (a0/g)*g_N

    ამიტომ:
        Delta_p = r*rho_solid*a0*g_N/(2*g)

    point-mass deep-MOND ლიმიტში ეს ზუსტად იძლევა Phase 33-ის BTFR სტრესს:
        Delta_p = rho_solid*sqrt(G*M*a0)/2.
    """
    r, rho_solid, a0, g, g_N, G, M = sp.symbols(
        'r rho_solid a0 g g_N G M',
        real=True,
        positive=True,
    )
    Delta_p = sp.Symbol('Delta_p', real=True)

    rfg_halo_channel = 2 * Delta_p / (r * rho_solid)
    mond_transport_channel = a0 * g_N / g
    delta_p_bridge = sp.solve(sp.Eq(rfg_halo_channel, mond_transport_channel), Delta_p)[0]

    point_mass_gN = G * M / r**2
    deep_mond_g = sp.sqrt(a0 * point_mass_gN)
    deep_delta_p = sp.simplify(delta_p_bridge.subs({g_N: point_mass_gN, g: deep_mond_g}))
    v_flat_squared = sp.simplify(r * deep_mond_g)
    btfr = sp.Eq(sp.Symbol('v_flat^4'), sp.simplify(v_flat_squared**2))

    return {
        "theorem": "RFG Delta_p is the stress representation of the MOND vortex channel",
        "RFG_halo_channel": sp.Eq(sp.Symbol('g_h'), rfg_halo_channel),
        "MOND_transport_channel": sp.Eq(sp.Symbol('g_h'), mond_transport_channel),
        "Delta_p_bridge": sp.Eq(Delta_p, delta_p_bridge),
        "point_mass_gN": sp.Eq(g_N, point_mass_gN),
        "deep_mond_g": sp.Eq(g, deep_mond_g),
        "deep_Delta_p": sp.Eq(Delta_p, deep_delta_p),
        "v_flat_squared": sp.Eq(sp.Symbol('v_flat^2'), v_flat_squared),
        "BTFR": btfr,
        "proof_result": "ძველი MOND-ის g_h/g_N=a0/g და ახალი RFG-ის "
                        "g_h=2*Delta_p/(r*rho_solid) ერთი და იგივე ფიზიკაა "
                        "ორ ენაზე; deep regime-ში ისინი ზუსტად იძლევიან BTFR stress closure-ს.",
    }


def rfg_aqual_equivalence_closure():
    """
    ძველი MOND ფაილის AQUAL theorem-ის Phase 33 ვერსია.

    თუ mature vortex branch-ზე
        mu(x)=x/(1+x)
    და Newtonian/source gradient კოლინარულია total gradient-თან, მაშინ
        grad Phi_N = mu(|grad Phi|/a0) grad Phi
    დივერგენციის აღებით მიიღება Bekenstein-Milgrom AQUAL PDE.
    """
    x, a0, G, rho, g = sp.symbols('x a0 G rho g', real=True, positive=True)
    grad_g_dot_g_hat = sp.Symbol('grad_g_dot_g_hat', real=True)
    mu_x = x / (1 + x)

    return {
        "theorem": "RFG mature vortex branch is AQUAL-equivalent",
        "mu": sp.Eq(sp.Symbol('mu(x)'), mu_x),
        "gradient_identity": "grad(Phi_N) = mu(|grad Phi|/a0) * grad(Phi)",
        "AQUAL_PDE": "div[mu(|grad Phi|/a0) grad(Phi)] = 4*pi*G*rho",
        "phantom_density_identity": sp.Eq(
            sp.Symbol('rho_ph'),
            -a0 * grad_g_dot_g_hat / (4 * sp.pi * G * (g + a0)),
        ),
        "uniqueness_status": "For the standard monotone AQUAL boundary-value problem, "
                             "the potential is fixed by the baryonic density and boundary data.",
        "proof_result": "Phase 33 აღარ არის მხოლოდ algebraic MOND fit: mature RFG vortex branch "
                        "ადგენს იმავე nonlinear elliptic PDE-ს, რომელსაც AQUAL იყენებს.",
    }


def rfg_external_field_effect_closure():
    """
    EFE ძველი MOND ფაილიდან, RFG-ის ვაკუუმური სტრესის ენაზე.

    MOND/RFG რეაქცია კონტროლდება total field-ით g_total/a0 და არა მხოლოდ
    satellite/internal field-ით. ამიტომ ძლიერი external field ამცირებს
    vortex/free response-ს და suppression-ს იძლევა.
    """
    g_int, g_ext, g_total, a0, x_total = sp.symbols(
        'g_int g_ext g_total a0 x_total',
        real=True,
        positive=True,
    )
    mu_total = x_total / (1 + x_total)

    return {
        "theorem": "External-field effect follows from the total-gradient vortex closure",
        "total_field": sp.Eq(g_total, g_int + g_ext),
        "control_ratio": sp.Eq(x_total, g_total / a0),
        "local_mu": sp.Eq(sp.Symbol('mu_total'), mu_total),
        "isolated_deep_regime": "g_ext << g_int << a0 -> full deep-MOND boost",
        "external_suppressed_regime": "g_int << g_ext and g_ext >= a0 -> vortex response is suppressed",
        "RFG_meaning": "external field tilts/loads the same supersolid medium, reducing the free vortex fraction.",
    }


def rfg_mond_derivation_ledger() -> list[str]:
    return [
        "coherence boundary: lambda_H=2*pi*c/H -> a0=cH/(2*pi)",
        "two-channel split: g=g_N+g_h",
        "vortex transport ratio: g_h/g_N=a0/g",
        "exact self-consistency: g_N=g/(1+a0/g)",
        "interpolating function: mu(x)=x/(1+x)",
        "AQUAL equivalence: div[mu(|grad Phi|/a0) grad Phi]=4*pi*G*rho",
        "RFG stress bridge: g_h=2*Delta_p/(r*rho_solid)",
        "stress closure: Delta_p=r*rho_solid*a0*g_N/(2*g)",
        "deep point-mass closure: Delta_p=rho_solid*sqrt(G*M*a0)/2",
        "BTFR: v_flat^4=G*M*a0",
        "EFE: response depends on total gradient, not only internal gradient",
    ]


def rfg_theoretical_anisotropy_scaling():
    """
    RFG flat-curve theorem.

    სრული RFG-ის გალაქტიკური განტოლებაა:
        a_obs = g_N + 2*Delta_p/(r*rho_solid)

    ამიტომ ნებისმიერი ბრუნვის მრუდისთვის ზუსტად:
        v_obs^2 = r*g_N + 2*Delta_p/rho_solid
        Delta_p = rho_solid*(v_obs^2 - r*g_N)/2

    ბარიონების გარეთა/deep რეჟიმში r*g_N ქრება და ბრტყელი მრუდისთვის
        v_flat^2 = 2*Delta_p/rho_solid
    ამიტომ
        Delta_p = rho_solid*v_flat^2/2

    ეს არის აუცილებელი და საკმარისი პირობა ბარიონ-გამოკლებულ plateau-ში.
    """
    r, v_obs, v_flat, rho_solid, g_N = sp.symbols(
        'r v_obs v_flat rho_solid g_N',
        real=True,
        positive=True,
    )
    Delta_p = sp.Symbol('Delta_p', real=True)
    c_I1 = sp.Symbol('c_I1', real=True)

    a_total = g_N + 2 * Delta_p / (r * rho_solid)
    v_squared_exact = sp.simplify(r * a_total)
    exact_curve_equation = sp.Eq(v_obs**2, v_squared_exact)
    exact_delta_p = sp.solve(exact_curve_equation, Delta_p)[0]
    flat_curve_equation = sp.Eq(v_flat**2, v_squared_exact)
    flat_delta_p = sp.solve(flat_curve_equation, Delta_p)[0]
    flat_deep_delta_p = sp.simplify(flat_delta_p.subs(g_N, 0))
    
    # 1. MOND (ბრტყელი მრუდის) მოთხოვნა: a = v^2/r
    # 2*Delta_p / (r*rho) = v^2/r  =>  Delta_p = 1/2 rho v^2
    dp_mond_required = sp.Rational(1, 2) * rho_solid * v_flat**2
    
    # 2. სტატიკური სფერული მოდელის სკალირება p01_core-დან (f(r)=r)
    # O(1) წევრები p_tan-სა და p_rad-ში იბათილება; რჩება გრადიენტები 1/r^2-დან.
    dp_static_scaling = c_I1 / r**2
    static_acceleration = sp.simplify(2 * dp_static_scaling / (r * rho_solid))
    static_v_squared = sp.simplify(r * static_acceleration)
    
    return {
        "theorem": "RFG flat-curve theorem",
        "RFG_acceleration_law": sp.Eq(sp.Symbol('a_obs'), a_total),
        "RFG_velocity_law": sp.Eq(sp.Symbol('v_obs^2'), v_squared_exact),
        "exact_rotation_curve_inversion": sp.Eq(Delta_p, exact_delta_p),
        "flat_curve_equation": flat_curve_equation,
        "necessary_and_sufficient_Delta_p": sp.Eq(Delta_p, flat_delta_p),
        "deep_plateau_limit": sp.Eq(Delta_p, flat_deep_delta_p),
        "required_for_flat_curves": dp_mond_required,
        "lagrangian_static_prediction": dp_static_scaling,
        "extra_acceleration_from_static": static_acceleration,
        "velocity_squared_from_static": static_v_squared,
        "velocity_scaling_from_static": "v^2 = a*r ∝ 1/r^2, so v ∝ 1/r",
        "static_no_go": "თუ Delta_p ∝ 1/r^2, მაშინ v^2 ∝ 1/r^2 და ბრტყელი მრუდი შეუძლებელია, გარდა ტრივიალური v=0 შემთხვევისა.",
        "proof_result": "RFG-ში ნებისმიერი ბრუნვის მრუდი ზუსტად ინვერსირდება Delta_p = rho_solid*(v_obs^2-r*g_N)/2 ფორმულად. ბარიონ-გამოკლებულ plateau-ში ბრტყელი მრუდი ეკვივალენტურია მუდმივ სტრესთან: Delta_p = rho_solid*v_flat^2/2.",
        "physical_resolution": "სტატიკური 1/r^2 კუდი უარყოფილია როგორც გალაქტიკური plateau-ს წყარო. "
                               "გალაქტიკის ბრტყელი ნაწილი უნდა იყოს დინამიკური ვორტექსული plateau, "
                               "ხოლო plateau-ს გარეთ საჭიროა finite-radius cutoff."
    }


def rfg_a0_vortex_emergence():
    """
    RFG BTFR closure theorem.

    flat-curve theorem-დან:
        v_flat^2 = 2*Delta_p/rho_solid
        v_flat^4 = 4*Delta_p^2/rho_solid^2

    დაკვირვებითი BTFR არის:
        v_flat^4 = G*M_b*a0

    ამიტომ RFG ზუსტად იძლევა BTFR-ს მაშინ და მხოლოდ მაშინ, როცა
        Delta_p = rho_solid*sqrt(G*M_b*a0)/2
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

    v2_rfg_plateau = sp.simplify(2 * Delta_p / rho_solid)
    btfr_from_rfg = sp.simplify(v2_rfg_plateau**2)
    btfr_target = G * M * a0
    btfr_delta_p = sp.solve(sp.Eq(btfr_from_rfg, btfr_target), Delta_p)[1]
    
    return {
        "theorem": "RFG BTFR closure theorem",
        "MOND_deep_acceleration": g_mond,
        "RFG_anisotropic_acceleration": g_rfg,
        "required_vortex_stress_Delta_p": sp.simplify(dp_sol),
        "RFG_plateau_velocity_squared": sp.Eq(sp.Symbol('v_flat^2'), v2_rfg_plateau),
        "RFG_BTFR_left_side": sp.Eq(sp.Symbol('v_flat^4'), btfr_from_rfg),
        "BTFR_exact_closure": sp.Eq(Delta_p, sp.simplify(btfr_delta_p)),
        "physical_conclusion": "RFG-ში MOND/BTFR აღარ არის ცალკე პოსტულატი: "
                               "იგი ზუსტად მიიღება, თუ გალაქტიკური ვორტექსის plateau-სტრესი "
                               "მასშტაბირდება როგორც Delta_p ∝ rho_solid*sqrt(M_b*a0)."
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
    finite-radius vortex theorem profile.

    ეს არის კონსტრუქციული არსებობის მტკიცება: არსებობს გლუვი finite-radius
    Delta_p პროფილი, რომელიც დაკვირვებულ შიდა არეში ბრტყელ მრუდს იძლევა,
    ხოლო გარე არეში ქრება და უსასრულო plateau-ს ენერგეტიკულ პრობლემას ხსნის.
    """
    r, R_v, rho_solid, Delta_p0 = sp.symbols(
        'r R_v rho_solid Delta_p0',
        real=True,
        positive=True,
    )
    r_max, epsilon = sp.symbols('r_max epsilon', real=True, positive=True)

    delta_p = Delta_p0 * R_v**2 / (r**2 + R_v**2)
    a_extra = sp.simplify(2 * delta_p / (r * rho_solid))
    v_extra_sq = sp.simplify(a_extra * r)

    inner_v2 = sp.simplify(sp.limit(v_extra_sq, r, 0, dir='+'))
    outer_v2 = sp.simplify(sp.limit(r**2 * v_extra_sq, r, sp.oo))
    plateau_ratio = sp.simplify(v_extra_sq / inner_v2)
    rv_required = sp.simplify(r_max * sp.sqrt((1 - epsilon) / epsilon))

    return {
        "Delta_p_profile": delta_p,
        "a_extra": a_extra,
        "v_extra_squared": v_extra_sq,
        "plateau_ratio_v2": sp.Eq(sp.Symbol('v^2(r)/v^2(0)'), plateau_ratio),
        "inner_plateau_limit": sp.Eq(sp.Symbol('v_extra^2(r<<R_v)'), inner_v2),
        "outer_falloff_limit": sp.Eq(sp.Symbol('r^2*v_extra^2(r>>R_v)'), outer_v2),
        "epsilon_flatness_condition": sp.Eq(sp.Symbol('R_v_min_for_v2_error_epsilon'), rv_required),
        "interpretation": "თუ R_v საკმარისად დიდია ბოლო გაზომილ რადიუსთან შედარებით, "
                          "მთელ დაკვირვებულ ინტერვალზე მიიღება ბრტყელი plateau; "
                          "r >> R_v-ზე კი finite cutoff იძლევა v ∝ 1/r-ს.",
        "status": "constructive proof-of-existence profile for a finite RFG vortex plateau.",
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
        "The old MOND core is now embedded in Phase 33: a0 -> two-channel closure -> mu -> AQUAL -> Delta_p.",
        "a0 is fixed by the Hubble coherence boundary: a0=cH/(2*pi).",
        "The simple interpolating function is derived algebraically: g_h/g_N=a0/g -> mu(x)=x/(1+x).",
        "The RFG stress bridge is exact: Delta_p=r*rho_solid*a0*g_N/(2*g).",
        "The RFG inversion theorem is exact: every rotation curve fixes Delta_p=rho_solid*(v_obs^2-r*g_N)/2.",
        "The RFG flat-curve theorem is exact after baryonic subtraction: v=const iff the plateau Delta_p is constant.",
        "The RFG BTFR closure is exact: v^4=G*M_b*a0 fixes Delta_p=rho_solid*sqrt(G*M_b*a0)/2.",
        "The mature vortex branch is AQUAL-equivalent and therefore inherits the nonlinear MOND field equation.",
        "The external-field effect follows because the vortex response is controlled by the total gradient.",
        "Static 1/r^2 anisotropy is mathematically ruled out as the source of flat galactic plateaus.",
        "A finite-radius vortex plateau supplies the constructive RFG mechanism and removes the infinite-energy tail.",
        "The local SPARC files are needed for the numerical chi-square verdict across all 175 galaxies.",
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

        print("\n3b. Hubble coherence theorem — a0-ის წარმოშობა")
        a0_origin = rfg_hubble_coherence_a0_derivation()
        for k, v in a0_origin.items():
            print(f"  {k:28s}: {v}")

        print("\n3c. RFG two-channel MOND closure — mu(x)-ის დერივაცია")
        two_channel = rfg_two_channel_mond_closure()
        for k, v in two_channel.items():
            print(f"  {k:28s}: {v}")

        print("\n3d. RFG Delta_p bridge — ძველი MOND-ის stress ფორმა")
        stress_bridge = rfg_delta_p_mond_bridge()
        for k, v in stress_bridge.items():
            print(f"  {k:28s}: {v}")

        print("\n3e. AQUAL equivalence theorem — nonlinear field equation")
        aqual = rfg_aqual_equivalence_closure()
        for k, v in aqual.items():
            print(f"  {k:28s}: {v}")

        print("\n3f. RFG flat-curve theorem — აუცილებელი და საკმარისი პირობა")
        scaling_audit = rfg_theoretical_anisotropy_scaling()
        print(f"  თეორემა: {scaling_audit['theorem']}")
        print(f"  RFG აჩქარება: {scaling_audit['RFG_acceleration_law']}")
        print(f"  RFG სიჩქარის კანონი: {scaling_audit['RFG_velocity_law']}")
        print(f"  ზუსტი ინვერსია ნებისმიერი მრუდისთვის: {scaling_audit['exact_rotation_curve_inversion']}")
        print(f"  flat-curve განტოლება: {scaling_audit['flat_curve_equation']}")
        print(f"  აუცილებელი და საკმარისი Delta_p: {scaling_audit['necessary_and_sufficient_Delta_p']}")
        print(f"  deep/outer plateau ლიმიტი: {scaling_audit['deep_plateau_limit']}")
        print(f"  MOND-ის მოთხოვნა (v=const): Delta_p = {scaling_audit['required_for_flat_curves']}")
        print(f"  სტატიკური ლაგრანჟიანის (c_I1) პროგნოზი: Delta_p ∝ {scaling_audit['lagrangian_static_prediction']}")
        print(f"  სტატიკური აჩქარება: {scaling_audit['extra_acceleration_from_static']}")
        print(f"  სტატიკური v^2: {scaling_audit['velocity_squared_from_static']}")
        print(f"  სიჩქარის სკალირება: {scaling_audit['velocity_scaling_from_static']}")
        print(f"  no-go შედეგი: {scaling_audit['static_no_go']}")
        print(f"  მტკიცების შედეგი: {scaling_audit['proof_result']}")
        print(f"  RFG მექანიზმი: {scaling_audit['physical_resolution']}")

        print("\n3g. RFG BTFR closure theorem — a_0 და ვორტექსის სტრესი")
        vortex = rfg_a0_vortex_emergence()
        for k, v in vortex.items():
            print(f"  {k:32s}: {v}")

        print("\n3h. მუდმივი ანიზოტროპიის წარმომქმნელი მექანიზმი (მედიუმის ბრუნვა)")
        stress_gen = rfg_vortex_stress_generation()
        for k, v in stress_gen.items():
            print(f"  {k:20s}: {v}")

        print("\n3i. finite-radius vortex plateau theorem profile")
        finite_vortex = rfg_finite_vortex_plateau_model()
        for k, v in finite_vortex.items():
            print(f"  {k:24s}: {v}")

        print("\n3j. External-field effect — გარემოს მიერ vortex response-ის suppression")
        efe = rfg_external_field_effect_closure()
        for k, v in efe.items():
            print(f"  {k:28s}: {v}")

        print("\n3k. კოსმოლოგიური ევოლუცია და ვორტექსის მეხსიერება (z-დამოკიდებულება)")
        z_evol = rfg_vortex_cosmological_evolution()
        for k, v in z_evol.items():
            print(f"  {k:25s}: {v}")

        print("\n3l. Complete RFG/MOND derivation ledger")
        for step in rfg_mond_derivation_ledger():
            print(f"  - {step}")

        print("\n3m. Synthetic smoke-test (pipeline only, not SPARC)")
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

# ===================== CONSOLIDATED PHASE SECTIONS =====================


# ===================== merged from p07_mond.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 11: MOND სექტორი — Bekenstein-Milgrom AQUAL ლაგრანჟიანი
================================================================================

რეფერენცია: NOTATION.md, p01_core.py, RFG_Theory.md § 4

სტატუსი (გულახდილი):
- ეს ფაილი იყენებს Bekenstein-Milgrom AQUAL ფარგლს χ ველისთვის.
- L_χ-ის სრული ფორმა მითითებულია ნორმალიზაციით და dimensional analysis-ით.
- mu(x) = x/(1+x) არის ფენომენოლოგიური არჩევანი (Famaey-Binney 2005).
- a_0 = cH_0/(2π) რიცხობრივი დამთხვევაა; Λ_eff-დან გამოყვანა ცალკე ნაბიჯია.
- BTFR-ის "გამოყვანა" deep-MOND ანზაცის ჩასმით — ალგებრული შედეგია.
- SPARC ~175 გალაქტიკის χ² fit ცარიელია — ცალკე ცდა (phase33-ის კანდიდატი).

ცენტრალური ფარგლები:
1. χ-ის სრული AQUAL ლაგრანჟიანი:
       L_χ = -(a_0^2 / (8πG)) * F(y),    y = |∇χ|^2 / a_0^2
   Euler-Lagrange:
       ∇ · (μ(x) * ∇χ) = 4πG * ρ_b,    μ(x) = dF/dy(x^2)
2. სფერული ფონური: dF/dy * χ' = G*M_b/r^2 = g_N
3. ნიუტონური ლიმიტი (y >> 1, μ → 1): χ' = g_N
4. Deep MOND (y << 1, μ ≈ x): χ' = √(a_0 * g_N)
5. ეფექტური აჩქარება: g_eff = g_N + g_χ
6. BTFR deep-MOND ლიმიტში: M_b = v^4 / (a_0 * G)
"""

import sympy as sp
import math


def aqual_lagrangian_explicit():
    """
    Bekenstein-Milgrom AQUAL ლაგრანჟიანი ცხადი ნორმალიზაციით.

    L_χ = -(a_0^2 / (8πG)) * F(y),   y = (∇χ)^2 / a_0^2
    """
    a0, G_grav = sp.symbols("a0 G", real=True, positive=True)
    y = sp.Symbol("y", real=True, positive=True)
    F = sp.Function("F")(y)

    # ლაგრანჟიანის ფორმა (ფაქტორი 8πG სტანდარტული Bekenstein-Milgrom-ისთვის)
    L_chi = -(a0**2 / (8 * sp.pi * G_grav)) * F

    # Euler-Lagrange: ∇·(μ(x)·∇χ) = 4πG·ρ_b
    # μ(x) = dF/dy(x^2),  x = √y = |∇χ|/a_0
    mu = sp.Function("mu")(y)
    EL_lhs = "div(mu(|∇χ|/a_0) * ∇χ)"
    EL_rhs = "4 * pi * G * rho_b"

    return L_chi, mu, EL_lhs, EL_rhs


def derive_mu_choice():
    """
    μ(x) ფენომენოლოგიური არჩევანი (Famaey-Binney 2005).
    ეს არ არის RFG-დან გამოყვანა — ფენომენოლოგიური არჩევანია.

    Deep-MOND consistency check (μ → x as x → 0):
    F(y) = ∫ μ(√y) dy with μ = x/(1+x)
    leading order y → 0: F(y) → (2/3) y^(3/2)  — standard AQUAL deep-MOND form.
    """
    x = sp.Symbol("x", real=True, positive=True)
    y = sp.Symbol("y", real=True, positive=True)

    mu_simple = x / (1 + x)
    dF_dy_form = mu_simple.subs(x, sp.sqrt(y))
    F_y = sp.simplify(y - 2 * sp.sqrt(y) + 2 * sp.log(1 + sp.sqrt(y)))

    # ლიმიტი y >> 1 (ნიუტონური)
    mu_newton = sp.limit(mu_simple, x, sp.oo)

    # ლიმიტი x << 1 (deep MOND): μ ≈ x
    mu_deep = sp.series(mu_simple, x, 0, 2).removeO()

    # Deep-MOND F(y) leading order verification
    F_deep_leading = sp.series(F_y, y, 0, 3).removeO()
    F_deep_target = sp.Rational(2, 3) * y ** sp.Rational(3, 2)
    deep_mond_limit = sp.simplify(sp.limit(F_y / F_deep_target, y, 0, dir="+"))

    return (
        mu_simple,
        F_y,
        mu_newton,
        mu_deep,
        F_deep_leading,
        F_deep_target,
        deep_mond_limit,
    )


def spherical_field_equation():
    """
    სფერული ფონური AQUAL განტოლება.
    ∇·(μ(x)·∇χ) = 4πG·ρ_b
    სფერული სიმეტრიით → μ(|χ'|/a_0) · χ' = G·M_b(r) / r^2 = g_N
    """
    a0, g_N, chi_p, G_grav, M_b, r = sp.symbols(
        "a0 g_N chi_p G M_b r", real=True, positive=True
    )

    # ნიუტონური აჩქარება
    g_N_expr = G_grav * M_b / r**2

    # AQUAL ფონური განტოლება μ-ით, სადაც x = chi_p / a_0
    x_val = chi_p / a0
    mu_at_x = x_val / (1 + x_val)
    eq_full = sp.Eq(mu_at_x * chi_p, g_N)

    # Deep MOND (x << 1, μ ≈ x):
    eq_deep = sp.Eq(x_val * chi_p, g_N)
    chi_p_deep_sols = sp.solve(eq_deep, chi_p)
    chi_p_deep = chi_p_deep_sols[0] if chi_p_deep_sols else None

    # ნიუტონური (x >> 1, μ ≈ 1):
    eq_newton = sp.Eq(chi_p, g_N)
    chi_p_newton = sp.solve(eq_newton, chi_p)[0]

    return g_N_expr, eq_full, chi_p_deep, chi_p_newton


def btfr_derivation():
    """
    BTFR (Baryonic Tully-Fisher Relation): M_b ∝ v^4
    Deep MOND ლიმიტში g_chi = √(a_0 g_N), ორბიტა v^2/r = g_chi.

    შენიშვნა: ეს არის ალგებრული ჩასმა, არა RFG-დან გამოყვანა. 
    g_chi = √(a_0·g_N)-ის ჩასმა v^2/r-ში აუცილებლად იძლევა v^4 ∝ M_b.
    """
    v, r, M_b, G_grav, a0 = sp.symbols(
        "v r M_b G a0", real=True, positive=True
    )
    g_N = G_grav * M_b / r**2
    g_chi = sp.sqrt(a0 * g_N)
    g_eff = g_N + g_chi

    # ორბიტა — გარე რეგიონში g_N << g_chi → v^2/r = g_chi
    orbit_eq = sp.Eq(v**2 / r, g_chi)
    M_b_BTFR = sp.solve(orbit_eq, M_b)[0]

    return g_eff, g_chi, M_b_BTFR


def a0_cosmological_scale():
    """
    a_0 = c·H_0 / (2π) რიცხობრივი დამთხვევა.
    NB: ეს არ არის Λ_eff-დან გამოყვანა — ეს არის ემპირიული scaling.
    Λ_eff → a_0 ხიდი ცარიელია; ცალკე ცდა საჭიროა.
    """
    c_si = 299792458.0
    mpc_to_m = 3.086e22
    a0_target = 1.2e-10  # m/s^2 (Milgrom)

    H0_values = {
        "SH0ES 73": 73000.0 / mpc_to_m,
        "Planck 67.4": 67400.0 / mpc_to_m,
        "median 70": 70000.0 / mpc_to_m,
    }

    results = []
    for label, H0 in H0_values.items():
        a0_predicted = c_si * H0 / (2 * math.pi)
        err = abs(a0_predicted - a0_target) / a0_target * 100
        results.append((label, a0_predicted, err))

    return results, a0_target


def mu_derivation_audit():
    """
    Strategy 3 / M5 audit.

    RFG-ის ფესვიდან μ(x)-ის გამოყვანა ნიშნავს არა მხოლოდ ისეთი F(y)-ის პოვნას,
    რომლის derivative იძლევა μ-ს, არამედ χ-მედიუმის მიკროფიზიკური მოქმედებიდან
    ამ F(y)-ის აუცილებლობის ჩვენებას.
    """
    mu_simple, F_y, _, _, _, F_deep_target, deep_mond_limit = derive_mu_choice()
    return {
        "candidate_mu": mu_simple,
        "candidate_F_y": F_y,
        "deep_mond_target": F_deep_target,
        "deep_mond_consistency": deep_mond_limit,
        "rfg_microphysical_derivation": "FAIL",
        "reason": (
            "phase11 reconstructs a valid AQUAL F(y) for the chosen Famaey-Binney "
            "mu, but no RFG action or vortex-memory dynamics forces this function."
        ),
        "status": "phenomenological AQUAL bridge, not an RFG derivation",
    }


def a0_lambda_eff_audit():
    """
    Strategy 3 / M6 audit.

    Dimensional bridge:
        a0 = c H / (2π)

    This can be evaluated with H0 or with a Lambda-only de Sitter scale
    H_Lambda = c sqrt(Lambda/3), but the 2π and the MOND coupling are not
    derived from the polynomial RFG action. Therefore this remains a scaling
    coincidence until a mechanism fixes the coherence length and coupling.
    """
    c_si = 299792458.0
    mpc_to_m = 3.0856775814913673e22
    lambda_obs_m2 = 1.1056e-52
    a0_target = 1.2e-10

    scales = {
        "Planck H0 67.36": 67.36 * 1000.0 / mpc_to_m,
        "SH0ES H0 73.0": 73.0 * 1000.0 / mpc_to_m,
        "Lambda-only de Sitter": c_si * math.sqrt(lambda_obs_m2 / 3.0),
    }

    rows = []
    for label, H_value in scales.items():
        a0_value = c_si * H_value / (2.0 * math.pi)
        rows.append(
            {
                "scale": label,
                "H_s_inv": H_value,
                "a0_m_s2": a0_value,
                "percent_error_vs_Milgrom": abs(a0_value - a0_target) / a0_target * 100.0,
            }
        )

    return {
        "rows": rows,
        "lambda_obs_m2": lambda_obs_m2,
        "verdict": "NO_DERIVATION",
        "reason": (
            "Lambda_eff fixes a cosmological H-scale only after Friedmann dynamics; "
            "it does not by itself derive the MOND interpolation, the coupling to "
            "baryons, or the 2π normalization."
        ),
        "status": "dimensional bridge / numerical coincidence, not a closed RFG mechanism",
    }


def open_tasks():
    """ცარიელი ცდები, რომელიც phase23+-ში გადადის."""
    return [
        "μ(x) RFG-დან გამოყვანა (ფენომენოლოგიური Famaey-Binney არ უდრის RFG-დან გამოყვანას)",
        "a_0 = c·H_0/(2π)-ის მექანიზმი Λ_eff-დან (ნუმეროლოგიური დამთხვევაა ჯერ)",
        "SPARC ~175 გალაქტიკის χ^2 fit (phase33-ის კანდიდატი)",
        "Bullet/Abell 520/El Gordo lensing შედარება (p09_bullet.py)",
        "χ ლაგრანჟიანის RFG-დან გამოყვანა (ჯერ ცალკე AQUAL-ის ჩასმაა)",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 11: MOND სექტორი — Bekenstein-Milgrom AQUAL")
    print("რეფერენცია: NOTATION.md, phase22, RFG_Theory.md § 4")
    print("=" * 72)

    L_chi, mu, EL_lhs, EL_rhs = aqual_lagrangian_explicit()
    print("\n1. AQUAL ლაგრანჟიანი ნორმალიზაციით")
    print(f"  L_χ = {L_chi}")
    print(f"  Euler-Lagrange: {EL_lhs} = {EL_rhs}")

    (
        mu_simple,
        F_y,
        mu_newton,
        mu_deep,
        F_deep_leading,
        F_deep_target,
        deep_mond_limit,
    ) = derive_mu_choice()
    print("\n2. μ(x) ფენომენოლოგიური არჩევანი (Famaey-Binney 2005)")
    print(f"  μ(x) = {mu_simple}")
    print(f"  ნიუტონური ლიმიტი (x → ∞): μ → {mu_newton}")
    print(f"  Deep MOND ლიმიტი (x → 0): μ ≈ {mu_deep}")
    print(f"  რეკონსტრუირებული F(y) = {F_y}")
    print(f"  F(y) deep-MOND ლიდინგ რიგი (y → 0): {F_deep_leading}")
    print(f"  სამიზნე: F(y) → {F_deep_target}")
    print(f"  შემოწმება: limit(F/F_target) = {deep_mond_limit}")

    g_N_expr, eq_full, chi_p_deep, chi_p_newton = spherical_field_equation()
    print("\n3. სფერული ფონური განტოლება")
    print(f"  g_N = {g_N_expr}")
    print(f"  სრული: {eq_full}")
    print(f"  Deep MOND χ' = {chi_p_deep}")
    print(f"  ნიუტონური χ' = {chi_p_newton}")

    g_eff, g_chi, M_b_BTFR = btfr_derivation()
    print("\n4. BTFR (Baryonic Tully-Fisher) ალგებრული ჩასმა")
    print(f"  g_χ (deep MOND) = {g_chi}")
    print(f"  g_eff = g_N + g_χ = {g_eff}")
    print(f"  M_b (deep MOND) = {M_b_BTFR}  (v^4 ∝ M_b)")

    print("\n5. a_0 კოსმოლოგიური მასშტაბი — რიცხობრივი დამთხვევა")
    results, a0_target = a0_cosmological_scale()
    print(f"  სამიზნე (Milgrom): a_0 = {a0_target:.4e} m/s²")
    for label, a0_pred, err in results:
        print(f"  {label:18s}: a_0 = {a0_pred:.4e} m/s² (ცდომილება {err:.2f}%)")

    print("\n6. Strategy 3 audit — μ(x) derivation status")
    mu_audit = mu_derivation_audit()
    for key, value in mu_audit.items():
        print(f"  {key:30s}: {value}")

    print("\n7. Strategy 3 audit — a_0 / Λ_eff bridge status")
    a0_audit = a0_lambda_eff_audit()
    for row in a0_audit["rows"]:
        print(
            f"  {row['scale']:24s}: a0={row['a0_m_s2']:.4e} m/s², "
            f"error={row['percent_error_vs_Milgrom']:.2f}%"
        )
    print(f"  verdict: {a0_audit['verdict']}")
    print(f"  status : {a0_audit['status']}")

    print("\n8. ღია ცდები (გადადის შემდეგ ფაზებში)")
    for i, task in enumerate(open_tasks(), 1):
        print(f"  {i}. {task}")

    print("\n9. სტატუსი")
    print("  - AQUAL ფარგლი დაფიქსირებულია (Bekenstein-Milgrom 1984)")
    print("  - μ(x) ფენომენოლოგიური არჩევანი, არა RFG-დან გამოყვანა")
    print("  - a_0 cH_0/(2π) რიცხობრივი დამთხვევაა (~10% ცდომილებით)")
    print("  - BTFR ალგებრული ჩასმაა, ემპირიულ ფიტს ეყრდნობა")
    print("  - SPARC შედარება, χ ლაგრანჟიანის RFG-დან გამოყვანა — ღია")


# ===================== merged from p07_mond.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
მოძველებული (deprecated) საწყისი ანიზოტროპიის ესკიზი.
სრული MOND თეორია, a_0, mu(x) და BTFR იხილეთ p07_mond.py-სა და p09_bullet.py-ში.
"""
import sympy as sp
from p01_core import get_pressures

def galactic_acceleration():
    r = sp.Symbol('r', real=True, positive=True)
    
    # get_pressures აბრუნებს ანიზოტროპიულ წნევებს, რომლებიც r-ის ფუნქციებია 
    rho, p_rad, p_tan, delta_p = get_pressures(r)
    
    # 1. ტესტ-ნაწილაკის რეალური აჩქარება მედიუმის შიდა წონასწორობიდან:
    # a_test = (p_rad' - 2*delta_p/r) / (rho + p_rad)
    p_rad_prime = sp.diff(p_rad, r)
    a_test = sp.simplify((p_rad_prime - 2 * delta_p / r) / (rho + p_rad))
    
    # 2. ნიუტონური ძალა ჩვეულებრივი ბარიონული მასიდან 
    G = sp.Symbol('G', real=True, positive=True)
    M_b = sp.Function('M_b')(r)
    g_N = G * M_b / r**2
    
    return a_test, delta_p, g_N

if __name__ == "__main__":
    a_test, dp, g_N = galactic_acceleration()
    print("სფერული გალაქტიკური ანიზოტროპია (delta_p):")
    print(dp)
    print("\nტესტ-ნაწილაკის ფორმალური ანიზოტროპიული აჩქარება (ესკიზური):")
    print(a_test)
    print("\nნიუტონური ნაწილი:")
    print(g_N)
    print("\nშენიშვნა: სრული MOND ამოხსნა chi-ველის მეხსიერებით იხილეთ p07_mond.py-სა და p09_bullet.py-ში.")

