# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 31: ნეიტრონული ვარსკვლავები — RFG ანიზოტროპიული TOV და M-R მრუდი
================================================================================

რეფერენცია: phase10_matter_coupling.py, phase1_tensor.py,
            phase22_full_stress_tensor.py, STRATEGY.md S3/E6

სტატუსი:
ეს არის TOV განტოლების შესამოწმებელი (executable) გარემო. ის არ წარმოადგენს სრულ 
nuclear-EOS მორგებას და არ არის საბოლოო NICER/GW170817 likelihood ცდა. 
RFG-ის ანიზოტროპია წარმოდგენილია ერთი ფენომენოლოგიური პარამეტრით (eta_delta):

    Delta p = p_tan - p_rad = eta_delta * p_rad * u,
    u = 2GM(r)/(r c^2).

eta_delta = 0 არის GR/იზოტროპული TOV ლიმიტი. დადებითი eta_delta ქმნის დამატებით 
ტანგენციალურ მხარდაჭერას და ხდის RFG-ის მსგავსი ანიზოტროპიული სტრესის გავლენას 
ხილულს. RFG-დან უშუალოდ გამოყვანილი ნეიტრონული ვარსკვლავის EoS ჯერჯერობით ღიაა.

რა არის იმპლემენტირებული:
- SI ერთეულების მქონე GR TOV ინტეგრატორი RK4 მეთოდით.
- RFG-ის ანიზოტროპიული TOV წევრი +2 Delta p/r.
- M-R (მასა-რადიუსი) მიმდევრობა ცენტრალური სიმკვრივის ცვლილებით.
- M_max >= 2.08 M_sun შემოწმება.
- Lambda_1.4 კომპაქტურობის პროქსი GW170817 ზღვრისთვის.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


M_SUN = 1.98847e30
G = 6.67430e-11
C = 299792458.0


OBSERVATIONS = {
    "NICER_J0030": {
        "mass_solar": 1.34,
        "mass_err": 0.16,
        "radius_km": 12.71,
        "radius_err_km": 1.14,
        "reference": "Miller et al. 2019",
    },
    "NICER_J0740": {
        "mass_solar": 2.08,
        "mass_err": 0.07,
        "radius_km": 13.7,
        "radius_err_km": 1.5,
        "reference": "Miller et al. 2021",
    },
    "GW170817_tidal_Lambda_1.4": {
        "upper_bound_90CL": 580,
        "reference": "Abbott et al. PRL 122:061104 reanalysis; earlier bound <800",
    },
    "PSR_J0348+0432": {
        "mass_solar": 2.01,
        "mass_err": 0.04,
        "reference": "Antoniadis et al. 2013",
    },
}


@dataclass(frozen=True)
class PolytropicEOS:
    """
    Rest-mass polytrope:
        p = K * rho0^gamma
        epsilon/c^2 = rho0 + p/((gamma - 1)c^2)

    K is chosen only as a controlled toy EOS. It is tuned to neutron-star
    scales but not claimed as SLy/APR/nuclear-matter inference.
    """

    name: str = "toy_polytrope_Gamma2.4"
    K: float = 1.20e-9
    gamma: float = 2.40

    def pressure_from_rest_density(self, rho0_kg_m3: float) -> float:
        return self.K * rho0_kg_m3**self.gamma

    def rest_density_from_pressure(self, pressure_pa: float) -> float:
        if pressure_pa <= 0:
            return 0.0
        return (pressure_pa / self.K) ** (1.0 / self.gamma)

    def mass_density_from_pressure(self, pressure_pa: float) -> float:
        rho0 = self.rest_density_from_pressure(pressure_pa)
        if rho0 <= 0:
            return 0.0
        return rho0 + pressure_pa / ((self.gamma - 1.0) * C**2)

    def sound_speed_sq(self, pressure_pa: float) -> float:
        rho0 = self.rest_density_from_pressure(pressure_pa)
        if rho0 <= 0:
            return 0.0
        dp_drho0 = self.gamma * self.K * rho0 ** (self.gamma - 1.0)
        drho_mass_drho0 = 1.0 + dp_drho0 / ((self.gamma - 1.0) * C**2)
        return dp_drho0 / drho_mass_drho0


@dataclass
class StarSolution:
    central_density_kg_m3: float
    eta_delta: float
    mass_solar: float
    radius_km: float
    compactness: float
    lambda_proxy: float
    max_sound_speed_over_c: float
    status: str


@dataclass
class SequenceSummary:
    label: str
    eta_delta: float
    n_models: int
    max_mass_solar: float
    radius_at_max_km: float
    rho_c_at_max: float
    radius_1p4_km: float | None
    lambda_1p4_proxy: float | None
    supports_2p08: bool
    lambda_bound_pass: bool | None
    max_sound_speed_over_c: float


def anisotropy_delta_p(pressure_pa: float, compactness_u: float, eta_delta: float) -> float:
    """Phenomenological RFG anisotropy: Delta p = eta * p * u."""
    return eta_delta * pressure_pa * compactness_u


def tov_derivatives(
    eos: PolytropicEOS,
    eta_delta: float,
    radius_m: float,
    mass_kg: float,
    pressure_pa: float,
) -> tuple[float, float] | None:
    """Return dm/dr and dp/dr for anisotropic TOV in SI units."""
    if pressure_pa <= 0:
        return 0.0, 0.0
    if radius_m <= 0:
        return None

    rho = eos.mass_density_from_pressure(pressure_pa)
    compact_factor = 1.0 - 2.0 * G * mass_kg / (radius_m * C**2)
    if compact_factor <= 1.0e-6:
        return None

    local_u = 2.0 * G * mass_kg / (radius_m * C**2)
    delta_p = anisotropy_delta_p(pressure_pa, local_u, eta_delta)

    dm_dr = 4.0 * math.pi * radius_m**2 * rho
    dp_gr = (
        -G
        * (rho + pressure_pa / C**2)
        * (mass_kg + 4.0 * math.pi * radius_m**3 * pressure_pa / C**2)
        / (radius_m**2 * compact_factor)
    )
    dp_dr = dp_gr + 2.0 * delta_p / radius_m

    return dm_dr, dp_dr


def rk4_step(
    eos: PolytropicEOS,
    eta_delta: float,
    radius_m: float,
    mass_kg: float,
    pressure_pa: float,
    dr_m: float,
) -> tuple[float, float] | None:
    """One RK4 step for (m, p)."""
    k1 = tov_derivatives(eos, eta_delta, radius_m, mass_kg, pressure_pa)
    if k1 is None:
        return None

    k2 = tov_derivatives(
        eos,
        eta_delta,
        radius_m + 0.5 * dr_m,
        mass_kg + 0.5 * dr_m * k1[0],
        pressure_pa + 0.5 * dr_m * k1[1],
    )
    if k2 is None:
        return None

    k3 = tov_derivatives(
        eos,
        eta_delta,
        radius_m + 0.5 * dr_m,
        mass_kg + 0.5 * dr_m * k2[0],
        pressure_pa + 0.5 * dr_m * k2[1],
    )
    if k3 is None:
        return None

    k4 = tov_derivatives(
        eos,
        eta_delta,
        radius_m + dr_m,
        mass_kg + dr_m * k3[0],
        pressure_pa + dr_m * k3[1],
    )
    if k4 is None:
        return None

    next_mass = mass_kg + dr_m * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0
    next_pressure = pressure_pa + dr_m * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0
    return next_mass, next_pressure


def dimensionless_lambda_proxy(mass_solar: float, radius_km: float, k2: float = 0.08) -> float:
    """
    Lambda proxy: Lambda = (2/3) k2 / C_compact^5.

    This is not a Love-number integration. It is a compactness diagnostic for
    the GW170817 Lambda_1.4 scale.
    """
    compactness = G * mass_solar * M_SUN / ((radius_km * 1000.0) * C**2)
    if compactness <= 0:
        return math.inf
    return (2.0 / 3.0) * k2 / compactness**5


def integrate_star(
    eos: PolytropicEOS,
    central_density_kg_m3: float,
    eta_delta: float = 0.0,
    dr_m: float = 50.0,
    max_radius_km: float = 50.0,
) -> StarSolution:
    """Integrate one stellar model until pressure reaches zero."""
    pressure = eos.pressure_from_rest_density(central_density_kg_m3)
    radius = dr_m
    rho_c = eos.mass_density_from_pressure(pressure)
    mass = 4.0 * math.pi * radius**3 * rho_c / 3.0
    previous = (radius, mass, pressure)
    max_sound_speed = math.sqrt(max(eos.sound_speed_sq(pressure), 0.0)) / C

    while radius < max_radius_km * 1000.0 and pressure > 0:
        step = rk4_step(eos, eta_delta, radius, mass, pressure, dr_m)
        if step is None:
            return StarSolution(
                central_density_kg_m3=central_density_kg_m3,
                eta_delta=eta_delta,
                mass_solar=mass / M_SUN,
                radius_km=radius / 1000.0,
                compactness=2.0 * G * mass / (radius * C**2),
                lambda_proxy=math.inf,
                max_sound_speed_over_c=max_sound_speed,
                status="STOP_COMPACTNESS",
            )

        previous = (radius, mass, pressure)
        mass, pressure = step
        radius += dr_m
        if pressure > 0:
            cs_over_c = math.sqrt(max(eos.sound_speed_sq(pressure), 0.0)) / C
            max_sound_speed = max(max_sound_speed, cs_over_c)

    if pressure <= 0:
        r0, m0, p0 = previous
        fraction = p0 / (p0 - pressure) if p0 != pressure else 0.0
        surface_radius = r0 + fraction * (radius - r0)
        surface_mass = m0 + fraction * (mass - m0)
        mass_solar = surface_mass / M_SUN
        radius_km = surface_radius / 1000.0
        compactness = 2.0 * G * surface_mass / (surface_radius * C**2)
        return StarSolution(
            central_density_kg_m3=central_density_kg_m3,
            eta_delta=eta_delta,
            mass_solar=mass_solar,
            radius_km=radius_km,
            compactness=compactness,
            lambda_proxy=dimensionless_lambda_proxy(mass_solar, radius_km),
            max_sound_speed_over_c=max_sound_speed,
            status="OK",
        )

    return StarSolution(
        central_density_kg_m3=central_density_kg_m3,
        eta_delta=eta_delta,
        mass_solar=mass / M_SUN,
        radius_km=radius / 1000.0,
        compactness=2.0 * G * mass / (radius * C**2),
        lambda_proxy=math.inf,
        max_sound_speed_over_c=max_sound_speed,
        status="STOP_MAX_RADIUS",
    )


def central_density_grid() -> list[float]:
    """Central rest-density sweep."""
    return [3.5e17 * (1.08**index) for index in range(30)]


def build_mr_sequence(eos: PolytropicEOS, eta_delta: float) -> list[StarSolution]:
    sequence = []
    for rho_c in central_density_grid():
        solution = integrate_star(eos, rho_c, eta_delta=eta_delta)
        if solution.status == "OK" and solution.mass_solar > 0:
            sequence.append(solution)
    return sequence


def stable_branch(sequence: list[StarSolution]) -> list[StarSolution]:
    if not sequence:
        return []
    max_index = max(range(len(sequence)), key=lambda idx: sequence[idx].mass_solar)
    return sequence[: max_index + 1]


def interpolate_at_mass(
    sequence: list[StarSolution],
    target_mass_solar: float,
) -> tuple[float, float] | None:
    branch = stable_branch(sequence)
    if len(branch) < 2:
        return None

    for left, right in zip(branch, branch[1:]):
        m0, m1 = left.mass_solar, right.mass_solar
        if (m0 <= target_mass_solar <= m1) or (m1 <= target_mass_solar <= m0):
            if abs(m1 - m0) < 1.0e-12:
                radius = 0.5 * (left.radius_km + right.radius_km)
            else:
                frac = (target_mass_solar - m0) / (m1 - m0)
                radius = left.radius_km + frac * (right.radius_km - left.radius_km)
            return radius, dimensionless_lambda_proxy(target_mass_solar, radius)
    return None


def summarize_sequence(
    label: str,
    eta_delta: float,
    sequence: list[StarSolution],
) -> SequenceSummary:
    if not sequence:
        raise ValueError("empty M-R sequence")

    max_solution = max(sequence, key=lambda item: item.mass_solar)
    onep4 = interpolate_at_mass(sequence, 1.4)
    lambda_bound = OBSERVATIONS["GW170817_tidal_Lambda_1.4"]["upper_bound_90CL"]

    return SequenceSummary(
        label=label,
        eta_delta=eta_delta,
        n_models=len(sequence),
        max_mass_solar=max_solution.mass_solar,
        radius_at_max_km=max_solution.radius_km,
        rho_c_at_max=max_solution.central_density_kg_m3,
        radius_1p4_km=onep4[0] if onep4 else None,
        lambda_1p4_proxy=onep4[1] if onep4 else None,
        supports_2p08=max_solution.mass_solar >= OBSERVATIONS["NICER_J0740"]["mass_solar"],
        lambda_bound_pass=(onep4[1] < lambda_bound if onep4 else None),
        max_sound_speed_over_c=max(item.max_sound_speed_over_c for item in sequence),
    )


def nearest_mass_model(sequence: list[StarSolution], target_mass_solar: float) -> StarSolution:
    return min(sequence, key=lambda item: abs(item.mass_solar - target_mass_solar))


def format_summary(summary: SequenceSummary) -> list[str]:
    return [
        f"label: {summary.label}",
        f"eta_delta: {summary.eta_delta:.3g}",
        f"models: {summary.n_models}",
        f"M_max: {summary.max_mass_solar:.3f} M_sun at R={summary.radius_at_max_km:.2f} km",
        f"rho_c(M_max): {summary.rho_c_at_max:.3e} kg/m^3",
        f"R_1.4: {summary.radius_1p4_km:.2f} km" if summary.radius_1p4_km else "R_1.4: not bracketed",
        (
            f"Lambda_1.4 proxy: {summary.lambda_1p4_proxy:.0f}"
            if summary.lambda_1p4_proxy
            else "Lambda_1.4 proxy: not bracketed"
        ),
        f"M_max >= 2.08 M_sun: {summary.supports_2p08}",
        f"Lambda_1.4 < 580 proxy: {summary.lambda_bound_pass}",
        f"max c_s/c in sequence: {summary.max_sound_speed_over_c:.3f}",
    ]


def model_scope_notes() -> list[str]:
    return [
        "The polytropic EOS is a controlled toy EOS, not a nuclear-matter fit.",
        "eta_delta is a phenomenological stand-in for phase1_tensor Delta p.",
        "Lambda_1.4 is a compactness proxy; real tidal deformability needs Love-number ODEs.",
        "A real NICER/GW170817 verdict needs EOS priors and Bayesian likelihoods.",
    ]


def main() -> None:
    print("=" * 72)
    print("PHASE 31: Neutron stars — RFG anisotropic TOV and M-R curve")
    print("=" * 72)

    eos = PolytropicEOS()
    models = {
        "GR_isotropic": 0.0,
        "RFG_anisotropic_eta0.5": 0.5,
    }

    print("\n1. Observational filters")
    for key, obs in OBSERVATIONS.items():
        print(f"  {key}: {obs}")

    print("\n2. EOS and anisotropy model")
    print(f"  EOS: {eos.name}, K={eos.K:.3e}, gamma={eos.gamma:.2f}")
    print("  anisotropy: Delta p = eta_delta * p_rad * 2GM/(r c^2)")

    summaries: list[SequenceSummary] = []
    sequences: dict[str, list[StarSolution]] = {}
    for label, eta_delta in models.items():
        sequence = build_mr_sequence(eos, eta_delta)
        sequences[label] = sequence
        summary = summarize_sequence(label, eta_delta, sequence)
        summaries.append(summary)

    print("\n3. M-R sequence summaries")
    for summary in summaries:
        print(f"\n  --- {summary.label} ---")
        for line in format_summary(summary):
            print(f"  {line}")

    print("\n4. Sequence samples")
    for label, sequence in sequences.items():
        print(f"\n  {label}")
        stride = max(1, len(sequence) // 5)
        for solution in sequence[::stride][:6]:
            print(
                f"    rho_c={solution.central_density_kg_m3:.2e} kg/m^3 | "
                f"M={solution.mass_solar:.3f} M_sun | R={solution.radius_km:.2f} km | "
                f"u={solution.compactness:.3f}"
            )

    print("\n5. NICER/GW170817 quick comparisons")
    for label, sequence in sequences.items():
        j0030 = nearest_mass_model(sequence, OBSERVATIONS["NICER_J0030"]["mass_solar"])
        j0740 = nearest_mass_model(sequence, OBSERVATIONS["NICER_J0740"]["mass_solar"])
        print(f"\n  {label}")
        print(
            f"    nearest J0030 mass: M={j0030.mass_solar:.3f}, R={j0030.radius_km:.2f} km "
            f"(obs R={OBSERVATIONS['NICER_J0030']['radius_km']}±{OBSERVATIONS['NICER_J0030']['radius_err_km']} km)"
        )
        print(
            f"    nearest J0740 mass: M={j0740.mass_solar:.3f}, R={j0740.radius_km:.2f} km "
            f"(mass target {OBSERVATIONS['NICER_J0740']['mass_solar']} M_sun)"
        )

    print("\n6. Scope notes")
    for note in model_scope_notes():
        print(f"  - {note}")

    print("\n7. Status")
    print("  - Strategy E6 TOV integrator: implemented.")
    print("  - M-R curve: generated for GR and RFG-anisotropic toy branches.")
    print("  - Delta p effect: quantified by eta_delta branch comparison.")
    print("  - Full EOS/Love-number/Bayesian fit: still open.")


if __name__ == "__main__":
    main()
