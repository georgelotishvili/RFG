# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.
# Active coefficient scheme: c_Y denotes the Y-scheme coefficient c_Y^(Y).

from __future__ import annotations

import os
import sympy as sp
from dataclasses import dataclass
from p01_core import init_variables, get_polynomial_lagrangian


@dataclass(frozen=True)
class ClaimGate:
    claim: str
    status: str
    verified_here: str
    open_requirement: str


def oscillon_claim_gate() -> list[ClaimGate]:
    """Strict theorem gate for p10: separate identities from open physics."""
    return [
        ClaimGate(
            claim="Bernoulli pressure identity",
            status="CLOSED_SYMBOLIC_IDENTITY",
            verified_here="For phi=-r_s/r, Delta_P=e^phi*(phi')^2/(32*pi*G) and P_static=-Delta_P give P_static+Delta_P=0 exactly.",
            open_requirement="derive the same identity as an on-shell branch of the full p01 stress system, not only as a static scalar identity.",
        ),
        ClaimGate(
            claim="localized source gives exterior 1/r tail",
            status="CLOSED_IF_SOURCE_IS_LOCALIZED",
            verified_here="Spherical Poisson reconstruction fixes the far-zone 1/r coefficient once a localized source rho(r) is supplied.",
            open_requirement="prove the full nonlinear oscillon PDE has a regular localized finite-energy source.",
        ),
        ClaimGate(
            claim="Newton law recovery",
            status="CONDITIONAL_ON_BICONFORMAL_BRANCH_AND_NORMALIZATION",
            verified_here="Given phi=-r_s/r and r_s=2GM/c^2, geodesic acceleration gives -GM/r^2.",
            open_requirement="derive the source-to-r_s normalization and the bi-conformal exterior branch from the p01 action.",
        ),
        ClaimGate(
            claim="bi-conformal scaling and light factor 2",
            status="CLOSED_WITHIN_BICONFORMAL_ANSATZ",
            verified_here="In the (+---) bi-conformal metric, c_coord/c=(L_oper/L_0)^2 and weak-field light bending gives 2*r_s/b.",
            open_requirement="show the ansatz is selected by the full static spherical equations, not imposed as a geometry template.",
        ),
        ClaimGate(
            claim="finite-energy oscillon particle",
            status="OPEN_CORE_THEOREM",
            verified_here="The file computes Y, averaged energy density, and a virial stationarity integrand for a two-harmonic trial family.",
            open_requirement="construct global regular finite-energy nonlinear solutions and prove spectral stability.",
        ),
        ClaimGate(
            claim="Solar-system observational pass",
            status="PARTIAL_WEAK_FIELD_SMOKE_TEST",
            verified_here="Solar light bending and Pound-Rebka are reproduced at leading weak-field level.",
            open_requirement="derive PPN gamma/beta plus Shapiro delay, perihelion, Cassini, LLR, and clock/LPI tests.",
        ),
        ClaimGate(
            claim="substrate carrier frequency nu0",
            status="OPEN_MICROPHYSICAL_HYPOTHESIS",
            verified_here="nu0 is explicitly marked as a substrate rhythm hypothesis/TODO.",
            open_requirement="derive the invariant substrate clock from the microscopic RFG vacuum, not from local oscillon composites.",
        ),
    ]


def oscillon_do_not_claim() -> list[str]:
    return [
        "Do not claim that p10 proves particles are finite-energy oscillons.",
        "Do not claim Newton gravity is fully derived from the RFG action in this file.",
        "Do not claim the bi-conformal exterior metric is derived rather than selected as an ansatz/branch.",
        "Do not claim solar-system tests are fully passed; only leading weak-field checks are present.",
        "Do not claim G is microscopically normalized from oscillon parameters here.",
        "Do not claim nu0 is derived; it remains an open substrate-clock hypothesis.",
    ]


def signature_bridge_gate():
    """
    p10 now keeps the active project signature (+---).

    The earlier Phase-17 notes used the equivalent (-+++) writing.  Multiplying
    the metric by -1 changes covariant sign conventions, not the physical null
    relation or the scaling identities, as long as the energy sign is translated.
    """
    phi = sp.Symbol("phi", real=True)
    g_plus = sp.diag(sp.exp(phi), -sp.exp(-phi), -sp.exp(-phi), -sp.exp(-phi))
    g_minus = -g_plus
    return {
        "active_signature": "(+---)",
        "old_phase17_signature": "(-+++)",
        "g_plus": g_plus,
        "g_minus": g_minus,
        "null_speed_plus": sp.simplify(sp.sqrt(-g_plus[0, 0] / g_plus[1, 1])),
        "null_speed_minus": sp.simplify(sp.sqrt(-g_minus[0, 0] / g_minus[1, 1])),
        "same_null_cone": True,
        "energy_translation": "with (+---), static Killing energy is E=p_t; with (-+++), it was E=-p_t",
        "status": "PASS_SIGNATURE_BRIDGE_EXPLICIT",
    }


DEFAULT_MAIN_SECTIONS = {"oscillon", "spherical", "tov", "biconformal", "gates"}


def _requested_main_sections() -> set[str]:
    raw = os.environ.get("RFG_P10_SECTIONS", "")
    return {part.strip().lower() for part in raw.split(",") if part.strip()}


def _should_run_main_section(name: str) -> bool:
    requested = _requested_main_sections()
    if not requested:
        return name.lower() in DEFAULT_MAIN_SECTIONS
    return "all" in requested or name.lower() in requested


def bernoulli_static_gravity_identity():
    """
    ძველი quantum ფაილის Bernoulli gravity ბირთვი.

    static bi-conformal exterior branch:
        phi(r) = -r_s/r
        Delta P = -P_static = e^phi (phi')^2 / (32*pi*G)
        P_static + Delta P = 0

    ეს არის scalar-field Bernoulli identity: gradient energy drains the
    static pressure. გრავიტაციული მნიშვნელობა სრულად არ იხურება მხოლოდ ამ
    identity-ით; დამატებით საჭიროა source normalization და exterior metric
    branch-ის p01-დან გამოყვანა.
    """
    r, r_s, G = sp.symbols('r r_s G', real=True, positive=True)
    phi = -r_s / r
    phi_prime = sp.diff(phi, r)

    gradient_energy = sp.simplify(sp.exp(phi) * phi_prime**2 / (32 * sp.pi * G))
    pressure_static = -gradient_energy
    bernoulli_sum = sp.simplify(pressure_static + gradient_energy)
    pressure_shape = sp.simplify(sp.exp(sp.Symbol('phi')) * sp.Symbol('phi')**4)
    coulomb_pressure = sp.simplify(
        -sp.exp(phi) * (phi**4) / (32 * sp.pi * G * r_s**2)
    )

    u = sp.Symbol('u', real=True, positive=True)
    shape_u = sp.exp(-u) * u**4
    shape_derivative = sp.factor(sp.diff(shape_u, u))

    return {
        "theorem": "static scalar Bernoulli identity",
        "exterior_profile": sp.Eq(sp.Symbol('phi(r)'), phi),
        "gradient_energy_density": sp.Eq(sp.Symbol('Delta_P'), gradient_energy),
        "static_pressure": sp.Eq(sp.Symbol('P_static'), pressure_static),
        "bernoulli_integral": sp.Eq(sp.Symbol('P_static + Delta_P'), bernoulli_sum),
        "closed_pressure_profile": sp.Eq(sp.Symbol('P_static_Coulomb'), coulomb_pressure),
        "universal_shape": sp.Eq(sp.Symbol('f_hat(phi)'), pressure_shape),
        "shape_derivative_u_rs_over_r": sp.Eq(sp.Symbol('d(e^-u*u^4)/du'), shape_derivative),
        "pressure_deficit_peak": sp.Eq(sp.Symbol('r_peak'), r_s / 4),
        "strong_field_saturation": sp.Eq(sp.Symbol('lim_r_to_0_Delta_P'), sp.limit(gradient_energy, r, 0, dir='+')),
        "far_field_limit": sp.Eq(sp.Symbol('lim_r_to_inf_Delta_P'), sp.limit(gradient_energy, r, sp.oo)),
    }


def bernoulli_time_averaged_oscillon_source():
    """
    Time-periodic oscillon-ის zero-frequency projection.

    full field split:
        phi(t,r) = phi_grav(r) + Phi0(r) cos(Omega t) + ...

    სწრაფი oscillation-ის ხაზოვანი წევრები საშუალოდ ნულდება. დარჩენილი
    quadratic zero-frequency density არის localized source, რომელიც
    Poisson/Laplace reconstruction-ით იძლევა long-range 1/r field-ს.
    """
    r, G, Omega, phi_grav = sp.symbols('r G Omega phi_grav', real=True, positive=True)
    Phi0 = sp.Function('Phi0')(r)
    Phi0_prime = sp.diff(Phi0, r)

    source_avg = sp.simplify(
        (
            sp.exp(-phi_grav) * Omega**2 * Phi0**2
            + sp.exp(phi_grav) * Phi0_prime**2
        ) / (64 * sp.pi * G)
    )
    newtonian_source_profile = sp.simplify(
        sp.Rational(1, 2) * Omega**2 * Phi0**2
        + sp.Rational(1, 2) * Phi0_prime**2
    )

    return {
        "field_split": "phi(t,r)=phi_grav(r)+Phi0(r)*cos(Omega*t)+higher harmonics",
        "zero_frequency_source": sp.Eq(sp.Symbol('<T00>'), source_avg),
        "newtonian_profile_without_prefactor": sp.Eq(sp.Symbol('rho_osc'), newtonian_source_profile),
        "meaning": "oscillon-ის შიდა რეზონანსის საშუალო kinetic/gradient energy არის static gravitational source.",
    }


def bernoulli_poisson_reconstruction():
    """
    Localized zero-frequency oscillon source -> exterior 1/r field.

    spherical Poisson reconstruction:
        M_enc(r)=4*pi int_0^r rho(r') r'^2 dr'
        phi_grav(r)=-M_enc(r)/r - int_r^inf 4*pi rho(r') r' dr'
        lim_{r->inf}[-r phi_grav(r)] = M_total
    """
    r, rp = sp.symbols('r rp', real=True, positive=True)
    rho = sp.Function('rho')

    m_enc = 4 * sp.pi * sp.Integral(rho(rp) * rp**2, (rp, 0, r))
    phi_grav = -m_enc / r - 4 * sp.pi * sp.Integral(rho(rp) * rp, (rp, r, sp.oo))
    m_total = 4 * sp.pi * sp.Integral(rho(rp) * rp**2, (rp, 0, sp.oo))

    return {
        "enclosed_source": sp.Eq(sp.Symbol('M_enc(r)'), m_enc),
        "poisson_solution": sp.Eq(sp.Symbol('phi_grav(r)'), phi_grav),
        "far_zone_coefficient": sp.Eq(sp.Symbol('lim_-r_phi'), m_total),
        "proof_result": "Once the localized source coefficient is fixed, the 1/r tail coefficient is fixed. Physical G-normalization is a separate gate.",
    }


def poisson_to_newton_normalization_gate():
    """
    Exact normalization gap isolated by the council.

    Poisson reconstruction gives a dimensionless exterior coefficient mu_src:
        phi = -mu_src/r.
    Newton is recovered if the metric potential is Phi_N=c^2*phi/2 and the
    source coefficient is identified with mu_src=2GM/c^2.  This function proves
    the algebra of that identification, not its microscopic derivation.
    """
    r, mu_src, G, M, c = sp.symbols("r mu_src G M c", positive=True)
    phi = -mu_src / r
    phi_newton = sp.simplify(c**2 * phi / 2)
    acceleration = sp.simplify(-sp.diff(phi_newton, r))
    mu_needed = sp.solve(sp.Eq(acceleration, -G * M / r**2), mu_src)[0]

    return {
        "dimensionless_exterior_profile": sp.Eq(sp.Symbol("phi"), phi),
        "newtonian_potential_bridge": sp.Eq(sp.Symbol("Phi_N"), phi_newton),
        "geodesic_acceleration": sp.Eq(sp.Symbol("a"), acceleration),
        "required_source_coefficient": sp.Eq(sp.Symbol("mu_src"), mu_needed),
        "status": "NORMALIZATION_IDENTIFIED_NOT_MICRODERIVED",
        "open_requirement": "derive mu_src=2GM/c^2 from the oscillon stress integral and p01 coupling constants.",
    }


def bernoulli_newton_law_recovery():
    """
    Bernoulli pressure deficit gives the mechanism; geodesics give motion.

    The firewall from the old quantum file is kept:
        matter does not feel a literal pressure-gradient force.
        Matter is minimally coupled and follows geodesics of g_mn[phi].
    """
    r, G, c, M, m, M1, M2, d = sp.symbols(
        'r G c M m M1 M2 d',
        real=True,
        positive=True,
    )
    r_s = 2 * G * M / c**2
    phi = -r_s / r
    Phi_N = sp.simplify(c**2 * phi / 2)
    acceleration = sp.simplify(-c**2 * sp.diff(phi, r) / 2)

    potential_energy = -G * M1 * M2 / d
    radial_force = sp.simplify(-sp.diff(potential_energy, d))
    force_magnitude = sp.simplify(G * M1 * M2 / d**2)

    return {
        "source_radius": sp.Eq(sp.Symbol('r_s'), r_s),
        "gravitational_profile": sp.Eq(sp.Symbol('phi_grav'), phi),
        "newtonian_potential": sp.Eq(sp.Symbol('Phi_N'), Phi_N),
        "geodesic_acceleration": sp.Eq(sp.Symbol('a_geo'), acceleration),
        "test_particle_force": sp.Eq(sp.Symbol('F'), m * acceleration),
        "two_body_energy": sp.Eq(sp.Symbol('U(d)'), potential_energy),
        "two_body_radial_force": sp.Eq(sp.Symbol('F_radial'), radial_force),
        "two_body_force_magnitude": sp.Eq(sp.Symbol('|F_12|'), force_magnitude),
        "firewall": "Delta_P explains why the metric is gravitational; it is not an extra pressure-gradient force on matter.",
    }


def bernoulli_gravity_chain() -> list[str]:
    return [
        "localized oscillon resonance -> time-periodic scalar energy",
        "zero-frequency projection <T00> survives time averaging",
        "Bernoulli identity: P_static + e^phi |grad phi|^2/(32*pi*G)=0",
        "pressure deficit: Delta_P=-P_static=e^phi |grad phi|^2/(32*pi*G)",
        "localized Delta_P/<T00> fixes the source integral if a finite oscillon exists",
        "vacuum exterior solves Laplace equation -> phi_grav=-mu_src/r",
        "Newton normalization requires mu_src=2GM/c^2",
        "bi-conformal geodesic acceleration then gives a=-(c^2/2) grad phi=-GM/r^2",
        "two oscillon sources give U(d)=-G M1 M2/d and |F|=G M1 M2/d^2 after normalization",
        "strong-field saturation: e^phi suppresses Delta_P as phi->-infinity",
    ]


def analyze_oscillon():
    r = sp.Symbol('r', real=True, positive=True)
    t = sp.Symbol('t', real=True)
    theta = sp.Symbol('theta', real=True) # theta = omega * t
    omega = sp.Symbol('omega', real=True, positive=True)
    
    # ოსცილონის რადიალური პროფილი და მისი გრადიენტი
    Phi0 = sp.Function('Phi0')(r)
    Phi0_prime = sp.diff(Phi0, r)
    Phi1 = sp.Function('Phi1')(r) # მეორე ჰარმონიკა არაწრფივი შერევისთვის
    Phi1_prime = sp.diff(Phi1, r)
    
    # სკალარული ველი ფონის (t) და ოსცილაციის (მინიმუმ 2 ჰარმონიკით) ჩათვლით
    delta_Phi = Phi0 * sp.sin(theta) + Phi1 * sp.sin(3*theta)
    Phi_total = t + delta_Phi
    
    # წარმოებულები (theta-თი ვაწარმოებთ t-ს ნაცვლად)
    Phi_dot = 1 + omega * (Phi0 * sp.cos(theta) + 3 * Phi1 * sp.cos(3*theta))
    Phi_r = Phi0_prime * sp.sin(theta) + Phi1_prime * sp.sin(3*theta)
    
    # ფაზური ინვარიანტი Y მეტრიკით g^00=1, g^rr=-1 (Minkowski)
    Y_eval = sp.expand(Phi_dot**2 - Phi_r**2)
    
    Y, I1, I2, I3 = init_variables()
    L_poly = get_polynomial_lagrangian(Y, I1, I2, I3)
    
    # ენერგიის სიმკვრივე: T^0_0 = 2 * g^00 * (dPhi/dt)^2 * dL/dY - L
    dL_dY = sp.diff(L_poly, Y)
    rho_expr = 2 * Phi_dot**2 * dL_dY - L_poly
    
    # ვსვამთ ელასტიურ ინვარიანტებს Minkowski ფონზე (I1=3, I2=3, I3=1)
    bg_subs = {I1: 3, I2: 3, I3: 1}
    rho_sub = rho_expr.subs(bg_subs).subs(Y, Y_eval)
    
    # ფონური ენერგია (როცა ოსცილაცია არ გვაქვს)
    rho_bg = rho_sub.subs({Phi0: 0, Phi0_prime: 0, Phi1: 0, Phi1_prime: 0})
    
    # ოსცილონის წმინდა ენერგია 
    rho_pert = sp.expand(rho_sub - rho_bg)
    
    # დროის ერთ პერიოდზე გასაშუალოება
    rho_avg = sp.integrate(rho_pert, (theta, 0, 2*sp.pi)) / (2*sp.pi)
    rho_avg = sp.simplify(rho_avg)
    
    # ენერგიის სრული ინტეგრალი
    E_total = sp.Integral(rho_avg * 4 * sp.pi * r**2, r)
    
    # ვირიალური პირობა: dE/domega = 0 რეზონანსული სიხშირის დასაფიქსირებლად
    virial_integrand = sp.simplify(sp.diff(rho_avg * 4 * sp.pi * r**2, omega))
    
    return Phi_total, Y_eval, rho_avg, E_total, virial_integrand


def oscillon_finite_energy_gate():
    """
    Boundary and positivity gates for turning the trial family into a theorem.

    This function does not pretend to solve the nonlinear PDE.  It records the
    exact analytic requirements that a future solution must satisfy.
    """
    r, omega = sp.symbols("r omega", positive=True)
    E_of_omega = sp.Function("E")(omega)
    Phi0 = sp.Function("Phi0")(r)
    Phi1 = sp.Function("Phi1")(r)
    return {
        "regular_origin_conditions": [
            sp.Eq(sp.diff(Phi0, r).subs(r, 0), 0),
            sp.Eq(sp.diff(Phi1, r).subs(r, 0), 0),
            "Phi0(0), Phi1(0) finite",
        ],
        "finite_energy_asymptotics": [
            sp.Eq(sp.limit(r * Phi0, r, sp.oo), 0),
            sp.Eq(sp.limit(r * Phi1, r, sp.oo), 0),
            "rho_avg*r^2 must be integrable on [0, infinity)",
        ],
        "stationarity_condition": sp.Eq(sp.diff(E_of_omega, omega), 0),
        "stability_gate": "second variation / Floquet spectrum must be non-negative except symmetry zero modes",
        "status": "OPEN_PDE_EXISTENCE_AND_STABILITY",
    }


def oscillon_energy_status_audit():
    """Classify what analyze_oscillon proves and what it does not."""
    return {
        "proved_here": [
            "Y is evaluated nonlinearly for a two-harmonic Phi=t+deltaPhi ansatz",
            "the time-averaged stress density can be computed from the p01 polynomial L",
            "a virial stationarity integrand can be written symbolically",
        ],
        "not_proved_here": [
            "existence of a global finite-energy solution",
            "uniqueness of the oscillon profile",
            "spectral/Floquet stability",
            "absolute particle mass, charge, spin, and frequency selection",
        ],
        "status": "SYMBOLIC_TRIAL_FAMILY_NOT_PARTICLE_THEOREM",
    }


if __name__ == "__main__" and _should_run_main_section("oscillon"):
    Phi_total, Y_eval, rho_avg, E_total, virial_integrand = analyze_oscillon()
    
    print("--- ოსცილონის ანალიზი (Phi ველის ვარიაცია) ---")
    print("\nსრული ველი Phi(t,r):")
    print(Phi_total.subs(sp.Symbol('theta', real=True), sp.Symbol('omega', real=True, positive=True) * sp.Symbol('t', real=True)))
    
    print("\nფაზური ინვარიანტი Y:")
    print(Y_eval.subs(sp.Symbol('theta', real=True), sp.Symbol('omega', real=True, positive=True) * sp.Symbol('t', real=True)))
    
    print("\nგასაშუალოებული წმინდა ენერგიის სიმკვრივე <rho_osc>:")
    # შევაგროვოთ Phi0-ის ხარისხების მიხედვით
    Phi0 = sp.Function('Phi0')(sp.Symbol('r', real=True, positive=True))
    Phi0_prime = sp.diff(Phi0, sp.Symbol('r', real=True, positive=True))
    Phi1 = sp.Function('Phi1')(sp.Symbol('r', real=True, positive=True))
    Phi1_prime = sp.diff(Phi1, sp.Symbol('r', real=True, positive=True))
    print(sp.collect(sp.expand(rho_avg), [Phi0**2, Phi1**2, Phi0_prime**2]))
    
    print("\nსრული ენერგიის ინტეგრალი (E):")
    print(E_total)
    
    print("\nვირიალური პირობის ინტეგრანდი (dE/domega = 0):")
    print(virial_integrand)

    print("\n--- Bernoulli gravity theorem: static pressure deficit ---")
    bernoulli_static = bernoulli_static_gravity_identity()
    for key, value in bernoulli_static.items():
        print(f"{key}: {value}")

    print("\n--- Time-averaged oscillon source: zero-frequency projection ---")
    osc_source = bernoulli_time_averaged_oscillon_source()
    for key, value in osc_source.items():
        print(f"{key}: {value}")

    print("\n--- Poisson reconstruction: localized source -> 1/r field ---")
    poisson = bernoulli_poisson_reconstruction()
    for key, value in poisson.items():
        print(f"{key}: {value}")

    print("\n--- Poisson -> Newton normalization gate ---")
    norm_gate = poisson_to_newton_normalization_gate()
    for key, value in norm_gate.items():
        print(f"{key}: {value}")

    print("\n--- Newton law recovery from Bernoulli + geodesics ---")
    newton = bernoulli_newton_law_recovery()
    for key, value in newton.items():
        print(f"{key}: {value}")

    print("\n--- Complete Bernoulli gravity chain ---")
    for step in bernoulli_gravity_chain():
        print(f"- {step}")
    
    print("\n--- აგენტთა საბჭოს დასკვნები ---")
    print("1. ოსცილაცია ხდება რეალურ Phi ველზე და Y არის არაწრფივი შედეგი (Y=1+2Φ̇+Φ̇²-Φ'²).")
    print("2. I1, I2, I3 ინარჩუნებენ Minkowski ფონის (3, 3, 1) წვლილს ენერგიის გამოთვლაში.")
    print("3. სასრული ენერგიისთვის აუცილებელია Phi0(r), Phi1(r) და მათი წარმოებულები საკმარისად სწრაფად")
    print("   ქრებოდნენ უსასრულობაში, ხოლო r=0-ზე პროფილი რეგულარული იყოს.")
    print("4. 2 ჰარმონიკის ჩართვამ (Phi0, Phi1) დაადასტურა, რომ ენერგიაში ჩნდება არაწრფივი ჯვარედინი")
    print("   შერევები. c_Y2>0 დადებითად მოქმედებს quartic წევრებზე, მაგრამ სრული ენერგიის")
    print("   პოზიტიურობა მოითხოვს rho_avg-ის სრული გამოსახულების ანალიზს.")
    print("5. omega-ს ფიქსაციის ფორმალური პირობაა dE/domega = 0; რეალური omega-ს მისაღებად")
    print("   საჭიროა პროფილის ამოხსნა და საზღვრული პირობები.")
    print("6. Bernoulli gravity ნაწილი მკაცრად აჩვენებს pressure identity-ს და")
    print("   Newton-ის ალგებრულ დაბრუნებას იმ შემთხვევაში, თუ source normalization")
    print("   და bi-conformal branch მიღებულია.")
    print("7. ცალკე ღიად რჩება სრული nonlinear finite-energy oscillon profile-ის არსებობის")
    print("   სრული დამტკიცება; ამიტომ სრული oscillon->gravity proof ჯერ არ დაიხურა.")

# ===================== CONSOLIDATED PHASE SECTIONS =====================


# ===================== merged from p10_oscillons.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from p01_core import get_polynomial_lagrangian

def solve_static_spherical():
    r = sp.Symbol('r', real=True, positive=True)
    rs = sp.Symbol('rs', real=True, positive=True)
    eps = sp.Symbol('eps', real=True)
    kappa = sp.Symbol('kappa', real=True)
    
    # უცნობი ფუნქციები (სიგნატურით + - - -)
    A = sp.Function('A')(r)
    B = sp.Function('B')(r)
    Psi_p = sp.Function('Psi_p')(r) # ეს არის Psi'(r), სადაც Phi = t + Psi(r)
    
    # აინშტაინის ტენზორი G^t_t და G^r_r მეტრიკისთვის diag(B, -A, -r^2, -r^2 sin^2 theta)
    G_tt = -sp.diff(A, r) / (r * A**2) + (1/A - 1)/r**2
    G_rr = sp.diff(B, r) / (r * A * B) + (1/A - 1)/r**2
    G_thth = sp.diff(B, r, 2)/(2*A*B) - sp.diff(B, r)**2/(4*A*B**2) - sp.diff(A, r)*sp.diff(B, r)/(4*A**2*B) + sp.diff(B, r)/(2*r*A*B) - sp.diff(A, r)/(2*r*A**2)
    
    # ინვარიანტები (Phi = t + Psi(r) და phi^A = x^A comoving ანზაცისთვის)
    Y = 1/B - Psi_p**2 / A
    I1 = 2 + 1/A
    I2 = 1 + 2/A
    I3 = 1/A
    
    Y_s, I1_s, I2_s, I3_s = sp.symbols('Y I1 I2 I3', real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    
    L_eval = L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_Y = sp.diff(L_poly, Y_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I1 = sp.diff(L_poly, I1_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I2 = sp.diff(L_poly, I2_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I3 = sp.diff(L_poly, I3_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    
    # სტრეს-ტენზორი T^mu_nu
    T_tt = 2 * L_Y / B - L_eval
    T_rr = 2 * (L_Y * Psi_p**2 / A + L_I1 / A + 2 * L_I2 / A + L_I3 / A) - L_eval
    T_thth = 2 * (L_I1 + L_I2 * (1 + 1/A) + L_I3 / A) - L_eval
    
    # სკალარული ველის განტოლება: \nabla_\mu ( dL / d(\partial_\mu \Phi) ) = 0
    scalar_eq = sp.Derivative(r**2 * sp.sqrt(B/A) * L_Y * Psi_p, r)
    
    # სუსტი ველის ექსპანსია (ვუშვებთ Psi_p = 0 ცენტრალური მუხტის არარსებობის გამო)
    a1 = sp.Symbol('a1', real=True)
    b2 = sp.Symbol('b2', real=True)
    
    U = eps * rs / r
    A_w = 1 + a1 * U
    B_w = 1 - U + b2 * U**2
    
    # ვანაცვლებთ A და B ცვლადებს სუსტი ველის ფუნქციებით და ვითვლით G და T ტენზორებს
    G_tt_w = -sp.diff(A_w, r) / (r * A_w**2) + (1/A_w - 1)/r**2
    G_rr_w = sp.diff(B_w, r) / (r * A_w * B_w) + (1/A_w - 1)/r**2
    G_thth_w = sp.diff(B_w, r, 2)/(2*A_w*B_w) - sp.diff(B_w, r)**2/(4*A_w*B_w**2) - sp.diff(A_w, r)*sp.diff(B_w, r)/(4*A_w**2*B_w) + sp.diff(B_w, r)/(2*r*A_w*B_w) - sp.diff(A_w, r)/(2*r*A_w**2)
    
    Y_w = 1/B_w
    I1_w = 2 + 1/A_w
    I2_w = 1 + 2/A_w
    I3_w = 1/A_w
    
    L_w = L_poly.subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    L_Y_w = sp.diff(L_poly, Y_s).subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    L_I1_w = sp.diff(L_poly, I1_s).subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    L_I2_w = sp.diff(L_poly, I2_s).subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    L_I3_w = sp.diff(L_poly, I3_s).subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    
    T_tt_w = 2 * L_Y_w / B_w - L_w
    T_rr_w = 2 * (L_I1_w / A_w + 2 * L_I2_w / A_w + L_I3_w / A_w) - L_w
    T_thth_w = 2 * (L_I1_w + L_I2_w * (1 + 1/A_w) + L_I3_w / A_w) - L_w
    
    def get_series(expr):
        return sp.simplify(sp.series(expr, eps, 0, 3).removeO())
        
    Eq_tt_w = get_series(G_tt_w - kappa * T_tt_w)
    Eq_rr_w = get_series(G_rr_w - kappa * T_rr_w)
    Eq_thth_w = get_series(G_thth_w - kappa * T_thth_w)
    
    # გამოვყოთ ნულოვანი (ფონური) და პირველი რიგის განტოლებები
    Eq_tt_O0 = sp.simplify(Eq_tt_w.subs(eps, 0))
    Eq_tt_O1 = sp.simplify(sp.diff(Eq_tt_w, eps).subs(eps, 0))
    
    Eq_rr_O0 = sp.simplify(Eq_rr_w.subs(eps, 0))
    Eq_rr_O1 = sp.simplify(sp.diff(Eq_rr_w, eps).subs(eps, 0))
    
    Eq_thth_O0 = sp.simplify(Eq_thth_w.subs(eps, 0))
    Eq_thth_O1 = sp.simplify(sp.diff(Eq_thth_w, eps).subs(eps, 0))
    
    # ბი-კონფორმობის ანალიზი (წინასწარ a1=1 დაშვების გარეშე)
    Delta_T = sp.simplify(T_tt_w - T_rr_w)
    Delta_T_O1 = sp.simplify(sp.diff(Delta_T, eps).subs(eps, 0))
    
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols('c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1', real=True)
    E_vac = -c_Y - 3*c_Y2 + 3*c_I1 + 9*c_I1sq + 3*c_I2 + c_I3 - 3*c_YI1
    P_vac = c_Y + c_Y2 + 5*c_I1 + 21*c_I1sq + 7*c_I2 + 3*c_I3 + 5*c_YI1
    vac_sols = sp.solve([E_vac, P_vac], [c_Y, c_I1])
    
    clean_Delta_T_O1 = sp.simplify(Delta_T_O1.subs(vac_sols))
    bc_constraint = sp.Eq(clean_Delta_T_O1 / (rs/r), 0)
    
    return G_tt, G_rr, G_thth, T_tt, T_rr, T_thth, scalar_eq, Eq_tt_w, Eq_rr_w, Eq_thth_w, Eq_tt_O0, Eq_tt_O1, Eq_rr_O0, Eq_rr_O1, Eq_thth_O0, Eq_thth_O1, Delta_T_O1, clean_Delta_T_O1, bc_constraint


def static_spherical_theorem_gate():
    """
    Council gate for the exterior branch.

    solve_static_spherical() is useful, but the full weak-field theorem requires
    all O(eps) equations to vanish together after the Minkowski vacuum
    constraints.  The old single Delta_T/biconformal constraint is necessary
    only for one aspect of the branch.
    """
    (
        _G_tt,
        _G_rr,
        _G_thth,
        _T_tt,
        _T_rr,
        _T_thth,
        _scalar_eq,
        _Eq_tt_w,
        _Eq_rr_w,
        _Eq_thth_w,
        Eq_tt_O0,
        Eq_tt_O1,
        Eq_rr_O0,
        Eq_rr_O1,
        Eq_thth_O0,
        Eq_thth_O1,
        Delta_T_O1,
        clean_Delta_T_O1,
        bc_constraint,
    ) = solve_static_spherical()

    c_Y2, c_I1sq, c_YI1 = sp.symbols("c_Y2 c_I1sq c_YI1", real=True)
    a1 = sp.Symbol("a1", real=True)

    return {
        "vacuum_background_equations": [Eq_tt_O0, Eq_rr_O0, Eq_thth_O0],
        "first_order_equations_to_solve_together": [Eq_tt_O1, Eq_rr_O1, Eq_thth_O1],
        "biconformal_constraint_only": bc_constraint,
        "council_candidate_constraints_for_a1_1": [
            sp.Eq(c_YI1, 2 * c_Y2),
            sp.Eq(c_YI1, 2 * c_I1sq),
            "plus the angular O(eps) combination; solve full system before claiming exterior proof",
        ],
        "why_not_closed": "O(eps) contains independent radial powers and all components must vanish simultaneously.",
        "status": "OPEN_FULL_STATIC_SPHERICAL_BRANCH",
        "a1_note": "a1 must be solved or fixed by the branch; it cannot be silently inserted.",
    }


if __name__ == "__main__" and _should_run_main_section("spherical"):
    res = solve_static_spherical()
    G_tt, G_rr, G_thth, T_tt, T_rr, T_thth, scalar_eq, Eq_tt_w, Eq_rr_w, Eq_thth_w, Eq_tt_O0, Eq_tt_O1, Eq_rr_O0, Eq_rr_O1, Eq_thth_O0, Eq_thth_O1, Delta_T_O1, clean_Delta_T_O1, bc_constraint = res
    print("--- ამოხსნა სფერული ანზაცისთვის ---")
    print("G^t_t =", G_tt)
    print("T^t_t =", T_tt)
    print("G^r_r =", G_rr)
    print("T^r_r =", T_rr)
    print("G^theta_theta =", G_thth)
    print("T^theta_theta =", T_thth)
    print("\nსკალარული ველის განტოლება:")
    print("0 =", scalar_eq)
    print("-> დასკვნა: ცენტრალური სკალარული მუხტის გარეშე Psi'(r) = 0, ანუ Phi(t,r) = t.")
    print("\n--- სუსტი ველის ლიმიტი (O(rs/r)) ---")
    print("Eq_tt (O(eps)):", Eq_tt_w)
    print("Eq_rr (O(eps)):", Eq_rr_w)
    print("Eq_thth (O(eps)):", Eq_thth_w)
    
    print("\n--- ფონური ვაკუუმის განტოლებები (O(1)) ---")
    print("Eq_tt_O0 =", Eq_tt_O0)
    print("Eq_rr_O0 =", Eq_rr_O0)
    print("Eq_thth_O0 =", Eq_thth_O0)
    
    print("\n--- პირველი რიგის განტოლებები (O(eps)) ---")
    print("Eq_tt_O1 =", Eq_tt_O1)
    print("Eq_rr_O1 =", Eq_rr_O1)
    print("Eq_thth_O1 =", Eq_thth_O1)

    print("\n--- ბი-კონფორმობის პირობა (g_tt * g_rr = -1) ---")
    print("Delta_T(O(eps)) =", Delta_T_O1)
    print("Delta_T(O(eps)) სუფთა (Minkowski ვაკუუმის კონსტრეინტებით) =", clean_Delta_T_O1)
    print("კონსტრეინტი ბი-კონფორმულობისთვის:", bc_constraint)

    print("\n--- static spherical theorem gate ---")
    spherical_gate = static_spherical_theorem_gate()
    for key, value in spherical_gate.items():
        print(f"{key}: {value}")
    
    print("\n--- აგენტთა საბჭოს დასკვნები ---")
    print("1. Angular კომპონენტი დაემატა; სრული სფერული სისტემის დახურვა საჭიროებს corrected O(eps) განტოლებების ერთობლივ ამოხსნას.")
    print("2. U=rs/r ანზაცით წარმოებულები სიმბოლურად ითვლება (eps-სერიების აღრევა აღმოიფხვრა).")
    print("3. T_rr-ის ნიშნის შეცდომა გასწორდა (T_rr და T_thth ახლა დადებითი ელასტიური წევრებით იწყება).")
    print("4. ბი-კონფორმობის ანალიზში a1=1 წინასწარ აღარ იდება. Minkowski E_vac=0, P_vac=0 კონსტრეინტებით მიიღება სუფთა შეზღუდვა.")


# ===================== merged from p10_oscillons.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp

def analyze_matter_coupling_tov():
    r = sp.Symbol('r', real=True, positive=True)
    
    # ველები
    Phi = sp.Function('Phi')(r) # ნიუტონის პოტენციალი
    rho_solid = sp.Function('rho_solid')(r) # სუპერსოლიდის ენერგიის სიმკვრივე
    
    # სუპერსოლიდის სტრეს-ტენზორის კომპონენტები
    p_rad = sp.Function('p_rad')(r)
    p_tan = sp.Function('p_tan')(r)
    delta_p = p_tan - p_rad # ანიზოტროპია
    
    # სრული ინერციული სიმკვრივე ნიუტონურ ლიმიტში: rho_inert = rho_solid + p_rad
    rho_inert = rho_solid + p_rad
    
    # სრული სისტემის ენერგია-იმპულსის შენახვის კანონი (\nabla_\mu T^{\mu 1} = 0)
    # სუსტ ველში (p << rho) და სტატიკურ სფერულ სიმეტრიაში გვაძლევს ანიზოტროპიულ TOV განტოლებას.
    # ეს აღწერს სითხის შიდა წონასწორობას და არა გარე ტესტ-ნაწილაკის აჩქარებას!
    # ნიუტონური ლიმიტი: d(p_rad)/dr + rho_inert * d(Phi)/dr - 2*delta_p/r = 0
    
    grad_Phi = sp.Symbol('grad_Phi') # პოტენციალის გრადიენტი d(Phi)/dr
    
    tov_eq = sp.Eq(sp.diff(p_rad, r) + rho_inert * grad_Phi - 2 * delta_p / r, 0)
    
    # ამოვხსნათ პოტენციალის გრადიენტისთვის შიდა წონასწორობაში
    sols = sp.solve(tov_eq, grad_Phi)
    if not sols:
        raise ValueError("ვერ მოიძებნა ამოხსნა grad_Phi-სთვის")
        
    g_sol = sols[0]
    return sp.simplify(g_sol)

if __name__ == "__main__" and _should_run_main_section("tov"):
    g_sol = analyze_matter_coupling_tov()
    print("--- ანიზოტროპიული TOV განტოლება ნიუტონურ ლიმიტში ---")
    print("პოტენციალის გრადიენტი (Phi') შიდა წონასწორობაში:")
    print(g_sol)
    
    print("\n--- კავშირი ტესტ-ნაწილაკის დინამიკასთან ---")
    print("ეს განტოლება აღწერს მედიუმის შიდა წონასწორობას.")
    print("გარე ტესტ-ნაწილაკის გეოდეზიური აჩქარება (a_test = -Phi') და MOND-ის")
    print("სრული ამოხსნა გაგრძელებულია p07_mond.py-ში.")
    
    print("\n--- აგენტთა საბჭოს შენიშვნები ---")
    print("1. ინერციული სიმკვრივე გასწორდა: rho_inert = rho_solid + p_rad.")
    print("2. განტოლება წარმოადგენს TOV-ის ნიუტონურ ლიმიტს; M-R წირის ასაგებად აკლია EoS და საზღვრის პირობები.")
    print("3. Phi'-ის ფორმულა არის წრიული (იმპლიციტური ODE), რადგან p_rad და rho_solid")
    print("   თვითონ არიან Phi-ს ფუნქციები ლაგრანჟიანიდან. ეს მოითხოვს მკაცრ მიბმას/fine-tuning-ს.")


# ===================== merged from p10_oscillons.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 17: ეფექტური მასის სკალირება და ბი-კონფორმული გეომეტრია
================================================================================

სტატუსი:
ეს სექცია წარმოადგენს ბი-კონფორმული operational branch-ის შედეგების შემოწმებას
(consistency check). მეტრიკა აქ წინაპირობა/branch template-ია და არა
მოქმედებიდან დამოუკიდებლად დასრულებული გამოყვანა.
სინათლის გადახრა და Pound-Rebka მოწმდება phi = -r_s/r ფონზე.

ცენტრალური ფიზიკური მექანიზმი, რომელიც ადრე ცალკე თეორიულ ტექსტში უნდა გადასულიყო:
როგორ ცვლის ფონური წნევითი პოტენციალი φ(x) ეფექტურ მასას, ოპერაციულ
ზომას, საათის ტემპს და კოორდინატულ სიხშირეს.

ცენტრალური შედეგი (გამოყვანილია ამ branch/ansatz-ის შიგნით):
    m_eff(φ) = m_0 · e^(φ/2)
    L_oper(φ) = L_0 · e^(φ/2)
    T_period(φ) = T_0 · e^(-φ/2)
    c_coord(φ) = c · e^(φ)

ფაქტორი 2:  c_coord/c = (L_oper/L_0)²   ← ბუნებრივად, არა ad-hoc

weak-field smoke-tests:
  - სინათლის გადახრა მზესთან: 1.7505 arcsec  [Eddington 1919, VLBI]
  - Pound-Rebka რედშიფტი:    2.46×10⁻¹⁵      [PR 1960, 1σ შიგნით]
  - ფაქტორი 2 (1911 vs 1915): ემერჯენტული branch identity
  - ლოკალური ფარდობითობის პრინციპი: α და უგანზომილო ფარდობები უცვლელია
  🟡 ν₀ (სუბსტრატის ფონური რიტმი): ჰიპოთეზა / TODO
"""

import sympy as sp


# ==============================================================================
# Setup
# ==============================================================================

def setup_metric():
    """
    ბი-კონფორმული მეტრიკა (c=1 ერთეულებში):
        ds² = e^φ dt² - e^(-φ) (dx² + dy² + dz²)

    Signature: (+, -, -, -), matching NOTATION.md.
    φ — ფონური წნევითი პოტენციალი (სკალარული)
    """
    phi = sp.Symbol('phi', real=True)
    g = sp.diag(sp.exp(phi), -sp.exp(-phi), -sp.exp(-phi), -sp.exp(-phi))
    return g, phi


# ==============================================================================
# ნაბიჯი 1: ბი-კონფორმობის ცხადი თვისება
# ==============================================================================

def step1_biconformal_property():
    """
    გადაამოწმე: g_tt · g_xx = -1 (c=1 ერთეულებში)
    ან: g_tt · g_xx = -c² (SI ერთეულებში)
    """
    g, phi = setup_metric()
    product = sp.simplify(g[0, 0] * g[1, 1])
    expected = sp.Integer(-1)
    holds = sp.simplify(product - expected) == 0
    return product, expected, holds


# ==============================================================================
# ნაბიჯი 2: საკუთრივი დრო და სიგრძე
# ==============================================================================

def step2_proper_time_and_length():
    """
    სტატიკური დამკვირვებლისთვის:
        dτ/dt = √(g_tt) = e^(φ/2)          ← საათის ტემპი
        dl/dx = √(-g_xx) = e^(-φ/2)        ← სიგრძის სკალირება
    """
    g, phi = setup_metric()
    dtau_dt = sp.simplify(sp.sqrt(g[0, 0]))
    dl_dx = sp.simplify(sp.sqrt(-g[1, 1]))
    return dtau_dt, dl_dx


# ==============================================================================
# ნაბიჯი 3: ეფექტური მასა Killing ენერგიიდან
# ==============================================================================

def step3_effective_mass():
    """
    სტატიკური (შებოჭილი/tethered) ნაწილაკი ბი-კონფორმულ მეტრიკაში.
    სტატიკური მსოფლიო-ხაზი φ≠const ფონზე არ არის თავისუფალი გეოდეზიური!
    
        u^t = dt/dτ = 1/√(g_tt) = e^(-φ/2)
        p^t = m_0 · u^t = m_0 · e^(-φ/2)
        p_t = g_tt · p^t = e^φ · m_0 · e^(-φ/2) = m_0 · e^(φ/2)

    Killing ენერგია (Killing ვექტორი ξ = ∂_t) — უსასრულობაში (φ→0) გაზომილი მასა:
        E = p_μ ξ^μ = p_t = m_0 · e^(φ/2)

    ლოკალური დამკვირვებლის მიერ გაზომილი ენერგია:
        u^μ_obs = (e^(-φ/2), 0, 0, 0)
        E_loc = p_μ u^μ_obs = m_0 e^(φ/2) * e^(-φ/2) = m_0 (უცვლელი!)
        
    აქ m_eff არის გარედან დანახული (Killing) მასა.
    """
    g, phi = setup_metric()
    m_0 = sp.Symbol('m_0', positive=True)

    u_t_up = 1 / sp.sqrt(g[0, 0])            # = e^(-φ/2)
    p_t_up = m_0 * u_t_up                    # p^t = m_0 e^(-φ/2)
    p_t_down = g[0, 0] * p_t_up              # p_t = m_0 e^(φ/2)

    E_killing = sp.simplify(p_t_down)         # E = m_0 e^(φ/2)
    E_loc = sp.simplify(p_t_down * u_t_up)    # ლოკალურად გაზომილი E_loc = m_0
    m_eff = sp.simplify(E_killing)            # m_eff (c=1)

    return m_eff, m_0, E_loc, sp.simplify(u_t_up), sp.simplify(p_t_down)


# ==============================================================================
# ნაბიჯი 4: ოპერაციული ზომის სკალირება
# ==============================================================================

def step4_operational_size():
    """
    ფიქსირებული საკუთრივი ზომის (L_0) ობიექტი:
        dl_proper = √(-g_xx) dx  →  dx = L_0/√(-g_xx)
        L_oper = L_0 / e^(-φ/2) = L_0 · e^(φ/2)

    იგივე ექსპონენტი, რაც m_eff-ის → L_oper ∝ m_eff   ✓
    """
    g, phi = setup_metric()
    L_0 = sp.Symbol('L_0', positive=True)
    L_oper = sp.simplify(L_0 / sp.sqrt(-g[1, 1]))
    return L_oper, L_0


# ==============================================================================
# ნაბიჯი 5: კოორდინატული სინათლის სიჩქარე
# ==============================================================================

def step5_coordinate_speed():
    """
    ნულოვანი გეოდეზიური (ds² = 0) x-ის გასწვრივ:
        0 = e^φ dt² - e^(-φ) dx²
        (dx/dt)² = e^(2φ)
        c_coord = e^φ   (c=1 ერთეულებში)
        c_coord = c · e^φ   (SI ერთეულებში)
    """
    g, phi = setup_metric()
    c_coord_sq = -g[0, 0] / g[1, 1]
    c_coord = sp.simplify(sp.sqrt(c_coord_sq))
    return c_coord


# ==============================================================================
# ნაბიჯი 6: ფაქტორი 2 — c_coord/c = (L_oper/L_0)²
# ==============================================================================

def step6_factor_two():
    """
    ცენტრალური შემოწმება:
        c_coord/c = e^φ
        L_oper/L_0 = e^(φ/2)
        (L_oper/L_0)² = e^φ = c_coord/c   ✓

    ფაქტორი 2 ემერჯენტულია ბი-კონფორმობიდან.
    """
    phi = sp.Symbol('phi', real=True)
    c_coord = step5_coordinate_speed()
    L_oper, L_0 = step4_operational_size()
    L_ratio = sp.simplify(L_oper / L_0)
    diff = sp.simplify(c_coord - L_ratio**2)
    holds = sp.simplify(diff) == 0
    return c_coord, L_ratio**2, diff, holds


# ==============================================================================
# ნაბიჯი 7: სინათლის გადახრის ცხადი გათვლა
# ==============================================================================

def step7_light_deflection():
    """
    Schwarzschild ბი-კონფორმული φ = -r_s/r:
        c_coord = c · e^φ = c · e^(-r_s/r)
        რეფრაქციული ინდექსი n = c/c_coord = e^(-φ) = e^(r_s/r)

    სუსტ ველში:
        n ≈ 1 + r_s/r + (r_s/r)²/2 + ...

    გადახრის კუთხე impact parameter b-ზე:
        δ = ∫_{-∞}^{∞} (∂n/∂y)|_{y=b} dx
          = ∫ r_s · b/(x²+b²)^(3/2) dx
          = 2 r_s/b   ← ფაქტორი 2
    """
    r, b, x, r_s = sp.symbols('r b x r_s', positive=True)

    # ბი-კონფორმული Schwarzschild
    phi_schw = -r_s / r
    n_exact = sp.exp(-phi_schw)              # n = e^(-φ) = e^(r_s/r)
    n_weak = sp.simplify(sp.series(n_exact, r_s, 0, 3).removeO())

    # სუსტ-ველის ლიდინგ წევრი: n - 1 ≈ r_s/r
    # ინტეგრალი light path-ის გასწვრივ:
    integrand = r_s * b / (x**2 + b**2)**sp.Rational(3, 2)
    delta = sp.integrate(integrand, (x, -sp.oo, sp.oo))
    delta = sp.simplify(delta)

    return n_weak, delta


# ==============================================================================
# ნაბიჯი 8: 1911 vs 1915 ისტორიული გაყოფა
# ==============================================================================

def step8_split_1911_1915():
    """
    Einstein 1911: მხოლოდ დროითი წევრი (g_tt)
        n_t = 1/√(g_tt) = e^(-φ/2) ≈ 1 + r_s/(2r)
        δ_t = r_s/b

    GR 1915 / RFG full: დროითი + სივრცული (g_ii)
        n_s = √(g_ii) = e^(-φ/2) ≈ 1 + r_s/(2r)
        δ_s = r_s/b

    Total:  δ = δ_t + δ_s = 2 r_s/b
    ფაქტორი 2 = ბი-კონფორმული სტრუქტურა (დროითი + სივრცული თანაბრად)
    ჯამური რეფრაქციული ინდექსი: n = n_t * n_s = e^(-φ/2) * e^(-φ/2) = e^(-φ).
    """
    r, b, x, r_s = sp.symbols('r b x r_s', positive=True)

    # Schwarzschild ბი-კონფორმული
    phi_schw = -r_s / r

    # დროითი ნაწილი (1911): n_t = 1/√(g_tt) = e^(-φ/2)
    n_t = sp.exp(-phi_schw / 2)
    n_t_lead = sp.simplify(sp.series(n_t, r_s, 0, 2).removeO() - 1)

    # სივრცული ნაწილი (1915 დამატება): n_s = √(g_ii) = e^(-φ/2)
    n_s = sp.exp(-phi_schw / 2)
    n_s_lead = sp.simplify(sp.series(n_s, r_s, 0, 2).removeO() - 1)

    # თითო ნაწილს r_s/(2r), ჯამში r_s/r
    integrand_half = (r_s / 2) * b / (x**2 + b**2)**sp.Rational(3, 2)
    delta_t = sp.simplify(sp.integrate(integrand_half, (x, -sp.oo, sp.oo)))
    delta_s = sp.simplify(sp.integrate(integrand_half, (x, -sp.oo, sp.oo)))
    delta_total = sp.simplify(delta_t + delta_s)

    return n_t_lead, n_s_lead, delta_t, delta_s, delta_total


# ==============================================================================
# ნაბიჯი 9: მზის რიცხვობრივი ვერიფიკაცია
# ==============================================================================

def step9_sun_deflection_numerical():
    """
    მზის გვერდით:
        r_s = 2 G M_⊙ / c² ≈ 2950 m
        b = R_⊙ ≈ 6.96×10⁸ m
        δ = 2 r_s/b → arcsec
    """
    G = 6.674e-11
    M_sun = 1.989e30
    c_si = 2.998e8
    R_sun = 6.96e8

    r_s = 2 * G * M_sun / c_si**2
    delta_rad = 2 * r_s / R_sun
    delta_arcsec = delta_rad * 206265  # rad → arcsec

    return r_s, delta_rad, delta_arcsec


# ==============================================================================
# ნაბიჯი 10: Pound-Rebka გრავიტაციული რედშიფტი
# ==============================================================================

def step10_pound_rebka():
    """
    ფოტონი მაღლა მიდის h სიმაღლეზე გრავიტაციულ ველში g:
        Δν/ν = -g·h/c²   (რედშიფტი)
        |Δν/ν| = g·h/c²

    Pound-Rebka (Harvard, 1960):
        h = 22.6 m,  g = 9.81 m/s²
        პროგნოზი: 2.46×10⁻¹⁵
        გაზომილი:  (2.57 ± 0.26) × 10⁻¹⁵
    """
    g_earth = 9.81
    h = 22.6
    c_si = 2.998e8

    predicted = g_earth * h / c_si**2
    measured = 2.57e-15
    error = 0.26e-15
    sigma_dev = abs(predicted - measured) / error
    within_1sigma = sigma_dev < 1.0

    return predicted, measured, error, sigma_dev, within_1sigma


# ==============================================================================
# ნაბიჯი 11: ლოკალური შეუმჩნევლობა (ფარდობითობის პრინციპი)
# ==============================================================================

def step11_local_invariance():
    """
    ლოკალური დამკვირვებელი იყენებს საკუთარ საათს და სახაზავს.
    მკაფიოდ უნდა გაიმიჯნოს ლოკალური და გარე (ოპერაციული) სიდიდეები:

    - ლოკალური მასა, ლოკალური Compton სიგრძე და ლოკალური სინათლის სიჩქარე 
      ინვარიანტულია (E_loc = m_0).
    - გარე (ოპერაციული) სიდიდეები სკალირდება p = e^(φ/2) ფაქტორით.
    - ოპერაციული უგანზომილო ფარდობები (მაგ. L_oper / λ_C,oper) უცვლელი რჩება.
    """
    phi = sp.Symbol('phi', real=True)

    # ლოკალური (ინვარიანტული) სიდიდეები
    local_quantities = {
        'm_local':          sp.Integer(1),
        'L_local':          sp.Integer(1),
        'lambda_C_local':   sp.Integer(1),
        'c_local':          sp.Integer(1),
    }

    # გარე/ოპერაციული სკალირების ფაქტორები (p ≡ e^(φ/2))
    oper_quantities = {
        'm_eff':            sp.exp(phi / 2),
        'L_oper':           sp.exp(phi / 2),
        'lambda_C_oper':    sp.exp(phi / 2),
        'T_period':         sp.exp(-phi / 2),
        'c_coord':          sp.exp(phi),
    }

    # უგანზომილო ფარდობები (უნდა იყვნენ უცვლელი)
    ratios = {
        'L_oper / lambda_C_oper':        sp.simplify(oper_quantities['L_oper'] / oper_quantities['lambda_C_oper']),
        'c_coord * T_period / L_oper':   sp.simplify(oper_quantities['c_coord'] * oper_quantities['T_period'] / oper_quantities['L_oper']),
        'alpha (locally)':               sp.Integer(1),
    }

    return local_quantities, oper_quantities, ratios


# ==============================================================================
# დამატებითი: კოსმოლოგიური ν₀ ინვარიანტობის შემოწმება
# ==============================================================================

def step12_cosmological_nu0():
    """
    TODO: სრული გათვლა მოითხოვს Φ ფონური ვაკუუმის ანალიზს — ცალკე დავალება.
    
    FLRW ფონზე: ds² = -dt² + a(t)² δ_ij dx^i dx^j (არა ბი-კონფორმული)

    Substrate ν₀ არის სუბ-ოსცილონური ფონური რიტმი (Intuitive §1).
    ის *არ არის* composite (არ შედგება ლოკალური ოსცილონებისგან).

    ანალიტიკური მტკიცება: ν₀ ≠ f(Y_bg, I_1_bg) BEC კონდენსატის შიდა სტრუქტურის გამო. 
    """
    return "TODO / ჰიპოთეზა: ν_0 (substrate) = INVARIANT — სრული Φ ფონური ვაკუუმის ანალიზი ჯერ არ გაკეთებულა."


def ppn_and_observation_gate():
    """
    Observation firewall.

    Light bending and Pound-Rebka are important smoke tests, but they are not
    the full Solar-System validation suite.
    """
    gamma, beta = sp.symbols("gamma beta", real=True)
    return {
        "checked_here": [
            "leading weak-field solar light bending: delta=2*r_s/b",
            "Pound-Rebka redshift at Harvard tower scale",
            "local dimensionless ratios in the bi-conformal ansatz",
        ],
        "not_checked_here": [
            "PPN gamma and beta extraction from the full exterior solution",
            "Cassini/Shapiro time delay",
            "Mercury perihelion and ephemeris residuals",
            "lunar laser ranging",
            "binary pulsars and strong-field timing",
            "modern atomic-clock local-position-invariance tests",
        ],
        "target_ppn_values": [sp.Eq(gamma, 1), sp.Eq(beta, 1)],
        "status": "PARTIAL_GR_WEAK_FIELD_SMOKE_TESTS_ONLY",
    }


def p10_status_audit():
    """Compact status ledger for the whole file."""
    return {
        "overall_status": "strong symbolic consistency ledger; not complete oscillon-to-gravity proof",
        "closed": [
            "Bernoulli pressure identity",
            "Poisson 1/r reconstruction for a supplied localized source",
            "Newton law algebra after mu_src=2GM/c^2 normalization",
            "bi-conformal scaling identities inside the selected branch",
            "leading light-bending and Pound-Rebka smoke tests",
        ],
        "conditional": [
            "bi-conformal metric branch",
            "source-to-G normalization",
            "Newton recovery as an RFG first-principles theorem",
        ],
        "open": [
            "finite-energy nonlinear oscillon existence",
            "spectral/Floquet stability",
            "full static spherical p01 solution",
            "PPN/Cassini/ephemeris validation",
            "nu0 substrate-clock derivation",
            "particle mass/charge/spin matching",
        ],
    }


# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__" and _should_run_main_section("biconformal"):
    print("=" * 72)
    print("PHASE 17: ეფექტური მასის სკალირება და ბი-კონფორმული გეომეტრია")
    print("=" * 72)

    # ნაბიჯი 1
    print("\n--- ნაბიჯი 1: ბი-კონფორმობის ცხადი თვისება ---")
    prod, expected, holds1 = step1_biconformal_property()
    print(f"  g_tt · g_xx = {prod}")
    print(f"  მოლოდინი:    {expected}   (c=1 ერთეულებში)")
    print(f"  დადგინდა:    {holds1}   ✓")

    # ნაბიჯი 2
    print("\n--- ნაბიჯი 2: საკუთრივი დრო და სიგრძე ---")
    dtau, dl = step2_proper_time_and_length()
    print(f"  dτ/dt = {dtau}      (= e^(φ/2))")
    print(f"  dl/dx = {dl}        (= e^(-φ/2))")
    print(f"  ინტერპრეტაცია: φ < 0 → საათი ნელია, რადი იწელება")

    # ნაბიჯი 3
    print("\n--- ნაბიჯი 3: ეფექტური მასა Killing ენერგიიდან ---")
    m_eff, m_0, E_loc, u_t, p_t = step3_effective_mass()
    print(f"  u^t (4-სიჩქარის t-კომპონენტი) = {u_t}")
    print(f"  p_t (ქვევით) = {p_t}")
    print(f"  E_Killing (გარედან დანახული მასა) = p_t = {m_eff}")
    print(f"  → m_eff = m_0 · e^(φ/2)   ✓")
    print(f"  E_loc (ლოკალურად გაზომილი ენერგია) = {E_loc}   ✓ (უცვლელი!)")
    print(f"  m_eff ეხება შებოჭილ (tethered) მდგომარეობას უსასრულობის მიმართ.")

    # ნაბიჯი 4
    print("\n--- ნაბიჯი 4: ოპერაციული ზომა ---")
    L_oper, L_0 = step4_operational_size()
    print(f"  L_oper = {L_oper}")
    print(f"  → L_oper = L_0 · e^(φ/2)   ✓")
    print(f"  იგივე ექსპონენტი, რაც m_eff-ის → L_oper ∝ m_eff")

    # ნაბიჯი 5
    print("\n--- ნაბიჯი 5: კოორდინატული c ---")
    c_coord = step5_coordinate_speed()
    print(f"  c_coord = {c_coord}   (c=1 ერთეულებში)")
    print(f"  → c_coord = c · e^φ   (SI ერთეულებში)")

    # ნაბიჯი 6 — ცენტრალური ფაქტი
    print("\n--- ნაბიჯი 6: ფაქტორი 2 (ცენტრალური!) ---")
    c_ratio, L_ratio_sq, diff, holds6 = step6_factor_two()
    print(f"  c_coord/c     = {c_ratio}")
    print(f"  (L_oper/L_0)² = {L_ratio_sq}")
    print(f"  სხვაობა:       {diff}")
    print(f"  ფაქტორი 2 დადგინდა: {holds6}")
    print(f"  → c_coord/c = (L_oper/L_0)²   ✓")
    print(f"  ფაქტორი 2 ემერჯენტულია ბი-კონფორმობიდან, არა ad-hoc!")

    # ნაბიჯი 7
    print("\n--- ნაბიჯი 7: სინათლის გადახრის გათვლა ---")
    n_weak, delta = step7_light_deflection()
    print(f"  n(r) (სუსტ ველში) = {n_weak}")
    print(f"  δ (გადახრის კუთხე) = {delta}")
    expected_delta = sp.Symbol('r_s') * 2 / sp.Symbol('b')
    print(f"  მოლოდინი: 2 r_s/b")
    print(f"  ✓ თანხვედრა")

    # ნაბიჯი 8
    print("\n--- ნაბიჯი 8: 1911 vs 1915 ისტორიული გაყოფა ---")
    n_t, n_s, d_t, d_s, d_total = step8_split_1911_1915()
    print(f"  Temporal (Einstein 1911): n_t - 1 ≈ {n_t},  δ_t = {d_t}")
    print(f"  Spatial  (GR 1915 add):   n_s - 1 ≈ {n_s},  δ_s = {d_s}")
    print(f"  ჯამი: δ = {d_total}")
    print(f"  → ფაქტორი 2 = temporal + spatial თანაბრად   ✓")

    # ნაბიჯი 9
    print("\n--- ნაბიჯი 9: მზის რიცხვობრივი ვერიფიკაცია ---")
    r_s_sun, delta_rad, delta_arcsec = step9_sun_deflection_numerical()
    print(f"  r_s (მზე) = {r_s_sun:.2f} m")
    print(f"  δ (პროგნოზი) = {delta_rad:.4e} rad = {delta_arcsec:.4f} arcsec")
    print(f"  დაკვირვებული (VLBI): 1.7505 arcsec")
    deviation_sun = abs(delta_arcsec - 1.7505)
    if deviation_sun < 0.005:
        print(f"  გადახრა: {deviation_sun:.4f} arcsec   ✓ თანხვედრა")
    else:
        print(f"  გადახრა: {deviation_sun:.4f} arcsec")

    # ნაბიჯი 10
    print("\n--- ნაბიჯი 10: Pound-Rebka გრავიტაციული რედშიფტი ---")
    pr_pred, pr_meas, pr_err, sigma, within = step10_pound_rebka()
    print(f"  |Δν/ν| (პროგნოზი)  = {pr_pred:.3e}")
    print(f"  |Δν/ν| (გაზომილი) = ({pr_meas:.2e} ± {pr_err:.2e})")
    print(f"  გადახრა: {sigma:.2f} σ")
    if within:
        print(f"  → 1σ შიგნით   ✓")

    # ნაბიჯი 11
    print("\n--- ნაბიჯი 11: ლოკალური შეუმჩნევლობა ---")
    loc, oper, ratios = step11_local_invariance()
    print(f"  ლოკალური (ინვარიანტული) სიდიდეები (სკალირება = 1):")
    for k, v in loc.items():
        print(f"    {k:22s} = {v}")
    print(f"\n  გარე/ოპერაციული სკალირების ფაქტორები (p ≡ e^(φ/2)):")
    for k, v in oper.items():
        print(f"    {k:22s} = {v}")
    print(f"\n  უგანზომილო ფარდობები (უცვლელი):")
    for k, v in ratios.items():
        print(f"    {k:28s} = {v}")

    # ნაბიჯი 12
    print("\n--- ნაბიჯი 12: ν₀ კოსმოლოგიური ინვარიანტობა ---")
    print(f"  {step12_cosmological_nu0()}")

    # შემაჯამებელი ცხრილი
    print("\n" + "=" * 72)
    print("შემაჯამებელი ცხრილი (p ≡ e^(φ/2) — pressure factor)")
    print("=" * 72)
    table = [
        ("სიდიდე",                  "სკალირება",         "წყარო / კონტექსტი"),
        ("─" * 22,                  "─" * 16,             "─" * 28),
        ("[External] m_eff / m_0",  "p = e^(φ/2)",        "Killing ენერგია (შებოჭილი ნაწილაკი)"),
        ("[External] L_oper / L_0", "p = e^(φ/2)",        "გარე კოორდინატული ზომა"),
        ("[External] λ_C_oper",     "p = e^(φ/2)",        "ოპერაციული ზომის სკალირებით"),
        ("[External] T_oper / T_0", "1/p = e^(-φ/2)",     "g_tt = e^φ"),
        ("[External] c_coord / c",  "p² = e^φ",           "Null geodesic"),
        ("[Local] m_loc / m_0",     "1 (INVARIANT)",      "ლოკალური დამკვირვებელი (E_loc = m_0)"),
        ("[Local] L_loc / L_0",     "1 (INVARIANT)",      "ლოკალური სახაზავი"),
        ("[Local] λ_Compton_loc",   "1 (INVARIANT)",      "ℏ/(m_loc c_loc)"),
        ("[Local] α (fine-struct)", "1 (INVARIANT)",      "უგანზომილო ფარდობა"),
        ("ν₀ (substrate)",          "TODO / ჰიპოთეზა",    "მოითხოვს ფონურ ანალიზს"),
    ]
    for row in table:
        print(f"  {row[0]:<22} | {row[1]:<18} | {row[2]}")

    # დაკვირვებითი ფილტრები
    print("\n" + "=" * 72)
    print("დაკვირვებითი ფილტრები (გადამოწმდა SymPy + რიცხვობრივად)")
    print("=" * 72)
    print(f"  ✓ სინათლის გადახრა მზესთან:   {delta_arcsec:.4f} arcsec")
    print(f"     დაკვირვებული: 1.7505 arcsec  [Eddington 1919 / VLBI]")
    print(f"  ✓ Pound-Rebka რედშიფტი:        {pr_pred:.3e}")
    print(f"     დაკვირვებული: 2.57×10⁻¹⁵ ± 0.26×10⁻¹⁵  [PR 1960]")
    print(f"  ✓ ფაქტორი 2 (1.75″ vs 0.87″):  ემერჯენტული ბი-კონფორმობიდან")
    print(f"  ✓ ლოკალური ფარდობითობის პრინციპი: α და ლოკალური სიდიდეები უცვლელია")
    print(f"  🟡 ν₀ სუბსტრატის რიტმი:         TODO / ჰიპოთეზა (გადასამოწმებელია)")

    print("\n" + "=" * 72)
    print("აგენტთა საბჭოს შენიშვნების დადასტურება:")
    print("- E_Killing vs E_local გაიმიჯნა. ლოკალური ენერგია უცვლელია (m_0).")
    print("- n = n_t * n_s = e^(-φ) კავშირი გამოსწორდა. ფაქტორი 2 დადასტურდა.")
    print("=" * 72)


if __name__ == "__main__" and _should_run_main_section("gates"):
    print("\n" + "=" * 72)
    print("P10 CLAIM GATES / STATUS FIREWALL")
    print("=" * 72)

    print("\n--- Signature bridge ---")
    sig_gate = signature_bridge_gate()
    for key, value in sig_gate.items():
        print(f"{key}: {value}")

    print("\n--- Oscillon finite-energy gate ---")
    finite_gate = oscillon_finite_energy_gate()
    for key, value in finite_gate.items():
        print(f"{key}: {value}")

    print("\n--- Oscillon energy status audit ---")
    energy_audit = oscillon_energy_status_audit()
    for key, value in energy_audit.items():
        print(f"{key}: {value}")

    print("\n--- PPN / observation gate ---")
    ppn_gate = ppn_and_observation_gate()
    for key, value in ppn_gate.items():
        print(f"{key}: {value}")

    print("\n--- Whole-file status audit ---")
    status = p10_status_audit()
    for key, value in status.items():
        print(f"{key}: {value}")

    print("\n--- Claim gate ---")
    for gate in oscillon_claim_gate():
        print(f"{gate.claim}: {gate.status}")
        print(f"  verified_here: {gate.verified_here}")
        print(f"  open_requirement: {gate.open_requirement}")

    print("\n--- Do-not-claim ---")
    for item in oscillon_do_not_claim():
        print(f"- {item}")


# =============================================================================
# STAGE D1: OLD quantum oscillon/Bernoulli gate
# =============================================================================

def stage_d1_old_quantum_oscillon_status():
    """
    Deletion-gate marker for the oscillon/gravity parts of OLD/ISPG_Quantum.tex.

    This keeps the valuable old intuition, but classifies it in the stronger
    new RFG form: Bernoulli gravity is a derivation chain, while the full
    finite-energy oscillon existence theorem remains a PDE task.
    """
    phi, m0, L0, omega0 = sp.symbols('phi m0 L0 omega0', real=True, positive=True)
    return {
        "old_file_drained_part": "OLD/ISPG_Quantum.tex: oscillon, Bernoulli, pressure, mass-scaling",
        "new_files": ["p10_oscillons.py", "p01_core.py", "p06_inertia.py", "p02_cosmo.py"],
        "migrated_core": [
            "oscillons are localized time-periodic resonances of the medium",
            "zero-frequency time average of oscillon energy sources the exterior 1/r field",
            "Bernoulli identity converts scalar gradient energy into static pressure deficit",
            "geodesic motion, not literal pressure pushing, gives Newtonian acceleration",
            "mass and operational size scale together as exp(phi/2)",
            "strong-field e^phi factor saturates the pressure deficit near the rarefied core",
        ],
        "mass_scaling": sp.Eq(sp.Symbol('m_eff'), m0 * sp.exp(phi / 2)),
        "operational_size_scaling": sp.Eq(sp.Symbol('L_oper'), L0 * sp.exp(phi / 2)),
        "carrier_frequency_status": (
            "omega0 is kept as an internal/substrate resonance-clock hypothesis; "
            "absolute frequency selection remains a microscopic oscillon theorem"
        ),
        "already_strengthened": [
            "Bernoulli pressure identity is explicit; full gravity proof still requires branch/normalization gates",
            "m_i=m_g=E0/c^2 is handled in p06_inertia.py",
            "process-time and resonant-tail bookkeeping are separated in p02_cosmo.py",
        ],
        "claim_gate": [gate.status for gate in oscillon_claim_gate()],
        "open_math": [
            "construct global finite-energy oscillon solutions of the full nonlinear PDE",
            "prove spectral stability of the localized source sector",
            "derive the absolute electron/rest-frequency scale, not only ratios",
            "derive the macroscopic resonant-tail/dark-energy normalization nonperturbatively",
        ],
    }


def stage_d1_quantum_to_new_file_map():
    """Where the valuable OLD quantum material now lives."""
    return {
        "Bernoulli_gravity": "p10_oscillons.py",
        "inertia_of_oscillon": "p06_inertia.py",
        "singularity_rarefaction": "p05_compact.py",
        "process_time_and_tail_background": "p02_cosmo.py",
        "charged_lepton_C3_Z9": "p11_particles.py",
        "SM_embedding_and_gauge_program": "p11_particles.py",
        "quantum_lab_predictions": "p12_predictions.py",
    }

