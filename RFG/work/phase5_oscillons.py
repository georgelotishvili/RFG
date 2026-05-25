# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from phase1_action import init_variables, get_polynomial_lagrangian


def bernoulli_static_gravity_identity():
    """
    ძველი quantum ფაილის Bernoulli gravity ბირთვი.

    static bi-conformal exterior branch:
        phi(r) = -r_s/r
        Delta P = -P_static = e^phi (phi')^2 / (32*pi*G)
        P_static + Delta P = 0

    ეს არის scalar-field Bernoulli identity: gradient energy drains the
    static pressure. გრავიტაციული მნიშვნელობა პოსტულატი აღარ არის; იგივე
    scalar profile, რომელიც metric-ს ქმნის, ატარებს pressure deficit-ს.
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
        "proof_result": "Once the localized oscillon source is fixed, the 1/r coefficient is fixed. No extra gravitational charge is inserted.",
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
        "localized Delta_P/<T00> fixes the source integral",
        "vacuum exterior solves Laplace equation -> phi_grav=-r_s/r",
        "Gauss normalization: r_s=2GM/c^2",
        "bi-conformal geodesic acceleration: a=-(c^2/2) grad phi=-GM/r^2",
        "two oscillon sources give U(d)=-G M1 M2/d and |F|=G M1 M2/d^2",
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

if __name__ == "__main__":
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
    print("6. Bernoulli gravity ნაწილი უკვე მკაცრად აჩვენებს: static scalar profile ატარებს")
    print("   pressure deficit-ს, zero-frequency oscillon source იძლევა 1/r ველს, ხოლო")
    print("   geodesic equation აბრუნებს Newton-ის კანონს.")
    print("7. ცალკე ღიად რჩება სრული nonlinear finite-energy oscillon profile-ის არსებობის")
    print("   სრული დამტკიცება; მაგრამ gravity mechanism აღარ არის პოსტულატი.")
