# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 18: RFG compact object, Bernoulli saturation, and singularity audit

ეს ფაილი არის compact-object სამუშაო ledger, არა theory-export proof.
აქ მოწმდება exponential exterior-ის algebraic/phenomenological branch:
1. exponential bi-conformal exterior:
       ds^2 = -exp(-r_s/r)c^2dt^2 + exp(r_s/r)(dr^2+r^2dOmega^2)
2. curvature invariants vanish at r -> 0;
3. Bernoulli pressure deficit is self-limiting;
4. finite-radius Killing horizon is absent;
5. areal throat, photon sphere, shadow radius, and golden-ratio ISCO follow;
6. exterior-only geodesic incompleteness is kept explicit;
7. a rarefaction cutoff plus C2 finite-core matching is a conditional
   regular-extension ansatz.

Open before theory export:
    close the p01 polynomial stress/source equations, define physical energy
    and junction stress, compute coupled perturbations/QNMs/echoes, and add
    rotating EHT ray tracing.
"""
import sympy as sp
from p01_core import get_polynomial_lagrangian


def compact_signature_bridge():
    """
    Sign-convention firewall.

    The global RFG work files use signature (+---).  The compact-object block
    below often writes a positive lapse/spatial-factor line element as
        ds^2 = -B(r)c^2dt^2 + A(r)(dr^2+r^2dOmega^2).

    Therefore B and A are positive metric functions, not direct (+---)
    components.  Any stress/Y calculation that touches p01 conventions must
    explicitly pass through this bridge.
    """
    r, r_s = sp.symbols('r r_s', positive=True)
    phi = -r_s / r
    B_lapse = sp.exp(phi)
    A_spatial = sp.exp(-phi)
    return {
        "global_RFG_signature": "(+---)",
        "compact_line_element_style": "ds^2=-B(r)c^2dt^2 + A(r)(dr^2+r^2dOmega^2)",
        "positive_lapse_B": sp.Eq(sp.Symbol('B'), B_lapse),
        "positive_spatial_A": sp.Eq(sp.Symbol('A'), A_spatial),
        "direct_component_warning": "B and A are positive functions; do not treat them as raw (+---) components without conversion.",
        "stress_bridge_status": "REQUIRED_BEFORE_EXPORTING_CENTER_STRESS_CLAIMS",
    }


def derive_exponential_exterior_from_phase_equation():
    """
    Exterior derivation at the RFG phase-potential level.

    Outside the compact source the static pressure/phase potential is harmonic:
        (r^2 phi')' = 0.

    Asymptotic flatness fixes the additive constant; the weak-field Newtonian
    limit g_tt ~= 1 + phi ~= 1 - r_s/r fixes the remaining integration
    constant.  The bi-conformal operational map then gives
        B=e^phi, A=e^-phi.

    This derives the static exponential exterior branch.  It does not yet prove
    that the p01 polynomial medium stress alone supplies the required source.
    """
    r, r_s, C1, C2 = sp.symbols('r r_s C1 C2', positive=True, real=True)
    phi_fn = sp.Function('phi')
    phi_general = C1 + C2 / r
    radial_laplace_residual = sp.simplify(sp.diff(r**2 * sp.diff(phi_general, r), r))
    phi_exterior = -r_s / r
    B_lapse = sp.exp(phi_exterior)
    A_spatial = sp.exp(-phi_exterior)

    weak_gtt_series = sp.series(B_lapse, r_s, 0, 2).removeO()
    weak_spatial_series = sp.series(A_spatial, r_s, 0, 2).removeO()

    return {
        "vacuum_phase_equation": sp.Eq(sp.diff(r**2 * sp.diff(phi_fn(r), r), r), 0),
        "general_spherical_solution": sp.Eq(sp.Symbol('phi'), phi_general),
        "laplace_residual_for_solution": radial_laplace_residual,
        "boundary_conditions": "phi(infinity)=0 and weak-field g_tt ~= 1-r_s/r",
        "exterior_phi": sp.Eq(sp.Symbol('phi_ext'), phi_exterior),
        "biconformal_lapse_B": sp.Eq(sp.Symbol('B'), B_lapse),
        "biconformal_spatial_A": sp.Eq(sp.Symbol('A'), A_spatial),
        "weak_lapse": sp.Eq(sp.Symbol('B_weak'), weak_gtt_series),
        "weak_spatial": sp.Eq(sp.Symbol('A_weak'), weak_spatial_series),
        "derivation_status": "DERIVED_FROM_VACUUM_PHASE_EQUATION_PLUS_BICONFORMAL_MAP",
        "remaining_gate": "p01 polynomial stress/source closure and compact-source matching remain open",
    }


def derive_exponential_effective_source_profile():
    """
    Effective Einstein-source profile of the exponential exterior.

    In the (+---) convention with ds^2=B dt^2 - A(dr^2+r^2dOmega^2),
    B=exp(-r_s/r), A=exp(+r_s/r), the mixed Einstein tensor is proportional
    to the same self-limiting profile as the Bernoulli deficit.
    """
    r, r_s, G = sp.symbols('r r_s G', positive=True, real=True)
    D = sp.simplify(r_s**2 * sp.exp(-r_s / r) / (4 * r**4))
    delta_p = sp.simplify(r_s**2 * sp.exp(-r_s / r) / (32 * sp.pi * G * r**4))
    G_mixed = {
        "G^t_t": -D,
        "G^r_r": D,
        "G^theta_theta": -D,
        "G^phi_phi": -D,
    }
    T_eff = {key.replace("G", "T_eff"): sp.simplify(value / (8 * sp.pi * G)) for key, value in G_mixed.items()}
    return {
        "G_mixed_plus_minus_minus": G_mixed,
        "Bernoulli_Delta_P": sp.Eq(sp.Symbol('Delta_P'), delta_p),
        "T_eff_if_G_eq_8piG_T": T_eff,
        "profile_match": sp.Eq(sp.Symbol('G^r_r/(8*pi*G)'), delta_p),
        "sign_note": "standard Einstein-sign reading gives T^t_t=-Delta_P and T^r_r=+Delta_P; physical energy interpretation is gated by the RFG sign/source convention.",
        "source_status": "EFFECTIVE_GEOMETRIC_SOURCE_PROFILE_DERIVED__PHYSICAL_MEDIUM_SOURCE_OPEN",
    }


def derive_black_hole_singularity_breaker_gate():
    """
    Geometry-level black-hole breaker gate.

    Compare the Schwarzschild finite-radius horizon / r=0 curvature blow-up
    with the RFG exponential exterior branch.  This is the strongest statement
    currently available before p01 source closure:

        Schwarzschild: B=1-r_s/r has B=0 at r=r_s and K~r^-6.
        RFG exponential: B=exp(-r_s/r)>0 for finite r and K->0 at r->0.

    The result excludes the Schwarzschild-type curvature singularity inside
    the exponential phase metric.  It does not by itself prove collapse
    completion, because exterior radial geodesics still reach the r=0 boundary
    in finite affine/proper parameter.
    """
    r, r_s, E_geo = sp.symbols('r r_s E_geo', positive=True, real=True)
    u = r_s / r
    B_exp = sp.exp(-u)
    A_exp = sp.exp(u)
    B_schw = 1 - r_s / r
    K_schw = 12 * r_s**2 / r**6
    K_exp = (
        r_s**2
        * (48 * r**2 - 32 * r * r_s + 7 * r_s**2)
        * sp.exp(-2 * r_s / r)
        / (4 * r**8)
    )
    delta_p = sp.simplify(r_s**2 * sp.exp(-r_s / r) / (32 * sp.pi * sp.Symbol('G', positive=True) * r**4))
    null_radial_velocity_affine = sp.Eq(sp.Symbol('dr_dlambda_null')**2, E_geo**2)
    timelike_radial_velocity_near_boundary = sp.limit(E_geo**2 - B_exp, r, 0, dir='+')

    return {
        "Schwarzschild_horizon_condition": sp.Eq(B_schw, 0),
        "Schwarzschild_horizon_radius": sp.Eq(sp.Symbol('r_h'), r_s),
        "Schwarzschild_Kretschmann": K_schw,
        "lim_r_to_0_K_Schwarzschild": sp.limit(K_schw, r, 0, dir='+'),
        "RFG_lapse": sp.Eq(sp.Symbol('B_exp'), B_exp),
        "RFG_finite_r_horizon_test": "exp(-r_s/r) is strictly positive for every finite r>0",
        "RFG_Kretschmann": K_exp,
        "lim_r_to_0_K_RFG": sp.limit(K_exp, r, 0, dir='+'),
        "RFG_Bernoulli_profile": sp.Eq(sp.Symbol('Delta_P'), delta_p),
        "lim_r_to_0_DeltaP_RFG": sp.limit(delta_p, r, 0, dir='+'),
        "spatial_proper_distance_to_r0": sp.Eq(
            sp.Integral(sp.sqrt(A_exp), (r, 0, sp.Symbol('r0', positive=True))),
            sp.oo,
        ),
        "radial_null_affine_equation": null_radial_velocity_affine,
        "timelike_radial_velocity_limit": sp.Eq(
            sp.Symbol('lim_r_to_0_dr_dtau_squared'),
            timelike_radial_velocity_near_boundary,
        ),
        "geometry_verdict": "SCHWARZSCHILD_CURVATURE_SINGULARITY_REMOVED_IN_EXPONENTIAL_BRANCH",
        "remaining_gate": "not a full black-hole replacement until p01 source closure, boundary/core completion, and perturbative stability are derived",
    }


def p01_polynomial_static_closure_gate():
    """
    Check whether the simplest p01 spherical medium f(r)=r supplies the
    exponential exterior source by itself.

    Result: the f=r branch is isotropic in the spatial mixed stresses, while
    the exponential geometry needs radial/tangential anisotropy.  Therefore the
    exterior is not closed by the simplest p01 polynomial background alone.
    A nontrivial radial deformation f(r), a derivative/Bernoulli sector, or an
    additional source derivation is still required.
    """
    r, r_s = sp.symbols('r r_s', positive=True, real=True)
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        'c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1',
        real=True,
    )
    u = r_s / r
    Tt = -(
        3*c_I1*sp.exp(2*u)
        + 9*c_I1sq*sp.exp(u)
        + 3*c_I2*sp.exp(u)
        + c_I3
        - c_Y*sp.exp(4*u)
        - 3*c_Y2*sp.exp(5*u)
        - 3*c_YI1*sp.exp(3*u)
    ) * sp.exp(-3*u)
    Tr = -(
        c_I1*sp.exp(2*u)
        - 3*c_I1sq*sp.exp(u)
        - c_I2*sp.exp(u)
        - c_I3
        + c_Y*sp.exp(4*u)
        + c_Y2*sp.exp(5*u)
        + c_YI1*sp.exp(3*u)
    ) * sp.exp(-3*u)
    Ttheta = Tr

    D = sp.simplify(r_s**2 * sp.exp(-r_s / r) / (4 * r**4))
    required_anisotropy = sp.simplify(D - (-D))

    return {
        "assumed_elastic_map": "SO(3) branch f(r)=r, so lambda_r=lambda_t=exp(-r_s/r)",
        "p01_T^t_t": sp.factor(Tt),
        "p01_T^r_r": sp.factor(Tr),
        "p01_T^theta_theta": sp.factor(Ttheta),
        "p01_spatial_anisotropy": sp.simplify(Tr - Ttheta),
        "required_geometry_G^r_r_minus_G^theta_theta": required_anisotropy,
        "closure_status": "FAILS_FOR_F_EQ_R__NEEDS_NONTRIVIAL_F_OR_DERIVATIVE_BERNOULLI_SOURCE",
    }


def derive_p01_anisotropic_deformation_route():
    """
    Exact p01 lever that can supply compact-object anisotropy.

    The f=r branch fails because lambda_r=lambda_t.  For a genuine SO(3)
    radial deformation phi^A=f(r)n^A, the p01 polynomial produces a spatial
    anisotropy proportional to lambda_r-lambda_t.  This is the missing degree
    of freedom needed by the exponential geometry.

    This function derives the algebraic anisotropy route; it is not yet the
    solved f(r) field equation.
    """
    r, r_s, A, C, Y, G = sp.symbols('r r_s A C Y G', positive=True, real=True)
    f = sp.Function('f')(r)
    lambda_r = sp.Symbol('lambda_r', positive=True)
    lambda_t = sp.Symbol('lambda_t', positive=True)
    c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        'c_I1 c_I1sq c_I2 c_I3 c_YI1',
        real=True,
    )

    I1 = lambda_r + 2 * lambda_t
    I2 = 2 * lambda_r * lambda_t + lambda_t**2
    I3 = lambda_r * lambda_t**2
    L_solid = (
        c_I1 * I1
        + c_I1sq * I1**2
        + c_I2 * I2
        + c_I3 * I3
        + c_YI1 * Y * I1
    )
    L_lr = sp.diff(L_solid, lambda_r)
    L_lt = sp.diff(L_solid, lambda_t) / 2
    pressure_anisotropy = sp.factor(2 * (lambda_r * L_lr - lambda_t * L_lt))

    lambda_r_map = sp.diff(f, r)**2 / A
    lambda_t_map = f**2 / C
    u = r_s / r
    required_mixed_geometry_anisotropy = sp.simplify(
        r_s**2 * sp.exp(-u) / (2 * r**4)
    )
    required_stress_anisotropy = sp.simplify(
        required_mixed_geometry_anisotropy / (8 * sp.pi * G)
    )
    anisotropy_quadratic = sp.Poly(
        sp.expand(pressure_anisotropy - required_stress_anisotropy),
        lambda_r,
    )
    lambda_r_roots = sp.solve(
        sp.Eq(pressure_anisotropy, required_stress_anisotropy),
        lambda_r,
    )
    deformation_lhs = pressure_anisotropy.subs({
        lambda_r: lambda_r_map,
        lambda_t: lambda_t_map,
    })
    deformation_equation = sp.Eq(deformation_lhs, required_stress_anisotropy)
    A_ext = sp.exp(u)
    C_ext = r**2 * sp.exp(u)
    Y_ext = sp.exp(u)
    exterior_deformation_equation = sp.Eq(
        sp.factor(sp.simplify(deformation_lhs.subs({A: A_ext, C: C_ext, Y: Y_ext}))),
        required_stress_anisotropy,
    )

    return {
        "lambda_r_map": sp.Eq(sp.Symbol('lambda_r'), lambda_r_map),
        "lambda_t_map": sp.Eq(sp.Symbol('lambda_t'), lambda_t_map),
        "solid_anisotropy_p_t_minus_p_r": pressure_anisotropy,
        "f_eq_r_failure_reason": "f=r with C=A*r^2 gives lambda_r=lambda_t, so this anisotropy vanishes",
        "required_geometry_G^r_r_minus_G^theta_theta": required_mixed_geometry_anisotropy,
        "required_stress_anisotropy": required_stress_anisotropy,
        "lambda_r_quadratic_coefficients": anisotropy_quadratic.all_coeffs(),
        "lambda_r_algebraic_roots": lambda_r_roots,
        "deformation_equation_to_solve": deformation_equation,
        "exponential_exterior_f_equation": exterior_deformation_equation,
        "route_status": "EXACT_ANISOTROPY_LEVER_AND_LAMBDA_R_BRANCH_DERIVED__SOLVE_F_R_AND_FULL_STRESS_NEXT",
    }


def derive_minimal_nontrivial_f_branch():
    """
    Minimal nontrivial radial-deformation branch.

    To make the p01 anisotropy concrete, freeze the anisotropy modulus to the
    simplest positive branch:

        c_I1sq = c_I2 = c_YI1 = 0,   c_I1 = K_A > 0.

    Then p_t-p_r = 2 K_A (lambda_r-lambda_t).  On the exponential exterior,
    lambda_r=e^-u f'^2 and lambda_t=e^-u f^2/r^2, so the required compact
    anisotropy becomes a first-order ODE for f(r).

    The closed form below is the asymptotic weak-deformation branch.  It proves
    that f=r is not forced: a nontrivial radial medium map can supply the
    required radial/tangential stress split at leading order.
    """
    r, r_s, G, K_A, eps = sp.symbols(
        'r r_s G K_A eps',
        positive=True,
        real=True,
    )
    chi = sp.Function('chi')(r)
    u = r_s / r
    f = sp.Function('f')(r)
    w = sp.Function('w')(r)

    lambda_r = sp.exp(-u) * sp.diff(f, r)**2
    lambda_t = sp.exp(-u) * f**2 / r**2
    required = r_s**2 * sp.exp(-u) / (16 * sp.pi * G * r**4)
    minimal_exact_ode = sp.Eq(
        sp.simplify(2 * K_A * (lambda_r - lambda_t)),
        required,
    )
    simplified_exact_ode = sp.Eq(
        sp.diff(f, r)**2 - f**2 / r**2,
        r_s**2 / (32 * sp.pi * G * K_A * r**4),
    )
    w_equation = sp.Eq(
        r * sp.diff(w, r),
        -w + sp.sqrt(w**2 + r_s**2 / (32 * sp.pi * G * K_A * r**4)),
    )

    f_eps = r * (1 + eps * chi)
    lambda_r_eps = sp.exp(-u) * sp.diff(f_eps, r)**2
    lambda_t_eps = sp.exp(-u) * f_eps**2 / r**2
    anisotropy_eps = sp.series(
        sp.simplify(2 * K_A * (lambda_r_eps - lambda_t_eps)),
        eps,
        0,
        2,
    ).removeO()
    linear_anisotropy = sp.simplify(sp.diff(anisotropy_eps, eps).subs(eps, 0))
    chi_prime_equation = sp.Eq(
        sp.diff(chi, r),
        r_s**2 / (64 * sp.pi * G * K_A * r**5),
    )
    chi_solution = -r_s**2 / (256 * sp.pi * G * K_A * r**4)
    f_asymptotic = sp.simplify(r * (1 + chi_solution))
    f_asymptotic_prime = sp.diff(f_asymptotic, r)
    asymptotic_lhs_first_order = sp.simplify(
        2
        * K_A
        * sp.exp(-u)
        * (
            2
            * r
            * sp.diff(chi_solution, r)
        )
    )
    linear_residual = sp.simplify(asymptotic_lhs_first_order - required)
    exact_residual_order = sp.simplify(
        2 * K_A * sp.exp(-u) * (f_asymptotic_prime**2 - f_asymptotic**2 / r**2)
        - required
    )

    return {
        "minimal_modulus_branch": "c_I1sq=c_I2=c_YI1=0, K_A=c_I1>0",
        "exact_minimal_f_ode": minimal_exact_ode,
        "simplified_exact_ode": simplified_exact_ode,
        "asymptotically_flat_w_equation_for_f_eq_r_w": w_equation,
        "linearized_f_definition": sp.Eq(sp.Symbol('f'), r * (1 + chi)),
        "linear_anisotropy": sp.Eq(sp.Symbol('Delta_p_linear'), linear_anisotropy),
        "chi_prime_equation": chi_prime_equation,
        "chi_solution_with_chi_infinity_0": sp.Eq(sp.Symbol('chi'), chi_solution),
        "nontrivial_f_asymptotic": sp.Eq(sp.Symbol('f_nontrivial'), f_asymptotic),
        "linear_source_residual": linear_residual,
        "exact_residual_after_linear_solution": exact_residual_order,
        "validity": "asymptotic/weak-deformation exterior branch; exact nonlinear f(r), full p01 stress, and core matching remain open",
        "branch_status": "NONTRIVIAL_F_R_BRANCH_DERIVED_AT_LINEAR_ORDER__EXACT_ODE_READY",
    }


def analyze_exponential_exterior_curvature():
    """
    Exponential compact-object exterior branch.

    metric:
        ds^2 = -exp(-r_s/r)c^2dt^2 + exp(r_s/r)(dr^2+r^2dOmega^2)

    Branch-level invariant check:
        R -> 0 and K -> 0 as r -> 0.
    """
    r, r_s = sp.symbols('r r_s', real=True, positive=True)
    u = sp.Symbol('u', real=True, positive=True)

    phi = -r_s / r
    g_tt = -sp.exp(phi)
    g_rr = sp.exp(-phi)

    ricci = -r_s**2 * sp.exp(-r_s / r) / (2 * r**4)
    ricci_sq = r_s**4 * sp.exp(-2 * r_s / r) / (4 * r**8)
    kretschmann = (
        r_s**2
        * (48 * r**2 - 32 * r * r_s + 7 * r_s**2)
        * sp.exp(-2 * r_s / r)
        / (4 * r**8)
    )

    k_shape = sp.simplify(u**6 * (48 - 32 * u + 7 * u**2) * sp.exp(-2 * u) / 4)
    k_extremum_poly = 7 * u**3 - 60 * u**2 + 160 * u - 144
    physical_roots = [
        root for root in sp.nroots(k_extremum_poly)
        if abs(sp.im(root)) < 1.0e-12 and sp.re(root) > 0
    ]
    u_peak = sp.N(sp.re(physical_roots[0]), 8) if physical_roots else sp.nan
    r_peak = sp.N(1 / u_peak, 8) * r_s if physical_roots else sp.nan
    k_peak_coeff = sp.N(k_shape.subs(u, u_peak), 8) if physical_roots else sp.nan

    return {
        "ansatz_result": "within the exponential exterior, curvature scalars do not blow up at r->0",
        "derivation_status": "DERIVED_AT_PHASE_METRIC_LEVEL__FULL_P01_SOURCE_CLOSURE_OPEN",
        "phi": sp.Eq(sp.Symbol('phi'), phi),
        "g_tt": sp.Eq(sp.Symbol('g_tt'), g_tt),
        "g_rr": sp.Eq(sp.Symbol('g_rr'), g_rr),
        "Ricci_scalar": ricci,
        "Ricci_squared": ricci_sq,
        "Kretschmann": kretschmann,
        "lim_r_to_0_R": sp.limit(ricci, r, 0, dir='+'),
        "lim_r_to_0_Ricci2": sp.limit(ricci_sq, r, 0, dir='+'),
        "lim_r_to_0_K": sp.limit(kretschmann, r, 0, dir='+'),
        "K_extremum_polynomial_u_rs_over_r": sp.Eq(k_extremum_poly, 0),
        "K_peak_u": u_peak,
        "K_peak_r": r_peak,
        "K_peak": sp.Eq(sp.Symbol('K_max'), k_peak_coeff / r_s**4),
    }


def analyze_bernoulli_singularity_saturation():
    """
    Bernoulli pressure-deficit branch check.

    Delta P = exp(phi)(phi')^2/(32*pi*G)
    for phi=-r_s/r gives
        Delta P = r_s^2 exp(-r_s/r)/(32*pi*G*r^4).

    It peaks at r_s/4 and vanishes at r->0.  This is not by itself an
    invariant total-energy proof.
    """
    r, r_s, G = sp.symbols('r r_s G', real=True, positive=True)
    phi = -r_s / r
    delta_p = sp.simplify(sp.exp(phi) * sp.diff(phi, r)**2 / (32 * sp.pi * G))
    p_static = -delta_p
    coordinate_energy = sp.simplify(
        sp.integrate(delta_p * 4 * sp.pi * r**2, (r, 0, sp.oo))
    )
    proper_volume_integrand = 4 * sp.pi * r**2 * sp.exp(3 * r_s / (2 * r))
    naive_proper_energy_integrand = sp.simplify(delta_p * proper_volume_integrand)
    u = sp.Symbol('u', real=True, positive=True)
    shape = sp.exp(-u) * u**4

    return {
        "ansatz_result": "Bernoulli pressure deficit has a finite peak and returns to zero at r->0",
        "Delta_P": sp.Eq(sp.Symbol('Delta_P'), delta_p),
        "P_static": sp.Eq(sp.Symbol('P_static'), p_static),
        "Bernoulli_identity_definition": sp.Eq(sp.Symbol('P_static + Delta_P'), 0),
        "dimensionless_shape": sp.Eq(sp.Symbol('shape(u)'), shape),
        "shape_derivative": sp.factor(sp.diff(shape, u)),
        "pressure_peak": sp.Eq(sp.Symbol('r_peak'), r_s / 4),
        "lim_r_to_0_Delta_P": sp.limit(delta_p, r, 0, dir='+'),
        "lim_r_to_inf_Delta_P": sp.limit(delta_p, r, sp.oo),
        "coordinate_measure_energy": sp.Eq(sp.Symbol('int_DeltaP_4pi_r2_dr'), coordinate_energy),
        "naive_proper_energy_integrand": naive_proper_energy_integrand,
        "lim_r_to_0_naive_proper_integrand": sp.limit(naive_proper_energy_integrand, r, 0, dir='+'),
        "energy_gate": "coordinate measure is finite, but invariant/proper/ADM energy is not established here",
        "meaning": "as phi -> -infinity, exp(phi) shuts off this branch's local gradient-energy density; export requires a physical energy measure.",
    }


def analyze_rarefaction_information_cutoff():
    """
    Phenomenological microscopic closure for the r=0 boundary.

    User intuition, written as a continuum-mechanics criterion:
    in the deepest Bernoulli deficit the active carrier density of the
    vacuum medium is rarefied, collisions become sparse, and the continuum
    no longer transmits information as a connected elastic fluid.

    Minimal closure ansatz:
        n_eff = n_0 exp(phi),       phi=-r_s/r
        c_eff = c exp(phi)
        ell_mfp = 1/(sigma n_eff)
        Gamma_coll = c_eff/ell_mfp.

    Then Gamma_coll -> 0, ell_mfp -> infinity, Kn=ell_mfp/r -> infinity.
    This converts the formal r=0 endpoint into an information-decoupled,
    dilute boundary/core instead of an infinite-density singularity.
    """
    r, r_s, G, n_0, sigma, c, a_osc = sp.symbols(
        'r r_s G n_0 sigma c a_osc',
        real=True,
        positive=True,
    )
    phi = -r_s / r
    n_eff = n_0 * sp.exp(phi)
    mean_spacing = n_eff ** (-sp.Rational(1, 3))
    ell_mfp = 1 / (sigma * n_eff)
    c_eff = c * sp.exp(phi)
    gamma_coll = sp.simplify(c_eff / ell_mfp)
    gamma_ratio = sp.simplify(gamma_coll / (c * sigma * n_0))

    gradient_length = sp.simplify(abs(phi / sp.diff(phi, r)))
    knudsen = sp.simplify(ell_mfp / r)
    carriers_in_gradient_cell = sp.simplify(n_eff * r**3)
    carriers_in_finite_oscillon = sp.simplify(n_eff * a_osc**3)

    delta_p = sp.simplify(sp.exp(phi) * sp.diff(phi, r)**2 / (32 * sp.pi * G))
    outward_pressure_force_density = sp.factor(sp.diff(delta_p, r))

    return {
        "closure_status": "PHENOMENOLOGICAL_KNUDSEN_CUTOFF_NOT_YET_DERIVED_FROM_ACTION",
        "closure_density": sp.Eq(sp.Symbol('n_eff'), n_eff),
        "mean_spacing": sp.Eq(sp.Symbol('d_eff'), mean_spacing),
        "mean_free_path": sp.Eq(sp.Symbol('ell_mfp'), ell_mfp),
        "effective_signal_speed": sp.Eq(sp.Symbol('c_eff'), c_eff),
        "collision_rate": sp.Eq(sp.Symbol('Gamma_coll'), gamma_coll),
        "collision_rate_ratio": sp.Eq(sp.Symbol('Gamma_coll/Gamma_0'), gamma_ratio),
        "lim_r_to_0_n_eff": sp.limit(n_eff, r, 0, dir='+'),
        "lim_r_to_0_mean_spacing": sp.limit(mean_spacing, r, 0, dir='+'),
        "lim_r_to_0_ell_mfp": sp.limit(ell_mfp, r, 0, dir='+'),
        "lim_r_to_0_collision_ratio": sp.limit(gamma_ratio, r, 0, dir='+'),
        "gradient_length": sp.Eq(sp.Symbol('L_grad'), gradient_length),
        "Knudsen_number": sp.Eq(sp.Symbol('Kn'), knudsen),
        "lim_r_to_0_Kn": sp.limit(knudsen, r, 0, dir='+'),
        "carriers_in_gradient_cell": sp.Eq(sp.Symbol('N_grad'), carriers_in_gradient_cell),
        "lim_r_to_0_N_grad": sp.limit(carriers_in_gradient_cell, r, 0, dir='+'),
        "carriers_in_finite_oscillon": sp.Eq(sp.Symbol('N_osc'), carriers_in_finite_oscillon),
        "lim_r_to_0_N_osc": sp.limit(carriers_in_finite_oscillon, r, 0, dir='+'),
        "medium_stress_gradient": sp.Eq(sp.Symbol('d_DeltaP_dr'), outward_pressure_force_density),
        "force_turning_radius": sp.Eq(sp.Symbol('r_turn'), r_s / 4),
        "inner_core_sign": "for r<r_s/4, d(Delta_P)/dr>0: medium-stress gradient reverses sign.",
        "force_firewall": "medium backreaction only; not a literal extra pressure force on matter.",
        "physical_meaning": "near r=0 the closure ansatz makes the medium dilute and non-communicating; this remains a kinetic derivation task.",
    }


def analyze_geodesic_completion_by_core_matching():
    """
    Conditional boundary-extension ansatz.

    The exponential exterior is not forced to run all the way to r=0 if the
    Knudsen cutoff is physically allowed.
    The continuum description self-terminates where Kn=ell_mfp/r reaches 1:

        exp(r_s/r_c)/(n_0 sigma r_c) = 1
        r_c = r_s / W(n_0 sigma r_s).

    For r<r_c use a regular kinetic/rarefied core.  Write

        ds^2 = -B(r)c^2dt^2 + A(r)(dr^2+r^2dOmega^2).

    Let q=r_s/r_c and x=r/r_c.  A C2 positive logarithmic core that matches
    the exterior A_+=exp(r_s/r), B_+=exp(-r_s/r) through second derivative is:

        log A_- = q(35x^2/8 - 21x^4/4 + 15x^6/8),
        log B_- = -q + q(-11x^2/8 + 9x^4/4 - 7x^6/8).

    This gives A(0)=1, B(0)=exp(-q)>0, first derivatives vanish at the center,
    and all center curvature scalars are finite.  In Cartesian coordinates the
    center is locally regular.  A full completion proof still requires the
    core stress tensor, junction conditions, energy bookkeeping and stability.
    """
    r, r_s, r_c, n_0, sigma, q, x = sp.symbols(
        'r r_s r_c n_0 sigma q x',
        real=True,
        positive=True,
    )

    r_kn = sp.simplify(r_s / sp.LambertW(n_0 * sigma * r_s))
    q_kn = sp.simplify(sp.LambertW(n_0 * sigma * r_s))
    pressure_turn_condition = sp.Ge(n_0 * sigma * r_s, 4 * sp.exp(4))

    log_a_core = q * (
        sp.Rational(35, 8) * x**2
        - sp.Rational(21, 4) * x**4
        + sp.Rational(15, 8) * x**6
    )
    log_b_core = -q + q * (
        -sp.Rational(11, 8) * x**2
        + sp.Rational(9, 4) * x**4
        - sp.Rational(7, 8) * x**6
    )
    log_a_ext = q / x
    log_b_ext = -q / x

    c2_match_a = [
        sp.simplify(sp.diff(log_a_core, x, order).subs(x, 1)
                    - sp.diff(log_a_ext, x, order).subs(x, 1))
        for order in range(3)
    ]
    c2_match_b = [
        sp.simplify(sp.diff(log_b_core, x, order).subs(x, 1)
                    - sp.diff(log_b_ext, x, order).subs(x, 1))
        for order in range(3)
    ]

    a2 = sp.Rational(35, 8) * q / r_c**2
    b2_over_b0 = -sp.Rational(11, 8) * q / r_c**2
    center_ricci = sp.simplify(6 * a2 - 6 * b2_over_b0)
    center_kretschmann = sp.simplify(12 * a2**2 + 12 * b2_over_b0**2)

    return {
        "conditional_ansatz": "C2 finite-core matching gives a locally regular extension if the Knudsen cutoff is allowed",
        "proof_status": "NOT_A_FIELD_EQUATION_DERIVATION__JUNCTION_STRESS_AND_STABILITY_OPEN",
        "Kn_cutoff_equation": sp.Eq(sp.exp(r_s / r_c) / (n_0 * sigma * r_c), 1),
        "core_radius": sp.Eq(sp.Symbol('r_c'), r_kn),
        "core_compactness_q": sp.Eq(sp.Symbol('q_c'), q_kn),
        "inside_pressure_reversal_condition": pressure_turn_condition,
        "x_definition": sp.Eq(sp.Symbol('x'), r / r_c),
        "log_A_core": sp.Eq(sp.Symbol('log_A_minus'), log_a_core),
        "log_B_core": sp.Eq(sp.Symbol('log_B_minus'), log_b_core),
        "log_A_exterior": sp.Eq(sp.Symbol('log_A_plus'), log_a_ext),
        "log_B_exterior": sp.Eq(sp.Symbol('log_B_plus'), log_b_ext),
        "C2_match_log_A_value_slope_curvature": c2_match_a,
        "C2_match_log_B_value_slope_curvature": c2_match_b,
        "center_A": sp.Eq(sp.Symbol('A_0'), 1),
        "center_B": sp.Eq(sp.Symbol('B_0'), sp.exp(-q)),
        "center_A_prime": sp.Eq(sp.Symbol("A'_0"), 0),
        "center_B_prime": sp.Eq(sp.Symbol("B'_0"), 0),
        "center_Ricci_scalar": sp.Eq(sp.Symbol('R_0'), center_ricci),
        "center_Kretschmann": sp.Eq(sp.Symbol('K_0'), center_kretschmann),
        "local_continuation_rule": "bounded center coefficients allow local geodesic continuation inside the ansatz core.",
        "global_completion_gate": "open until field-equation source, junction stress, ADM/proper energy and perturbative stability are checked.",
        "physical_meaning": "the formal endpoint is replaced by a candidate dilute kinetic core fixed by Kn=1.",
    }


def derive_c2_core_matching_coefficients():
    """
    Derive the polynomial coefficients used in the conditional C2 core.

    This proves only the matching algebra.  It is not a field-equation source
    derivation for the core material.
    """
    x, q = sp.symbols('x q', positive=True, real=True)
    a2, a4, a6, b2, b4, b6 = sp.symbols('a2 a4 a6 b2 b4 b6', real=True)

    log_a_poly = q * (a2*x**2 + a4*x**4 + a6*x**6)
    log_a_ext = q / x
    a_solution = sp.solve(
        [
            sp.Eq(log_a_poly.subs(x, 1), log_a_ext.subs(x, 1)),
            sp.Eq(sp.diff(log_a_poly, x).subs(x, 1), sp.diff(log_a_ext, x).subs(x, 1)),
            sp.Eq(sp.diff(log_a_poly, x, 2).subs(x, 1), sp.diff(log_a_ext, x, 2).subs(x, 1)),
        ],
        [a2, a4, a6],
        dict=True,
    )[0]

    log_b_poly = -q + q * (b2*x**2 + b4*x**4 + b6*x**6)
    log_b_ext = -q / x
    b_solution = sp.solve(
        [
            sp.Eq(log_b_poly.subs(x, 1), log_b_ext.subs(x, 1)),
            sp.Eq(sp.diff(log_b_poly, x).subs(x, 1), sp.diff(log_b_ext, x).subs(x, 1)),
            sp.Eq(sp.diff(log_b_poly, x, 2).subs(x, 1), sp.diff(log_b_ext, x, 2).subs(x, 1)),
        ],
        [b2, b4, b6],
        dict=True,
    )[0]

    return {
        "log_A_coefficients": a_solution,
        "log_B_coefficients": b_solution,
        "log_A_core_derived": sp.Eq(sp.Symbol('log_A_minus'), sp.simplify(log_a_poly.subs(a_solution))),
        "log_B_core_derived": sp.Eq(sp.Symbol('log_B_minus'), sp.simplify(log_b_poly.subs(b_solution))),
        "derivation_status": "C2_MATCHING_COEFFICIENTS_DERIVED__CORE_FIELD_EQUATIONS_OPEN",
    }


def analyze_horizon_throat_and_boundary():
    """
    Horizonless exterior-only boundary ledger.

    finite r>0-ზე g_tt never vanishes; r=0 is a boundary, not a finite
    horizon. Proper distance to r=0 diverges, while captured geodesics
    reach the boundary at finite affine/proper parameter before the
    conditional rarefied-core extension is imposed.
    """
    r, r_s, r_0, c, E_geo = sp.symbols('r r_s r_0 c E_geo', real=True, positive=True)
    phi = -r_s / r
    g_tt_abs = sp.exp(phi)
    c_coord = c * sp.exp(phi)
    proper_integrand = sp.exp(r_s / (2 * r))
    coordinate_null_integrand = sp.exp(r_s / r) / c

    areal_radius = r * sp.exp(r_s / (2 * r))
    d_areal = sp.diff(areal_radius, r)
    throat = r_s / 2
    areal_min = sp.simplify(areal_radius.subs(r, throat))

    return {
        "exterior_result": "no finite-radius Killing horizon in the static exponential exterior",
        "g_tt_abs": sp.Eq(sp.Symbol('|g_tt|'), g_tt_abs),
        "finite_r_horizon_test": "exp(-r_s/r)>0 for every finite r>0",
        "coordinate_light_speed": sp.Eq(sp.Symbol('dr_dt_null'), c_coord),
        "lim_r_to_0_c_coord": sp.limit(c_coord, r, 0, dir='+'),
        "proper_distance_integrand": proper_integrand,
        "proper_distance_to_boundary": sp.Eq(sp.Symbol('L_prop'), sp.oo),
        "external_coordinate_time_to_boundary": sp.Eq(sp.Symbol('t_external'), sp.oo),
        "radial_null_affine_arrival": sp.Eq(sp.Symbol('lambda_0'), r_0 * c / E_geo),
        "areal_radius": sp.Eq(sp.Symbol('R_areal'), areal_radius),
        "dR_dr": d_areal,
        "throat_coordinate": sp.Eq(sp.Symbol('r_throat'), throat),
        "throat_areal_radius": sp.Eq(sp.Symbol('R_min'), areal_min),
        "boundary_label": "r=0 is an infinite-redshift exterior boundary, not yet a proven collapse endpoint.",
        "geodesic_status": "exterior-only captured geodesics reach r=0 in finite affine parameter; the Knudsen core is a conditional extension ansatz.",
    }


def analyze_photon_shadow_isco():
    """
    Strong-field observables of the exponential exterior.

    Timelike geodesics use
        rdot^2 = E^2 - V_eff,
        V_eff = exp(-r_s/r) + L^2 exp(-2r_s/r)/r^2.

    Circular orbit:
        dV_eff/dr = 0.
    Marginal stability:
        d^2V_eff/dr^2 = 0.

    photon sphere:
        d/dr [exp(-2r_s/r)/r^2]=0 -> r=r_s
    shadow:
        b_c = e r_s
    massive ISCO:
        r^2 - 3 r_s r + r_s^2 = 0 -> r_ISCO = phi_golden^2 r_s
    """
    r, r_s, c, L = sp.symbols('r r_s c L', real=True, positive=True)
    phi_golden = (1 + sp.sqrt(5)) / 2

    photon_barrier = sp.exp(-2 * r_s / r) / r**2
    photon_condition = sp.factor(sp.diff(photon_barrier, r))

    v_eff = sp.exp(-r_s / r) + L**2 * sp.exp(-2 * r_s / r) / r**2
    circular_condition = sp.factor(sp.diff(v_eff, r))
    specific_l2 = sp.simplify(r_s * r**2 * sp.exp(r_s / r) / (2 * (r - r_s)))
    specific_e2 = sp.simplify(
        sp.exp(-r_s / r) + specific_l2 * sp.exp(-2 * r_s / r) / r**2
    )

    v_eff_second = sp.diff(v_eff, r, 2)
    stability_second_derivative = sp.factor(
        sp.simplify(v_eff_second.subs(L**2, specific_l2))
    )
    stability_polynomial = r**2 - 3 * r_s * r + r_s**2
    isco_poly = r**2 - 3 * r_s * r + r_s**2
    isco_roots = sp.solve(sp.Eq(isco_poly, 0), r)
    isco_physical = sp.simplify(isco_roots[1])

    omega_sq = c**2 * r_s * sp.exp(-2 * r_s / r) / (r**2 * (2 * r - r_s))
    omega_isco_sq = sp.simplify(omega_sq.subs(r, isco_physical))
    omega_gr_sq = c**2 / (54 * r_s**2)
    omega_ratio = sp.N(sp.sqrt(sp.simplify(omega_isco_sq / omega_gr_sq)), 8)

    local_speed_sq = sp.simplify(c**2 * r_s / (2 * r - r_s))
    local_speed_isco = sp.simplify(sp.sqrt(local_speed_sq.subs(r, isco_physical)))
    e2_isco = sp.simplify(specific_e2.subs(r, isco_physical))
    e_isco = sp.sqrt(e2_isco)
    binding_efficiency = sp.N(1 - e_isco, 8)
    gr_binding_efficiency = sp.N(1 - sp.sqrt(sp.Rational(8, 9)), 8)

    gr_isco_iso = (5 + 2 * sp.sqrt(6)) * r_s / 4
    gr_isco_areal = 3 * r_s
    isco_iso_ratio = sp.N(isco_physical / gr_isco_iso, 8)
    shadow_rfg = sp.E * r_s
    shadow_gr = 3 * sp.sqrt(3) * r_s / 2
    shadow_ratio = sp.N(shadow_rfg / shadow_gr, 8)

    return {
        "effective_potential_timelike": sp.Eq(sp.Symbol('V_eff'), v_eff),
        "circular_orbit_condition": sp.Eq(circular_condition, 0),
        "specific_L_squared_circular": sp.Eq(sp.Symbol('L_circ^2'), specific_l2),
        "specific_E_squared_circular": sp.Eq(sp.Symbol('E_circ^2'), specific_e2),
        "massive_orbit_existence": "L_circ^2>0 requires r>r_s; r=r_s is the photon boundary.",
        "local_orbital_speed_squared": sp.Eq(sp.Symbol('v_local^2'), local_speed_sq),
        "photon_barrier": photon_barrier,
        "photon_condition": sp.Eq(photon_condition, 0),
        "photon_sphere": sp.Eq(sp.Symbol('r_ph'), r_s),
        "critical_impact_parameter": sp.Eq(sp.Symbol('b_c'), shadow_rfg),
        "GR_shadow_reference": sp.Eq(sp.Symbol('b_c_GR'), shadow_gr),
        "shadow_ratio_RFG_over_GR": shadow_ratio,
        "shadow_size_shift": f"{float((shadow_ratio - 1) * 100):.2f}%",
        "stability_second_derivative": stability_second_derivative,
        "stability_polynomial": sp.Eq(stability_polynomial, 0),
        "ISCO_polynomial": sp.Eq(isco_poly, 0),
        "ISCO_roots": isco_roots,
        "ISCO_physical": sp.Eq(sp.Symbol('r_ISCO'), isco_physical),
        "golden_ratio_identity": sp.Eq(sp.Symbol('r_ISCO'), phi_golden**2 * r_s),
        "ISCO_local_speed": sp.Eq(sp.Symbol('v_ISCO'), local_speed_isco),
        "ISCO_specific_energy_squared": sp.Eq(sp.Symbol('E_ISCO^2'), e2_isco),
        "ISCO_binding_efficiency": binding_efficiency,
        "GR_binding_efficiency_reference": gr_binding_efficiency,
        "GR_areal_ISCO_reference": sp.Eq(sp.Symbol('R_ISCO_GR'), gr_isco_areal),
        "GR_isotropic_ISCO_reference": sp.Eq(sp.Symbol('r_ISCO_GR_iso'), gr_isco_iso),
        "ISCO_radius_ratio_RFG_over_GR_iso": isco_iso_ratio,
        "Omega_ISCO_squared": sp.Eq(sp.Symbol('Omega_ISCO^2'), omega_isco_sq),
        "GR_Omega_ISCO_squared_reference": sp.Eq(sp.Symbol('Omega_GR^2'), omega_gr_sq),
        "frequency_ratio_RFG_over_GR": omega_ratio,
        "frequency_proxy": "f_ISCO = 0.931 f_ISCO_GR for the same total mass",
        "mechanism": "golden-ratio ISCO is the marginal-stability root of the exponential-vacuum geodesic potential, not a fitted number.",
    }


def singularity_strength_ledger() -> list[str]:
    return [
        "Geometry-level breaker: Schwarzschild has B=0 at r_s and K->infinity at r=0; the RFG exponential branch has B>0 for every finite r>0 and K->0 at r=0.",
        "Within the exponential exterior ansatz: R->0, Ricci^2->0, K->0 at r->0.",
        "Within the Bernoulli branch: Delta_P peaks at r_s/4 and returns to 0 at r->0.",
        "The r=0 endpoint is not a curvature singularity in the exponential branch, but exterior radial geodesics still expose a boundary unless a derived core/boundary law is added.",
        "The exact p01 anisotropy lever is lambda_r-lambda_t; f=r kills it, so a nontrivial radial deformation f(r) is the next closure target.",
        "A minimal nontrivial branch gives f=r*(1-r_s^2/(256*pi*G*K_A*r^4)) at linear order and exactly cancels the required anisotropy at that order.",
        "The rarefaction closure ansatz gives n_eff->0, mean free path->infinity, and collision/information rate->0 at r->0.",
        "The Knudsen number Kn=ell_mfp/r diverges, so the continuum model flags its own breakdown.",
        "Inside r_s/4 the medium-stress gradient changes sign; this is backreaction language, not an extra force on matter.",
        "If Kn=1 is the physical cutoff, the candidate matching radius is r_c=r_s/W(n_0*sigma*r_s).",
        "The C2 logarithmic core matches the exponential exterior through value, slope, and curvature at r_c.",
        "The matched core has A(0)=1, B(0)>0, A'(0)=B'(0)=0 and finite R_0, K_0 in the ansatz.",
        "There is no finite-radius Killing horizon in the static exterior: exp(-r_s/r)>0 for every r>0.",
        "The areal radius has a throat at r_s/2, with R_min=e*r_s/2.",
        "Static photon sphere is r_s and the static critical shadow impact parameter is b_c=e*r_s.",
        "Massive circular orbits exist only for r>r_s; r=r_s is a static photon-orbit boundary.",
        "Static massive-particle ISCO follows from V_eff''=0 and is r_ISCO=phi_golden^2*r_s.",
        "The static ISCO frequency proxy is f_ISCO=0.931 f_ISCO_GR for the same total mass.",
        "External observers see infinite redshift/coordinate-time freezing toward r=0 in the exterior.",
        "A full compact-object claim remains blocked by field-equation derivation, physical energy, junction stress, rotation, QNMs/echoes and EHT ray tracing.",
    ]


def analyze_regular_center():
    r = sp.Symbol('r', real=True, positive=True)
    a_2, b_2 = sp.symbols('a_2 b_2', real=True)
    B_0 = sp.Symbol('B_0', real=True, positive=True) # B(0) > 0 აუცილებელია
    
    # რეგულარული ცენტრის ანზაცი (A(0) = 1 აუცილებელია რეგულარულობისთვის)
    A_core = 1 + a_2 * r**2
    B_core = B_0 + b_2 * r**2
    
    # G^mu_nu გეომეტრიული ნაწილები
    G_tt = -sp.diff(A_core, r) / (r * A_core**2) + (1/A_core - 1)/r**2
    G_rr = sp.diff(B_core, r) / (r * A_core * B_core) + (1/A_core - 1)/r**2
    G_thth = sp.diff(B_core, r, 2)/(2*A_core*B_core) - sp.diff(B_core, r)**2/(4*A_core*B_core**2) - sp.diff(A_core, r)*sp.diff(B_core, r)/(4*A_core**2*B_core) + sp.diff(B_core, r)/(2*r*A_core*B_core) - sp.diff(A_core, r)/(2*r*A_core**2)
    
    # რიჩის სკალარი R = -G^\mu_\mu
    R_scalar = -(G_tt + G_rr + 2*G_thth)
    
    # ლიმიტები r -> 0
    G_tt_0 = sp.simplify(sp.limit(G_tt, r, 0))
    G_rr_0 = sp.simplify(sp.limit(G_rr, r, 0))
    G_thth_0 = sp.simplify(sp.limit(G_thth, r, 0))
    R_0 = sp.simplify(sp.limit(R_scalar, r, 0))
    
    # Kretschmann სკალარი (K = R_{abcd}R^{abcd}) რეალური სინგულარობის შესამოწმებლად
    A_p, B_p = sp.diff(A_core, r), sp.diff(B_core, r)
    B_pp = sp.diff(B_p, r)
    term1 = (B_pp/(2*B_core) - B_p**2/(4*B_core**2) - A_p*B_p/(4*A_core*B_core))**2
    term2 = (B_p/(r*B_core))**2
    term3 = (A_p/(r*A_core))**2
    term4 = (1 - 1/A_core)**2 / r**4
    K_scalar = (4 / A_core**2) * term1 + (2 / A_core**2) * term2 + (2 / A_core**2) * term3 + 4 * term4
    K_0 = sp.simplify(sp.limit(K_scalar, r, 0))
    
    # სუპერსოლიდის სტრეს-ტენზორი r -> 0-სას
    # ვუშვებთ ცენტრალური სკალარული მუხტის არარსებობას (Psi'=0)
    Y = 1/B_core
    I1 = 2 + 1/A_core
    I2 = 1 + 2/A_core
    I3 = 1/A_core
    
    Y_s, I1_s, I2_s, I3_s = sp.symbols('Y I1 I2 I3', real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    
    L_eval = L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_Y = sp.diff(L_poly, Y_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I1 = sp.diff(L_poly, I1_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I2 = sp.diff(L_poly, I2_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I3 = sp.diff(L_poly, I3_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    
    T_tt = 2 * L_Y / B_core - L_eval
    T_rr = 2 * (L_I1 / A_core + 2 * L_I2 / A_core + L_I3 / A_core) - L_eval
    
    T_tt_0 = sp.simplify(sp.limit(T_tt, r, 0))
    T_rr_0 = sp.simplify(sp.limit(T_rr, r, 0))
    
    return G_tt_0, G_rr_0, G_thth_0, R_0, K_0, T_tt_0, T_rr_0

if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 18: RFG compact object and singularity audit")
    print("=" * 72)

    print("\n0. Sign-convention bridge")
    for key, value in compact_signature_bridge().items():
        print(f"  {key:36s}: {value}")

    print("\n1. Exterior derivation from vacuum phase equation")
    phase_derivation = derive_exponential_exterior_from_phase_equation()
    for key, value in phase_derivation.items():
        print(f"  {key:36s}: {value}")

    print("\n1b. Effective source profile of the exponential exterior")
    source_profile = derive_exponential_effective_source_profile()
    for key, value in source_profile.items():
        print(f"  {key:36s}: {value}")

    print("\n1c. Black-hole singularity breaker gate")
    bh_breaker = derive_black_hole_singularity_breaker_gate()
    for key, value in bh_breaker.items():
        print(f"  {key:36s}: {value}")

    print("\n1d. p01 polynomial static closure gate")
    p01_gate = p01_polynomial_static_closure_gate()
    for key, value in p01_gate.items():
        print(f"  {key:36s}: {value}")

    print("\n1e. p01 anisotropic deformation route")
    anisotropic_route = derive_p01_anisotropic_deformation_route()
    for key, value in anisotropic_route.items():
        print(f"  {key:36s}: {value}")

    print("\n1f. Minimal nontrivial f(r) branch")
    nontrivial_f = derive_minimal_nontrivial_f_branch()
    for key, value in nontrivial_f.items():
        print(f"  {key:36s}: {value}")

    print("\n1g. Exponential exterior curvature branch")
    exterior = analyze_exponential_exterior_curvature()
    for key, value in exterior.items():
        print(f"  {key:36s}: {value}")

    print("\n2. Bernoulli saturation of the pressure deficit")
    bernoulli = analyze_bernoulli_singularity_saturation()
    for key, value in bernoulli.items():
        print(f"  {key:36s}: {value}")

    print("\n3. Rarefaction and information-decoupling cutoff")
    rarefaction = analyze_rarefaction_information_cutoff()
    for key, value in rarefaction.items():
        print(f"  {key:36s}: {value}")

    print("\n4. Conditional finite-core matching")
    completion = analyze_geodesic_completion_by_core_matching()
    for key, value in completion.items():
        print(f"  {key:36s}: {value}")

    print("\n4b. C2 core matching coefficient derivation")
    core_coeffs = derive_c2_core_matching_coefficients()
    for key, value in core_coeffs.items():
        print(f"  {key:36s}: {value}")

    print("\n5. Horizonless exterior, throat, and boundary status")
    boundary = analyze_horizon_throat_and_boundary()
    for key, value in boundary.items():
        print(f"  {key:36s}: {value}")

    print("\n6. Photon sphere, shadow, and golden-ratio ISCO")
    strong_observables = analyze_photon_shadow_isco()
    for key, value in strong_observables.items():
        print(f"  {key:36s}: {value}")

    print("\n7. Singularity-strength ledger")
    for item in singularity_strength_ledger():
        print(f"  - {item}")

    print("\n8. Regular-center ansatz cross-check")
    G_t, G_r, G_th, R_0, K_0, T_t, T_r = analyze_regular_center()
    print(f"G^t_t (გეომეტრიული სიმრუდე ცენტრში) = {G_t}")
    print(f"G^r_r (გეომეტრიული წნევა ცენტრში) = {G_r}")
    print(f"G^th_th (კუთხური სიმრუდე) = {G_th}")
    print(f"R (რიჩის სკალარი ცენტრში) = {R_0}")
    print(f"K (Kretschmann სკალარი ცენტრში) = {K_0}")
    
    print(f"\nT^t_t (ენერგიის სიმკვრივე r=0-ზე) = {T_t}")
    print(f"T^r_r (რადიალური წნევა r=0-ზე) = {T_r}")
    
    print("\n--- სამუშაო დასკვნები / export gate ---")
    print("1. რეგულარულობის ანზაცით (A=1+O(r^2), B>0) სიმრუდის ინვარიანტები, მათ შორის")
    print("   Kretschmann (K), ცენტრში სასრულია ამ local ansatz-ის ფარგლებში.")
    print("2. მედიუმის სტრეს-ტენზორი T^t_t და T^r_r ასევე სასრულია.")
    print("3. T^t_t და T^r_r ზოგად შემთხვევაში არ არიან ტოლი, ამიტომ ეს არ არის სუფთა")
    print("   de Sitter-ის ვაკუუმი (w=-1). ტერმინი MD ტექსტში შეიცვალა 'სასრულ-ენერგიული ბირთვით'.")
    print("4. Static exponential exterior branch იძლევა horizonless benchmark-ებს:")
    print("   r_ph=r_s, b_c=e*r_s, r_ISCO=phi_golden^2*r_s.")
    print("   ISCO გამოდის V_eff-ის მარგინალური სტაბილურობიდან, არა fitting-ით.")
    print("5. ფუძის ნაწილაკების rarefaction closure აჩვენებს: n_eff->0, ell_mfp->infinity,")
    print("   Gamma_coll->0 და Kn->infinity, თუ closure ansatz ფიზიკურად სწორია.")
    print("6. r_s/4-ის შიგნით medium-stress gradient ნიშანს იცვლის; ეს არის")
    print("   backreaction ledger, არა დამატებითი წნევითი ძალა მატერიაზე.")
    print("7. Kn=1 ზედაპირზე მიიღება finite core radius r_c=r_s/W(n_0*sigma*r_s).")
    print("8. C2 matching core აკერებს exponential exterior-ს value/slope/curvature დონეზე.")
    print("9. სრული compact-object proof ჯერ ბლოკირებულია field-equation source,")
    print("   junction stress, physical energy, rotation, QNM/echo და EHT ray tracing-ით.")

# ===================== CONSOLIDATED PHASE SECTIONS =====================


# ===================== merged from p05_compact.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 29: EHT shadow - static spherical compact-object benchmark in RFG.

Phase 18 now derives the exponential exterior:

    ds^2 = -exp(-r_s/r)c^2dt^2 + exp(r_s/r)(dr^2+r^2dOmega^2).

For null circular orbits:

    d/dr [exp(-2r_s/r)/r^2] = 0 -> r_ph = r_s.

The critical impact parameter is therefore:

    b_c^RFG = e*r_s,

while Schwarzschild gives:

    b_c^GR = (3*sqrt(3)/2)*r_s.

Thus the static spherical RFG benchmark gives a shadow diameter larger by
2e/(3sqrt(3))-1 = 4.63%.  This is not a full EHT model until rotation,
plasma/accretion emission and ray-traced image fitting are added.
"""

import math


M_SUN = 1.98847e30
G = 6.67430e-11
C = 299792458.0
MUAS_TO_RAD = math.pi / (180 * 3600 * 1e6)
PC = 3.0856775814913673e16
KPC = 1.0e3 * PC
MPC = 1.0e6 * PC


EHT_OBSERVATIONS = {
    "M87": {
        "shadow_diameter_uas": 42.0,
        "shadow_error_uas": 3.0,
        "mass_solar": 6.5e9,
        "distance_Mpc": 16.8,
    },
    "SgrA": {
        "shadow_diameter_uas": 51.8,
        "shadow_error_uas": 2.3,
        "mass_solar": 4.154e6,
        "distance_pc": 8178,
        "mass_distance_source": "GRAVITY Collaboration priors (2019/2022)",
    },
}


def schwarzschild_shadow_prediction(m_solar, distance_m):
    """GR Schwarzschild critical-curve diameter: theta = 3*sqrt(3)*r_s/D."""
    mass = m_solar * M_SUN
    r_s = 2 * G * mass / C**2
    diameter = 3 * math.sqrt(3) * r_s
    theta_uas = (diameter / distance_m) / MUAS_TO_RAD
    return {
        "r_s_meters": r_s,
        "critical_impact_parameter_m": 0.5 * diameter,
        "shadow_diameter_meters": diameter,
        "theta_uas": theta_uas,
    }


def rfg_shadow_prediction(m_solar, distance_m):
    """RFG exponential-exterior critical-curve diameter: theta = 2*e*r_s/D."""
    mass = m_solar * M_SUN
    r_s = 2 * G * mass / C**2
    b_c = math.e * r_s
    diameter = 2 * b_c
    theta_uas = (diameter / distance_m) / MUAS_TO_RAD
    return {
        "r_s_meters": r_s,
        "photon_sphere_r": r_s,
        "critical_impact_parameter_m": b_c,
        "shadow_diameter_meters": diameter,
        "theta_uas": theta_uas,
    }


def compare_with_observation(name, m_solar, distance_m, theta_obs, theta_err):
    """Compare GR and RFG static spherical shadow benchmarks to one observation."""
    gr = schwarzschild_shadow_prediction(m_solar, distance_m)
    rfg = rfg_shadow_prediction(m_solar, distance_m)
    ratio = rfg["theta_uas"] / gr["theta_uas"]

    return {
        "name": name,
        "r_s_meters": gr["r_s_meters"],
        "GR_prediction_uas": gr["theta_uas"],
        "RFG_prediction_uas": rfg["theta_uas"],
        "RFG_over_GR_ratio": ratio,
        "RFG_shadow_shift_percent": (ratio - 1.0) * 100.0,
        "EHT_observation_uas": theta_obs,
        "EHT_error_uas": theta_err,
        "GR_deviation_sigma": abs(gr["theta_uas"] - theta_obs) / theta_err,
        "RFG_deviation_sigma": abs(rfg["theta_uas"] - theta_obs) / theta_err,
    }


def distance_to_meters(obs):
    """Convert observation distance fields to meters."""
    if "distance_Mpc" in obs:
        return obs["distance_Mpc"] * MPC
    if "distance_kpc" in obs:
        return obs["distance_kpc"] * KPC
    if "distance_pc" in obs:
        return obs["distance_pc"] * PC
    raise KeyError("distance field not found")


def rfg_shadow_derivation_ledger():
    return [
        "exponential exterior: g_tt=-exp(-r_s/r), g_rr=exp(r_s/r)",
        "null barrier: V_null proportional to exp(-2r_s/r)/r^2",
        "photon sphere: dV_null/dr=0 -> r_ph=r_s",
        "critical impact parameter: b_c=r*exp(r_s/r) at r=r_s -> e*r_s",
        "GR reference: b_c=(3*sqrt(3)/2)*r_s",
        "static spherical benchmark: RFG shadow diameter is +4.63% relative to GR",
    ]


def predictions_summary():
    """RFG vs GR shadow status."""
    ratio = 2.0 * math.e / (3.0 * math.sqrt(3.0))
    return {
        "current_status": "static spherical benchmark derived inside the exponential branch; full EHT verdict open",
        "RFG_b_c": "e*r_s",
        "GR_b_c": "3*sqrt(3)*r_s/2",
        "RFG_over_GR": ratio,
        "shift_percent": (ratio - 1.0) * 100.0,
        "needed_for_decisive_test": "rotating RFG exterior, accretion/plasma emission, mass-distance priors, and ray-traced image modelling",
        "ngEHT_BHEX_window": "few-percent shadow/ring precision can test the +4.63% benchmark",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 29: EHT shadow - RFG b_c=e*r_s static benchmark")
    print("=" * 72)

    print("\n1. დაკვირვება (EHT priors used in this local script)")
    for name, obs in EHT_OBSERVATIONS.items():
        print(f"\n  {name}")
        for key, val in obs.items():
            print(f"    {key:25s}: {val}")

    print("\n2. Derivation ledger")
    for item in rfg_shadow_derivation_ledger():
        print(f"  - {item}")

    print("\n3. GR vs RFG static benchmark compared to EHT numbers")
    for name, obs in EHT_OBSERVATIONS.items():
        result = compare_with_observation(
            name,
            obs["mass_solar"],
            distance_to_meters(obs),
            obs["shadow_diameter_uas"],
            obs["shadow_error_uas"],
        )
        print(f"\n  {result['name']}")
        print(f"    r_s = {result['r_s_meters']:.3e} m")
        print(f"    GR theta  = {result['GR_prediction_uas']:.2f} microas")
        print(f"    RFG theta = {result['RFG_prediction_uas']:.2f} microas")
        print(f"    RFG/GR    = {result['RFG_over_GR_ratio']:.8f}")
        print(f"    shift     = {result['RFG_shadow_shift_percent']:.2f}%")
        print(f"    observed  = {result['EHT_observation_uas']:.1f} +/- {result['EHT_error_uas']:.1f} microas")
        print(f"    GR sigma  = {result['GR_deviation_sigma']:.2f}")
        print(f"    RFG sigma = {result['RFG_deviation_sigma']:.2f}")

    print("\n4. Predictions summary")
    for key, val in predictions_summary().items():
        print(f"  {key:26s}: {val}")

    print("\n5. სტატუსი")
    print("  - +4.6% shadow shift static exponential benchmark-ად გამოდის.")
    print("  - მიმდინარე EHT რიცხვები არ არის საკმარისი სუფთა GR/RFG გარჩევისთვის.")
    print("  - decisive test მოითხოვს rotating RFG ray tracing-ს, plasma model-ს და ngEHT/BHEX კლასის სიზუსტეს.")


# ===================== merged from p05_compact.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 31: ნეიტრონული ვარსკვლავები — RFG ანიზოტროპიული TOV და M-R მრუდი
================================================================================

რეფერენცია: p10_oscillons.py, p01_core.py,
            p01_core.py, STRATEGY.md S3/E6

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

# merged future import removed: from __future__ import annotations

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
    supports_2p08_toy: bool
    lambda_bound_toy: bool | None
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
        supports_2p08_toy=max_solution.mass_solar >= OBSERVATIONS["NICER_J0740"]["mass_solar"],
        lambda_bound_toy=(onep4[1] < lambda_bound if onep4 else None),
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
        f"toy M_max >= 2.08 M_sun diagnostic: {summary.supports_2p08_toy}",
        f"toy Lambda_1.4 proxy < 580 diagnostic: {summary.lambda_bound_toy}",
        f"max c_s/c in sequence: {summary.max_sound_speed_over_c:.3f}",
    ]


def model_scope_notes() -> list[str]:
    return [
        "The polytropic EOS is a controlled toy EOS, not a nuclear-matter fit.",
        "eta_delta is a phenomenological stand-in for p01_core Delta p.",
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


# ===================== OLD COMPACT INTEGRATION =====================

"""
STAGE A3: OLD compact-object boundary/interior/verification ledger

Sources drained from:
    OLD/12. ISPG_Geodesics.tex
    OLD/14. ISPG_Interior.tex
    OLD/15. ISPG_Verification_Task1a.tex

OLD material is represented here as candidate ledger material, not as final
authority.  The exponential exterior, rarefaction/core matching, horizonless
boundary, photon/shadow and golden-ISCO blocks still carry the open checks
listed in the central gate.
"""


def stage_a3_old_geodesic_interior_drain():
    return {
        "OLD/12_geodesics": {
            "status": "represented_as_candidate_ledger",
            "new_functions": [
                "analyze_horizon_throat_and_boundary",
                "analyze_photon_shadow_isco",
                "analyze_geodesic_completion_by_core_matching",
            ],
            "core_results": [
                "photon sphere r_ph = r_s",
                "critical impact parameter b_c = e*r_s",
                "proper distance to old r=0 boundary diverges",
                "curvature invariants vanish as r -> 0",
                "old exterior geodesic incompleteness is conditionally addressed by Knudsen-core matching",
            ],
        },
        "OLD/14_interior": {
            "status": "represented_as_candidate_ledger",
            "new_functions": [
                "analyze_rarefaction_information_cutoff",
                "analyze_bernoulli_singularity_saturation",
                "singularity_strength_ledger",
            ],
            "core_results": [
                "clock/freezing: d tau/dt = exp(-r_s/(2r)) -> 0",
                "proper radial element exp(r_s/(2r)) dr diverges",
                "effective signal speed c_eff = c exp(-r_s/r) -> 0",
                "mass-as-pressure-deficit is retained as a candidate interpretation",
                "deep deficit is self-perpetuating through divergent relaxation time",
            ],
        },
    }


def stage_a3_scalar_perturbation_verification():
    """
    OLD/15 Task 1a probe-level ledger in executable symbolic form.

    Probe scalar perturbations on phi0=-r_s/r reduce to a Schrodinger-form
    equation after the tortoise coordinate and first-derivative removal.
    The dangerous l=0 well exists only in 0<r<r_s/4 and is exponentially
    suppressed.  This block is still a probe scalar test; the coupled
    even/polar metric-medium sector is open.
    """
    r, r_s, c, ell = sp.symbols("r r_s c ell", positive=True)
    x = sp.Symbol("x", positive=True)
    phi0 = -r_s / r
    phi0_p = sp.diff(phi0, r)
    phi0_pp = sp.diff(phi0_p, r)

    c_eff_sq = c**2 * sp.exp(2 * phi0)
    drstar_dr = sp.exp(-phi0) / c
    dr_drstar = c * sp.exp(phi0)
    P = -c * sp.exp(phi0) * phi0_p
    Q = c**2 * sp.exp(2 * phi0) * ell * (ell + 1) / r**2
    V_eff = sp.simplify(
        c**2
        * sp.exp(2 * phi0)
        * (ell * (ell + 1) / r**2 - sp.Rational(1, 2) * phi0_pp - sp.Rational(1, 4) * phi0_p**2)
    )
    V_l0 = sp.factor(V_eff.subs(ell, 0))
    V_l0_dimensionless = sp.simplify(V_l0.subs({c: 1, r_s: 1, r: x}))
    drstar_dx_dimensionless = sp.exp(1 / x)
    well_integrand_dr = sp.simplify(-V_l0_dimensionless)
    well_integrand_drstar = sp.simplify(well_integrand_dr * drstar_dx_dimensionless)
    well_integral_dr = sp.integrate(well_integrand_dr, (x, 0, sp.Rational(1, 4)))
    well_integral_drstar = sp.integrate(well_integrand_drstar, (x, 0, sp.Rational(1, 4)))

    return {
        "phi0_prime": sp.Eq(sp.Symbol("phi0_prime"), phi0_p),
        "phi0_second": sp.Eq(sp.Symbol("phi0_second"), phi0_pp),
        "c_eff_squared": sp.Eq(sp.Symbol("c_eff^2"), c_eff_sq),
        "c_eff_positive": "c_eff^2 = c^2 exp(-2 r_s/r) > 0 for r>0; no gradient instability",
        "tortoise_drstar_dr": sp.Eq(sp.Symbol("drstar_dr"), drstar_dr),
        "tortoise_dr_drstar": sp.Eq(sp.Symbol("dr_drstar"), dr_drstar),
        "first_derivative_coefficient": sp.Eq(sp.Symbol("P"), P),
        "angular_Q": sp.Eq(sp.Symbol("Q"), Q),
        "V_eff": sp.Eq(sp.Symbol("V_eff"), V_eff),
        "V_l0_factorized": sp.Eq(sp.Symbol("V_l0"), V_l0),
        "negative_region_l0": "0 < r < r_s/4",
        "V_l0_dimensionless_c_rs_1": V_l0_dimensionless,
        "well_integrand_absV_dr_c_rs_1": well_integrand_dr,
        "well_integral_absV_dr_c_rs_1": sp.Eq(sp.Symbol("I_dr"), well_integral_dr),
        "well_integrand_absV_drstar_c_rs_1": well_integrand_drstar,
        "well_integral_absV_drstar_c_rs_1": sp.Eq(sp.Symbol("I_drstar"), well_integral_drstar),
        "well_integral_absV_drstar_numeric": sp.N(well_integral_drstar, 8),
        "probe_verdict": "small l=0 well at probe level; not a coupled stability theorem",
        "ell_ge_1_status": "centrifugal barrier is favorable, but no hard-coded bound is exported here",
        "remaining_scope": "even/polar coupled metric-scalar sector remains a future verification task",
    }


def stage_a3_compact_old_file_status():
    return {
        "geodesic_interior": stage_a3_old_geodesic_interior_drain(),
        "scalar_perturbations": stage_a3_scalar_perturbation_verification(),
        "integration_status": "Stage A3 represented as candidate ledger in p05_compact.py; theory export remains gated",
    }


def compact_central_claim_gate():
    """
    One-place export gate for p05_compact.py.

    This keeps the useful algebraic results while blocking the stronger
    compact-object claim until the missing physical/dynamical checks are done.
    """
    signature = compact_signature_bridge()
    scalar_probe = stage_a3_scalar_perturbation_verification()
    return {
        "file_export_status": "NOT_READY_FOR_RFG_THEORY_EXPORT",
        "signature_bridge": signature["stress_bridge_status"],
        "exterior_status": "DERIVED_FROM_VACUUM_PHASE_EQUATION_PLUS_BICONFORMAL_MAP",
        "black_hole_breaker_status": "SCHWARZSCHILD_CURVATURE_SINGULARITY_REMOVED_AT_GEOMETRY_LEVEL__GEODESIC_BOUNDARY_STILL_OPEN",
        "effective_source_status": "GEOMETRIC_SOURCE_PROFILE_MATCHES_BERNOULLI_DELTA_P",
        "p01_source_closure": "F_EQ_R_BRANCH_FAILS__LINEAR_NONTRIVIAL_F_R_BRANCH_DERIVED__EXACT_FULL_STRESS_SOLUTION_OPEN",
        "core_status": "C2_MATCHING_CONDITIONAL_ANSATZ__JUNCTION_STRESS_AND_STABILITY_OPEN",
        "energy_status": "COORDINATE_ENERGY_FINITE__PHYSICAL_PROPER_ADM_ENERGY_OPEN",
        "scalar_stability": scalar_probe["probe_verdict"],
        "shadow_status": "STATIC_SPHERICAL_BENCHMARK_ONLY__ROTATION_PLASMA_RAYTRACING_OPEN",
        "neutron_star_status": "TOY_POLYTROPE_AND_LAMBDA_PROXY_ONLY__EOS_LOVE_BAYESIAN_FIT_OPEN",
        "observational_blockers": [
            "rotating compact-object exterior",
            "EHT ray-traced images with mass-distance and plasma/accretion priors",
            "full coupled QNM/ringdown and echo transfer function",
            "surface/absorption luminosity for horizonless core candidates",
            "RFG-derived neutron-star EOS/anisotropy and Love-number ODEs",
        ],
        "do_not_claim": [
            "do not claim a derived RFG black-hole replacement from the static ansatz alone",
            "do not claim geodesic completion from C2 matching alone",
            "do not claim no-horizon observational viability before QNM/echo/surface tests",
            "do not claim EHT support from the static +4.63% benchmark alone",
            "do not claim NICER/GW170817 pass from toy TOV diagnostics",
        ],
    }


if __name__ == "__main__":
    print("=" * 72)
    print("STAGE A3: OLD compact/interior/verification integration ledger")
    print("=" * 72)

    status = stage_a3_compact_old_file_status()

    print("\n1. OLD/12 and OLD/14 drain status")
    for old_file, info in status["geodesic_interior"].items():
        print(f"  {old_file}: {info['status']}")
        for item in info["core_results"]:
            print(f"    - {item}")

    print("\n2. OLD/15 scalar perturbation verification")
    for key, value in status["scalar_perturbations"].items():
        print(f"  {key:28s}: {value}")

    print("\n3. Integration verdict")
    print("  - OLD/12, OLD/14 and OLD/15 are represented here as candidate ledger")
    print("    material. They are not final theory authority until the central")
    print("    compact-object gates are closed.")

    print("\n4. Central compact-object claim gate")
    for key, value in compact_central_claim_gate().items():
        print(f"  {key:30s}: {value}")

