# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 18: RFG compact object, Bernoulli saturation, and singularity audit

ეს ფაილი აძლიერებს ძველი თეორიის სინგულარობის ბირთვს:
1. exponential bi-conformal exterior:
       ds^2 = -exp(-r_s/r)c^2dt^2 + exp(r_s/r)(dr^2+r^2dOmega^2)
2. curvature invariants vanish at r -> 0;
3. Bernoulli pressure deficit is self-limiting;
4. finite-radius Killing horizon is absent;
5. areal throat, photon sphere, shadow radius, and golden-ratio ISCO follow;
6. exterior-only geodesic incompleteness is kept explicit;
7. a rarefaction cutoff plus C2 finite-core matching closes the endpoint.

No new action, no new field, no extra force law.
"""
import sympy as sp
from phase1_action import get_polynomial_lagrangian


def analyze_exponential_exterior_curvature():
    """
    ძველი compact-object branch-ის curvature theorem.

    metric:
        g_tt = -exp(-r_s/r)
        g_rr = exp(+r_s/r)

    ძველი თეორიის closed-form invariants:
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
        "theorem": "exponential exterior removes the curvature blow-up",
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
    Bernoulli pressure-deficit explanation of singularity avoidance.

    Delta P = exp(phi)(phi')^2/(32*pi*G)
    for phi=-r_s/r gives
        Delta P = r_s^2 exp(-r_s/r)/(32*pi*G*r^4).

    It peaks at r_s/4 and vanishes at r->0: the deficit saturates.
    """
    r, r_s, G = sp.symbols('r r_s G', real=True, positive=True)
    phi = -r_s / r
    delta_p = sp.simplify(sp.exp(phi) * sp.diff(phi, r)**2 / (32 * sp.pi * G))
    p_static = -delta_p
    coordinate_energy = sp.simplify(
        sp.integrate(delta_p * 4 * sp.pi * r**2, (r, 0, sp.oo))
    )
    u = sp.Symbol('u', real=True, positive=True)
    shape = sp.exp(-u) * u**4

    return {
        "theorem": "Bernoulli pressure deficit is self-limiting",
        "Delta_P": sp.Eq(sp.Symbol('Delta_P'), delta_p),
        "P_static": sp.Eq(sp.Symbol('P_static'), p_static),
        "Bernoulli_identity": sp.Eq(sp.Symbol('P_static + Delta_P'), 0),
        "dimensionless_shape": sp.Eq(sp.Symbol('shape(u)'), shape),
        "shape_derivative": sp.factor(sp.diff(shape, u)),
        "pressure_peak": sp.Eq(sp.Symbol('r_peak'), r_s / 4),
        "lim_r_to_0_Delta_P": sp.limit(delta_p, r, 0, dir='+'),
        "lim_r_to_inf_Delta_P": sp.limit(delta_p, r, sp.oo),
        "finite_coordinate_energy": sp.Eq(sp.Symbol('int_DeltaP_d3x'), coordinate_energy),
        "meaning": "as phi -> -infinity, exp(phi) shuts off the gradient-energy transfer; no curvature blow-up forms.",
    }


def analyze_rarefaction_information_cutoff():
    """
    Microscopic closure for the r=0 boundary.

    User intuition, written as a continuum-mechanics criterion:
    in the deepest Bernoulli deficit the active carrier density of the
    vacuum medium is rarefied, collisions become sparse, and the continuum
    no longer transmits information as a connected elastic fluid.

    Minimal closure:
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
        "theorem": "rarefaction closes the information channel before a material singularity forms",
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
        "outward_pressure_force_density": sp.Eq(sp.Symbol('f_pressure'), outward_pressure_force_density),
        "force_turning_radius": sp.Eq(sp.Symbol('r_turn'), r_s / 4),
        "inner_core_sign": "for r<r_s/4, d(Delta_P)/dr>0: the pressure-gradient reaction points outward.",
        "physical_meaning": "near r=0 the medium becomes dilute and non-communicating; oscillons contain vanishing carrier number in any finite unresolved cell.",
    }


def analyze_geodesic_completion_by_core_matching():
    """
    Boundary-completion theorem.

    The exponential exterior is not forced to run all the way to r=0.
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
    center is regular, so geodesics hitting the old boundary continue through
    the finite core.
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
        "theorem": "finite-core matching completes the old r=0 boundary",
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
        "geodesic_completion_rule": "replace the formal r=0 endpoint by r<=r_c regular core; bounded Christoffels imply local continuation of null/timelike geodesics.",
        "physical_meaning": "the singular endpoint is not part of the continuum phase; it is pre-empted by a dilute kinetic core fixed by Kn=1.",
    }


def analyze_horizon_throat_and_boundary():
    """
    Horizonless exterior-only boundary ledger.

    finite r>0-ზე g_tt never vanishes; r=0 is a boundary, not a finite
    horizon. Proper distance to r=0 diverges, while captured geodesics
    reach the boundary at finite affine/proper parameter before the
    rarefied-core completion is imposed.
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
        "theorem": "no finite-radius Killing horizon; r=0 is a frozen boundary",
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
        "geodesic_status": "exterior-only captured geodesics reach r=0 in finite affine parameter; section 4 replaces that endpoint by a regular Knudsen core.",
    }


def analyze_photon_shadow_isco():
    """
    Strong-field observables of the exponential exterior.

    photon sphere:
        d/dr [exp(-2r_s/r)/r^2]=0 -> r=r_s
    shadow:
        b_c = e r_s
    massive ISCO:
        r^2 - 3 r_s r + r_s^2 = 0 -> r_ISCO = phi_golden^2 r_s
    """
    r, r_s, c = sp.symbols('r r_s c', real=True, positive=True)
    phi_golden = (1 + sp.sqrt(5)) / 2

    photon_barrier = sp.exp(-2 * r_s / r) / r**2
    photon_condition = sp.factor(sp.diff(photon_barrier, r))

    isco_poly = r**2 - 3 * r_s * r + r_s**2
    isco_roots = sp.solve(sp.Eq(isco_poly, 0), r)
    isco_physical = sp.simplify(isco_roots[1])

    omega_sq = c**2 * r_s * sp.exp(-2 * r_s / r) / (r**2 * (2 * r - r_s))
    omega_isco_sq = sp.simplify(omega_sq.subs(r, isco_physical))

    gr_isco_iso = (5 + 2 * sp.sqrt(6)) * r_s / 4
    shadow_rfg = sp.E * r_s
    shadow_gr = 3 * sp.sqrt(3) * r_s / 2
    shadow_ratio = sp.N(shadow_rfg / shadow_gr, 8)

    return {
        "photon_barrier": photon_barrier,
        "photon_condition": sp.Eq(photon_condition, 0),
        "photon_sphere": sp.Eq(sp.Symbol('r_ph'), r_s),
        "critical_impact_parameter": sp.Eq(sp.Symbol('b_c'), shadow_rfg),
        "GR_shadow_reference": sp.Eq(sp.Symbol('b_c_GR'), shadow_gr),
        "shadow_ratio_RFG_over_GR": shadow_ratio,
        "shadow_size_shift": f"{float((shadow_ratio - 1) * 100):.2f}%",
        "ISCO_polynomial": sp.Eq(isco_poly, 0),
        "ISCO_roots": isco_roots,
        "ISCO_physical": sp.Eq(sp.Symbol('r_ISCO'), isco_physical),
        "golden_ratio_identity": sp.Eq(sp.Symbol('r_ISCO'), phi_golden**2 * r_s),
        "GR_isotropic_ISCO_reference": sp.Eq(sp.Symbol('r_ISCO_GR_iso'), gr_isco_iso),
        "Omega_ISCO_squared": sp.Eq(sp.Symbol('Omega_ISCO^2'), omega_isco_sq),
        "frequency_proxy": "old-theory benchmark: f_ISCO approximately 0.931 f_ISCO_GR for the same total mass",
    }


def singularity_strength_ledger() -> list[str]:
    return [
        "Curvature blow-up is removed: R->0, Ricci^2->0, K->0 at r->0.",
        "Bernoulli pressure deficit saturates: Delta_P peaks at r_s/4 and returns to 0 at r->0.",
        "Microscopic rarefaction gives n_eff->0, mean free path->infinity, and collision/information rate->0 at r->0.",
        "The Knudsen number Kn=ell_mfp/r diverges, so the continuum approximation self-terminates before an infinite-density core can form.",
        "Inside r_s/4 the Bernoulli pressure-gradient reaction points outward, giving a built-in core self-regulation channel.",
        "The formal boundary is replaced at Kn=1 by r_c=r_s/W(n_0*sigma*r_s), a finite positive matching radius.",
        "A positive C2 logarithmic core matches the exponential exterior through value, slope, and curvature at r_c.",
        "The matched core has A(0)=1, B(0)>0, A'(0)=B'(0)=0 and finite R_0, K_0, so geodesics continue through a regular center.",
        "There is no finite-radius Killing horizon: exp(-r_s/r)>0 for every r>0.",
        "The areal radius has a throat at r_s/2, with R_min=e*r_s/2.",
        "Photon sphere is r_s and the critical shadow impact parameter is b_c=e*r_s.",
        "Massive-particle ISCO is r_ISCO=phi_golden^2*r_s.",
        "External observers see infinite redshift/coordinate-time freezing toward r=0.",
        "Captured geodesics that reached the old r=0 endpoint now enter the finite rarefied core instead.",
        "Therefore the strengthened claim is curvature-regular, horizonless, Bernoulli-saturated, and geodesically completed by Knudsen-core matching.",
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

    print("\n1. Exponential exterior curvature theorem")
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

    print("\n4. Geodesic completion by finite-core matching")
    completion = analyze_geodesic_completion_by_core_matching()
    for key, value in completion.items():
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
    
    print("\n--- აგენტთა საბჭოს დასკვნები ---")
    print("1. რეგულარულობის ანზაცით (A=1+O(r^2), B>0) სიმრუდის ინვარიანტები, მათ შორის")
    print("   უმკაცრესი Kretschmann (K) სკალარი, ცენტრში აბსოლუტურად სასრულია (არა სინგულარული)!")
    print("2. მედიუმის სტრეს-ტენზორი T^t_t და T^r_r ასევე სასრულია.")
    print("3. T^t_t და T^r_r ზოგად შემთხვევაში არ არიან ტოლი, ამიტომ ეს არ არის სუფთა")
    print("   de Sitter-ის ვაკუუმი (w=-1). ტერმინი MD ტექსტში შეიცვალა 'სასრულ-ენერგიული ბირთვით'.")
    print("4. ძველი exponential exterior უკვე იძლევა horizonless compact-object პროგნოზებს:")
    print("   r_ph=r_s, b_c=e*r_s, r_ISCO=phi_golden^2*r_s.")
    print("5. ფუძის ნაწილაკების rarefaction closure აჩვენებს: n_eff->0, ell_mfp->infinity,")
    print("   Gamma_coll->0 და Kn->infinity, ამიტომ ინფორმაციის/დარტყმითი გადაცემა ითიშება.")
    print("6. r_s/4-ის შიგნით ბერნულის წნევის გრადიენტის რეაქცია გარეთკენ არის, რაც")
    print("   ბუნებრივ core self-regulation/back-reaction არხს იძლევა.")
    print("7. Kn=1 ზედაპირზე მიიღება finite core radius r_c=r_s/W(n_0*sigma*r_s).")
    print("8. აშენებული C2 matching core ზუსტად აკერებს exponential exterior-ს და ცენტრში")
    print("   იძლევა A(0)=1, B(0)>0, სასრულ R_0-ს და სასრულ K_0-ს.")
    print("9. ამიტომ ძველი r=0 endpoint იცვლება რეგულარული rarefied core-ით:")
    print("   geodesic incompleteness იხურება core-matching completion-ით.")
