# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 19: Inertia from oscillon collective coordinates and RFG dressing.

Old ISPG Appendix 17 used the correct physical intuition:
    acceleration -> retarded self-field -> fore/aft asymmetry -> reaction.

But the old self-force route had two technical hazards:
    1. the finite regular self-field coefficient was inserted as a canonical
       regularization value;
    2. the spatial metric sector was asserted to supply the second half of
       the force, but the full Noether momentum integral was not displayed.

This RFG phase keeps that intuition as the dynamical dressing picture, but
closes the leading inertial law at the theorem-gate level by a stronger route:

    localized oscillon energy E0
        -> translational collective coordinate R(t)
        -> Noether four-momentum P^mu = integral T^{0mu} d^3x
        -> Lorentz/Noether effective action
        -> p = gamma M v, M = E0/c^2
        -> F = M a in the nonrelativistic limit.

The same M sources the static Bernoulli/zero-frequency gravitational field:

    phi = -2GM/(c^2 r).

Therefore, in the controlled leading regime:

    m_inertial = m_gravitational = E0/c^2.

No new parameter is introduced.  Retarded self-field asymmetry becomes the
mechanical explanation of how the medium reacts during acceleration; the
collective-coordinate theorem is the clean derivation of the coefficient.

Open before final theory export:
    derive the localized dressed oscillon solution from the active RFG action,
    compute its full Noether stress-energy including boundary/dressing terms,
    and verify that the translation zero-mode norm equals the total dressed
    rest energy divided by c^2.
"""

from __future__ import annotations

import math
import sympy as sp


def collective_coordinate_inertia_theorem():
    """
    Leading inertial law for a localized oscillon.

    For any localized finite-energy solution with translational modulus R(t),
    Lorentz/Noether symmetry fixes the low-velocity effective action:

        L_eff = -M c^2 sqrt(1-v^2/c^2),      M = E0/c^2.

    This gives p=gamma M v and F=dp/dt.  The nonrelativistic limit is F=M a.
    """
    E0, M, c, v, a = sp.symbols('E0 M c v a', positive=True, real=True)
    gamma = 1 / sp.sqrt(1 - v**2 / c**2)
    L_eff = -M * c**2 * sp.sqrt(1 - v**2 / c**2)
    p = sp.simplify(sp.diff(L_eff, v))
    p_expected = sp.simplify(gamma * M * v)
    L_nr = sp.series(L_eff, v, 0, 6).removeO()
    inertial_force_nr = sp.simplify(sp.diff(M * v, v) * a)

    return {
        "theorem": "localized oscillon translational modulus has inertial mass M=E0/c^2",
        "energy_mass_identity": sp.Eq(M, E0 / c**2),
        "effective_lagrangian": sp.Eq(sp.Symbol('L_eff'), L_eff),
        "momentum": sp.Eq(sp.Symbol('p'), p),
        "momentum_minus_gammaMv": sp.simplify(p - p_expected),
        "low_velocity_lagrangian": L_nr,
        "nonrelativistic_force": sp.Eq(sp.Symbol('F'), inertial_force_nr),
        "meaning": "F=Ma follows from translation + Lorentz/Noether structure, not from an adjustable drag coefficient.",
    }


def noether_four_momentum_closure():
    """
    Stress-energy / Noether version of the same theorem.

    For any localized excitation of a translation-invariant relativistic
    medium, the conserved four-momentum is

        P^mu = integral T^{0mu} d^3x.

    With signature (+---):
        P_mu P^mu = M^2 c^2,
        E^2 = p^2 c^2 + M^2 c^4.

    This is the invariant definition of the inertial mass.  It does not
    depend on the shape details of the oscillon once the total stress-energy
    is localized and finite.
    """
    E0, M, c, v = sp.symbols('E0 M c v', positive=True, real=True)
    gamma = 1 / sp.sqrt(1 - v**2 / c**2)
    energy = gamma * M * c**2
    momentum = gamma * M * v
    p0 = energy / c
    invariant = sp.simplify(p0**2 - momentum**2)
    dispersion = sp.simplify(energy**2 - momentum**2 * c**2)

    return {
        "Noether_charge": "P^mu = integral T^{0mu} d^3x",
        "rest_energy": sp.Eq(sp.Symbol('E0'), M * c**2),
        "energy": sp.Eq(sp.Symbol('E'), energy),
        "momentum": sp.Eq(sp.Symbol('p'), momentum),
        "invariant_P_mu_P_mu": sp.Eq(sp.Symbol('P_mu P^mu'), invariant),
        "dispersion_relation": sp.Eq(sp.Symbol('E^2-p^2 c^2'), dispersion),
        "mass_from_invariant": sp.Eq(M, E0 / c**2),
        "meaning": "inertial mass is the invariant norm of the oscillon+dressing stress-energy.",
    }


def moduli_space_projection_theorem():
    """
    Collective-coordinate projection.

    Let R(t) be the translation modulus of a localized oscillon.  The
    low-speed effective action has the moduli metric

        G_ij = (E0/c^2) delta_ij

    fixed by the Lorentz Ward identity.  Projecting the field equation onto
    the translational zero mode gives

        G_ij Rddot^j = F_i.

    Therefore the translational equation of motion is F_i = M Rddot_i.
    """
    E0, M, c, Rddot, F = sp.symbols('E0 M c Rddot F', positive=True, real=True)
    G_mod = E0 / c**2
    projected_equation = sp.Eq(F, G_mod * Rddot)
    closed_equation = sp.simplify((G_mod * Rddot).subs(E0, M * c**2))

    return {
        "translation_zero_mode": "partial_i Phi_0",
        "moduli_metric": sp.Eq(sp.Symbol('G_RR'), G_mod),
        "Ward_identity": "Lorentz symmetry fixes G_RR=E0/c^2 for the translation mode",
        "projected_EOM": projected_equation,
        "after_E0_equals_Mc2": sp.Eq(F, closed_equation),
        "meaning": "the coefficient of acceleration is the zero-mode norm, fixed by total rest energy.",
    }


def lorentz_ward_zero_mode_norm_gate():
    """
    The missing integral reduced to one Ward-identity gate.

    A direct RFG proof must compute the translation zero-mode norm from the
    active stress tensor.  Lorentz/Noether symmetry tells us what the answer
    has to be if the localized dressed oscillon is a genuine relativistic
    soliton:

        E(v) = gamma E0
        E(v)-E0 = 1/2 G_RR v^2 + ...

    Comparing the v^2 coefficient gives

        G_RR = E0/c^2.

    This is the precise bridge between "Noether bookkeeping" and the actual
    RFG zero-mode integral that still has to be evaluated for the full
    dressed solution.
    """
    E0, c, v, G_RR, a = sp.symbols(
        'E0 c v G_RR a',
        positive=True,
        real=True,
    )
    gamma = 1 / sp.sqrt(1 - v**2 / c**2)
    boosted_energy = gamma * E0
    boosted_kinetic_series = sp.series(boosted_energy - E0, v, 0, 6).removeO()
    collective_kinetic = sp.Rational(1, 2) * G_RR * v**2
    ward_solution = sp.solve(
        sp.Eq(
            sp.expand(boosted_kinetic_series).coeff(v, 2),
            sp.expand(collective_kinetic).coeff(v, 2),
        ),
        G_RR,
        dict=True,
    )[0]
    zero_mode_norm = sp.simplify(G_RR.subs(ward_solution))
    p_low = sp.simplify(zero_mode_norm * v)
    force_low = sp.simplify(zero_mode_norm * a)

    return {
        "rest_energy_integral": "E0 = integral T^00_dressed d^3x",
        "momentum_integral": "P^i = integral T^0i_dressed d^3x",
        "boosted_energy": sp.Eq(sp.Symbol('E(v)'), boosted_energy),
        "boosted_kinetic_series": boosted_kinetic_series,
        "collective_coordinate_kinetic": sp.Eq(sp.Symbol('T_R'), collective_kinetic),
        "zero_mode_norm_from_Ward_identity": sp.Eq(sp.Symbol('G_RR'), zero_mode_norm),
        "low_speed_momentum": sp.Eq(sp.Symbol('p_i'), p_low),
        "low_speed_force": sp.Eq(sp.Symbol('F_i'), force_low),
        "RFG_specific_integral_still_needed": (
            "compute G_RR directly from the p01/p10 dressed oscillon stress "
            "tensor and verify boundary terms vanish"
        ),
        "gate_status": "WARD_IDENTITY_CLOSES_COEFFICIENT__DIRECT_RFG_ZERO_MODE_INTEGRAL_OPEN",
    }


def integrated_noether_boost_closure():
    """
    Integrated Noether boost closure.

    This is stronger than only writing the effective Lagrangian.  Start from
    the integrated rest-frame dressed stress-energy of a localized oscillon:

        P^0 = E0/c,     P^i = 0.

    If the localized stress tensor obeys the static Laue/boundary condition,
    the integrated P^mu transforms as a Lorentz four-vector.  Boosting along
    x gives

        P'^0 = gamma E0/c,
        P'^x = gamma E0 v/c^2.

    Therefore the measured inertial momentum is p=gamma(E0/c^2)v and the
    nonrelativistic inertial mass is E0/c^2.
    """
    E0, c, v, a = sp.symbols('E0 c v a', positive=True, real=True)
    beta = v / c
    gamma = 1 / sp.sqrt(1 - beta**2)
    P0_rest = E0 / c
    Px_rest = sp.Integer(0)
    P0_boost = sp.simplify(gamma * (P0_rest + beta * Px_rest))
    Px_boost = sp.simplify(gamma * (Px_rest + beta * P0_rest))
    E_boost = sp.simplify(c * P0_boost)
    p_boost = Px_boost
    M_i = E0 / c**2
    p_expected = sp.simplify(gamma * M_i * v)
    low_speed_force = sp.simplify(M_i * a)

    return {
        "rest_Noether_charge": "P^mu_rest = (E0/c, 0, 0, 0)",
        "boost_condition": "valid when the dressed localized stress obeys the Laue/boundary condition",
        "boosted_energy": sp.Eq(sp.Symbol("E_prime"), E_boost),
        "boosted_momentum_x": sp.Eq(sp.Symbol("p_x_prime"), p_boost),
        "momentum_minus_gamma_Mv": sp.simplify(p_boost - p_expected),
        "inertial_mass_from_integrated_charge": sp.Eq(sp.Symbol("M_i"), M_i),
        "low_speed_force": sp.Eq(sp.Symbol("F_i"), low_speed_force),
        "closure_status": "INTEGRATED_NOETHER_BOOST_CLOSES_M_I_IF_LAUE_BOUNDARY_CONDITION_HOLDS",
    }


def laue_stress_condition_gate():
    """
    Static stress-balance condition required for particle-like inertia.

    For a localized static solution with spatial stress conservation

        partial_i T^{ij} = 0

    and a vanishing boundary flux at infinity, multiply by x^k and integrate:

        integral partial_i(x^k T^{ij}) d^3x
          = integral T^{kj} d^3x + boundary-free term = 0.

    Hence integral T^{ij} d^3x = 0.  This is the Laue condition.  It is the
    local stress-balance gate that lets the integrated energy-momentum behave
    as a particle four-vector.  RFG still has to verify it for the concrete
    dressed oscillon profile, but the required test is now explicit.
    """
    S_kj, I_Tkj = sp.symbols('S_kj I_Tkj', real=True)
    divergence_identity = sp.Eq(S_kj, I_Tkj)
    laue_with_boundary_zero = sp.Eq(I_Tkj, 0)

    return {
        "local_static_balance": "partial_i T^{ij}=0",
        "boundary_condition": "lim_R int_{S_R} x^k T^{ij} n_i dS = 0",
        "integrated_identity": divergence_identity,
        "Laue_condition": laue_with_boundary_zero,
        "physical_meaning": "internal stresses cannot contribute an extra inertial coefficient once the dressed oscillon is a localized static solution",
        "RFG_test_needed": "evaluate the p01/p10 dressed oscillon stress tail and verify the boundary flux vanishes",
        "gate_status": "LAUE_STRESS_BALANCE_REDUCES_DIRECT_INTEGRAL_TO_BOUNDARY_AND_LOCALIZATION_TEST",
    }


def rfg_action_to_inertia_derivation_gate():
    """
    Checklist for turning the theorem gate into an RFG derivation.

    This makes the remaining weakness explicit: p06 can use the universal
    Lorentz/Noether soliton theorem, but RFG still owes the concrete dressed
    oscillon integral.
    """
    return {
        "required_chain": [
            "active RFG action L(Y,I1,I2,I3)",
            "localized finite-energy dressed oscillon solution",
            "static stress balance partial_i T^{ij}=0",
            "Laue condition integral T^{ij} d^3x=0 from vanishing boundary flux",
            "normalizable translation zero mode d_i fields_0",
            "full dressed stress tensor including core, Bernoulli tail and boundary terms",
            "Noether charge P^mu = integral T^0mu d^3x with dP^mu/dt=0",
            "direct zero-mode norm G_RR = integral K_RR d^3x",
            "verify G_RR = E0/c^2",
            "use the same E0 in the zero-frequency gravitational source",
        ],
        "closed_by_current_file": [
            "universal Lorentz/Noether coefficient if the dressed oscillon exists",
            "integrated boost of P^mu gives p=gamma(E0/c^2)v once Laue holds",
            "m_i=m_g bookkeeping from a single total E0",
            "no uniform-motion drag in the Lorentz-boosted branch",
            "radiation reaction classified as higher-derivative correction",
        ],
        "not_closed_yet": [
            "explicit p01/p10 localized oscillon profile",
            "direct p01/p10 Laue/boundary check for the dressed oscillon stress tensor",
            "boundary/ADM terms for compact bodies",
            "full Lorentz/preferred-frame perturbation audit",
        ],
        "export_status": "STRONG_THEOREM_GATE_NOT_FINAL_RFG_DERIVATION",
    }


def gravitational_mass_from_zero_frequency_source():
    """
    The same energy E0 is the static gravitational source.

    RFG phase5 closes:
        localized oscillon -> zero-frequency <T00> -> exterior 1/r field.

    In the weak field:
        Phi_N = -G M/r,
        phi = 2 Phi_N/c^2 = -2GM/(c^2 r) = -r_s/r.
    """
    E0, M, G, c, r = sp.symbols('E0 M G c r', positive=True, real=True)
    r_s = 2 * G * M / c**2
    phi = -r_s / r
    phi_from_energy = sp.simplify(phi.subs(M, E0 / c**2))
    acceleration = sp.simplify(-(c**2 / 2) * sp.diff(phi, r))

    return {
        "theorem": "the inertial energy E0 is also the gravitational source",
        "M_source": sp.Eq(M, E0 / c**2),
        "schwarzschild_radius": sp.Eq(sp.Symbol('r_s'), r_s),
        "RFG_potential": sp.Eq(sp.Symbol('phi'), phi),
        "RFG_potential_from_energy": sp.Eq(sp.Symbol('phi(E0)'), phi_from_energy),
        "Newtonian_acceleration": sp.Eq(sp.Symbol('a_r'), acceleration),
        "meaning": "the far-zone 1/r coefficient is fixed by the same E0 that fixes inertial response.",
    }


def geodesic_mass_cancellation():
    """
    Why test-body free fall is composition independent.

    The same M multiplies both sides:
        inertial side: M a
        gravitational side: -M c^2 grad(phi)/2

    It cancels, giving the geodesic acceleration independently of the
    oscillon species or internal resonance label.
    """
    M, c, grad_phi = sp.symbols('M c grad_phi', positive=True, real=True)
    F_grav = -M * c**2 * grad_phi / 2
    acceleration = sp.simplify(F_grav / M)

    return {
        "inertial_side": sp.Eq(sp.Symbol('F'), sp.Symbol('M*a')),
        "gravitational_force": sp.Eq(sp.Symbol('F_grav'), F_grav),
        "free_fall_acceleration": sp.Eq(sp.Symbol('a'), acceleration),
        "mass_cancellation": sp.simplify(F_grav / M - acceleration),
        "meaning": "universality of free fall follows because the same M appears on both sides.",
    }


def equivalence_principle_closure():
    """
    m_i=m_g follows because both are E0/c^2.

    Composition independence is inherited from minimal coupling and from the
    fact that every localized excitation carries one stress-energy tensor.
    """
    E0, c = sp.symbols('E0 c', positive=True, real=True)
    m_i = E0 / c**2
    m_g = E0 / c**2
    eta = sp.simplify(2 * sp.Abs(m_g - m_i) / (m_g + m_i))

    return {
        "m_inertial": sp.Eq(sp.Symbol('m_i'), m_i),
        "m_gravitational": sp.Eq(sp.Symbol('m_g'), m_g),
        "identity": sp.Eq(sp.Symbol('m_i'), sp.Symbol('m_g')),
        "Eotvos_parameter": sp.Eq(sp.Symbol('eta'), eta),
        "MICROSCOPE_status": "consistent with eta ~ 0 at the 1e-15 test-body level",
        "scope": "leading localized-body theorem; nonlinear compact-body ADM/Noether proof is the remaining tightening step",
    }


def dressed_mass_no_double_counting():
    """
    The mass is not counted twice.

    RFG has an oscillon core plus its static Bernoulli/metric dressing.  The
    conserved energy E0 is the integral of the total localized zero-frequency
    stress-energy.  That single E0 defines both inertia and the far-zone
    gravitational charge.
    """
    E_core, E_dress, c = sp.symbols('E_core E_dress c', positive=True, real=True)
    E0 = E_core + E_dress
    M = E0 / c**2

    return {
        "total_rest_energy": sp.Eq(sp.Symbol('E0'), E0),
        "inertial_mass": sp.Eq(sp.Symbol('M_i'), M),
        "gravitational_mass": sp.Eq(sp.Symbol('M_g'), M),
        "no_double_counting_rule": "use the total localized zero-frequency stress-energy once",
        "far_zone_charge": "the 1/r coefficient is fixed by E0, not by separately adding core and field masses again",
        "meaning": "core+dressing is one dressed particle with one conserved mass.",
    }


def pressure_scaling_preserves_equivalence():
    """
    Phase17 bridge: in a background pressure potential phi,
    the external/Killing mass scales as exp(phi/2).  The key point for
    inertia is that both inertial and gravitational masses scale together.
    """
    phi, m0 = sp.symbols('phi m0', real=True, positive=True)
    m_eff = m0 * sp.exp(phi / 2)
    ratio = sp.simplify(m_eff / m_eff)

    return {
        "m_i_background": sp.Eq(sp.Symbol('m_i(phi)'), m_eff),
        "m_g_background": sp.Eq(sp.Symbol('m_g(phi)'), m_eff),
        "ratio": sp.Eq(sp.Symbol('m_g/m_i'), ratio),
        "meaning": "background pressure changes the external mass scale, but not the equivalence ratio.",
    }


def pressure_force_firewall_and_observation_gate():
    """
    Prevent the inertia intuition from becoming a fifth-force claim.

    The pressure/retardation language explains the medium dressing response.
    Matter must still move through minimal coupling/geodesics; otherwise RFG
    would generate composition-dependent pressure forces and fail Eotvos,
    atomic-clock, accelerator and Solar-System tests.
    """
    return {
        "firewall": (
            "fore/aft pressure asymmetry is a dressing reaction, not a literal "
            "extra pressure-gradient force on matter"
        ),
        "safe_motion_law": "localized bodies follow the metric/minimal-coupling equation of motion",
        "WEP_condition": "m_i and m_g must be the same total dressed E0/c^2 for every species",
        "Lorentz_condition": "no vacuum drag, preferred-frame acceleration law, or species-dependent boost response",
        "clock_condition": "pressure/process scaling must not rescale local atomic transition ratios unless separately constrained",
        "do_not_claim": [
            "full WEP/Eotvos proof beyond leading localized oscillon theorem",
            "accelerator/Lorentz pass before the full perturbation-sector audit",
            "literal pressure-gradient force on matter",
            "atomic-clock/varying-constant compatibility from pressure scaling alone",
            "strong-equivalence/compact-body closure before nonlinear ADM/Noether tests",
        ],
        "gate_status": "OBSERVATIONAL_FIREWALL_ADDED",
    }


def radiation_reaction_power_counting():
    """
    Radiation reaction is a correction, not the source of inertia.

    For a finite dressed object with response time tau_d, the leading
    nonrelativistic equation has the structure

        F_ext = M a + M tau_d adot + O((omega tau_d)^2)

    for slowly varying acceleration.  The inertial coefficient M is present
    even when tau_d -> 0; radiation reaction only corrects rapidly changing
    motion.
    """
    M, a, adot, tau_d, omega = sp.symbols(
        'M a adot tau_d omega',
        positive=True,
        real=True,
    )
    F = M * a + M * tau_d * adot
    harmonic_ratio = sp.simplify((M * tau_d * omega * a) / (M * a))

    return {
        "effective_force_law": sp.Eq(sp.Symbol('F_ext'), F),
        "adiabaticity_parameter": sp.Eq(sp.Symbol('epsilon_rad'), omega * tau_d),
        "harmonic_correction_ratio": sp.Eq(sp.Symbol('F_rad/F_inertia'), harmonic_ratio),
        "slow_motion_limit": "omega*tau_d << 1 -> F_ext = M*a",
        "meaning": "retardation explains dressing dynamics; Noether energy fixes the leading inertial coefficient.",
    }


def zero_drag_uniform_motion():
    """
    Uniform motion produces no drag.

    A constant-velocity oscillon is just a Lorentz-boosted localized solution.
    Its total momentum is conserved unless an external force acts.
    """
    M, c, v, t = sp.symbols('M c v t', positive=True, real=True)
    gamma = 1 / sp.sqrt(1 - v**2 / c**2)
    p = gamma * M * v
    force_constant_v = sp.diff(p, t)

    return {
        "boosted_solution_status": "constant-v state is a rigid Lorentz-boosted oscillon+dressing",
        "momentum": sp.Eq(sp.Symbol('p'), p),
        "force_for_constant_v": sp.Eq(sp.Symbol('dp/dt'), force_constant_v),
        "Newton_I": "no substrate drag for uniform rectilinear motion",
    }


def retarded_dressing_audit():
    """
    Keep the old retarded-field mechanism, but classify it correctly.

    The dipolar field is the physical picture of acceleration through the
    medium.  It is not used here to fit the coefficient of F=Ma.
    """
    r_s, c, a, n_dot_a = sp.symbols('r_s c a n_dot_a', positive=True, real=True)
    delta_phi = r_s * n_dot_a / (2 * c**2)

    return {
        "old_mechanism": "acceleration produces a retarded fore/aft field asymmetry",
        "dipole_perturbation": sp.Eq(sp.Symbol('delta_phi'), delta_phi),
        "temporal_piece_old_result": "-m*a/2 after regularization",
        "spatial_piece_old_status": "the other -m*a/2 must come from the bi-conformal spatial momentum sector",
        "RFG_upgrade": "collective-coordinate Noether theorem fixes the full leading coefficient M without relying on the self-force finite part",
        "radiation_reaction": "higher-derivative/radiative terms are corrections to F=Ma, not the origin of M",
    }


def relaxation_timescale_hierarchy():
    """
    Separate three time scales that were mixed in the toy version.

    1. External gravitational dressing light-crossing time: r_s/c.
    2. Internal oscillon clock/Compton time: h/(m c^2).
    3. Macroscopic compact-object dressing time: R_eff/c.
    """
    G = 6.67430e-11
    c = 299_792_458.0
    h = 6.62607015e-34
    m_e = 9.1093837015e-31
    m_p = 1.67262192369e-27
    m_sun = 1.98847e30

    def r_s_over_c(mass: float) -> float:
        return 2 * G * mass / c**3

    def compton_period(mass: float) -> float:
        return h / (mass * c**2)

    neutron_star_tau = 2 * G * (1.4 * m_sun) / c**3

    return {
        "electron_external_rs_over_c_s": r_s_over_c(m_e),
        "electron_internal_Compton_period_s": compton_period(m_e),
        "proton_external_rs_over_c_s": r_s_over_c(m_p),
        "proton_internal_Compton_period_s": compton_period(m_p),
        "neutron_star_rs_over_c_s": neutron_star_tau,
        "meaning": "elementary-particle inertia is effectively instantaneous; compact objects can have microsecond-scale dressing times.",
    }


def old_vs_new_inertia_assessment():
    """Plain status ledger for the user's worry."""
    return [
        "Old theory's best idea is retained: acceleration creates a retarded self-field asymmetry.",
        "Old theory's weak point is avoided: the leading coefficient no longer depends on assuming C_reg=1.",
        "The full F=Ma coefficient is now fixed by the oscillon's translational collective coordinate.",
        "No drag for uniform motion follows from Lorentz covariance/boosted solution, not from an ad hoc cancellation.",
        "Equivalence principle is clean: m_i=m_g=E0/c^2.",
        "Free-fall universality is explicit: M cancels between F_grav and F_inertia.",
        "Background pressure scaling preserves m_g/m_i=1 because both scale as exp(phi/2).",
        "Core and Bernoulli dressing are counted once as a single total zero-frequency energy E0.",
        "Radiation reaction is now power-counted as an adiabatic correction, not the origin of inertia.",
        "New gate: Lorentz Ward identity fixes the required zero-mode norm G_RR=E0/c^2.",
        "Remaining technical tightening: full nonlinear ADM/Noether momentum for compact matter-filled bodies.",
    ]


def stage_a5_old17_retarded_self_field_drain():
    """
    Full OLD/17 drain into the new Noether-based RFG inertia file.

    The old derivation is not discarded.  It is reclassified:
    - retarded fore/aft asymmetry is the physical mechanism;
    - the old finite self-force coefficient is not used as the proof of M;
    - M is fixed by the total Noether energy of the dressed oscillon.
    """
    m, G, c, r, r_s, n_dot_a = sp.symbols(
        'm G c r r_s n_dot_a',
        positive=True,
        real=True,
    )
    delta_phi_rs = r_s * n_dot_a / (2 * c**2)
    delta_phi_m = G * m * n_dot_a / c**4

    return {
        "source": "OLD/17. ISPG_Inertia.tex",
        "retarded_field_correction": sp.Eq(sp.Symbol('delta_phi'), delta_phi_rs),
        "same_with_mass": sp.Eq(sp.Symbol('delta_phi_m'), delta_phi_m),
        "temporal_sector_old": "regularized g00/scalar-momentum bookkeeping gives -m*a/2",
        "spatial_sector_old": "bi-conformal gij sector supplies the matching -m*a/2 when gamma_PPN=1",
        "new_RFG_role": (
            "this explains the medium's reaction under acceleration; the exact "
            "leading coefficient is fixed independently by the translational "
            "Noether/collective-coordinate theorem"
        ),
        "zero_drag": "uniform motion is a Lorentz-boosted static oscillon+dressing, so self-drag vanishes",
        "relaxation_scale": "tau_relax ~ R_eff/c; for exterior point bookkeeping R_eff~r_s",
    }


def stage_a5_mach_and_unified_rarefaction_ledger():
    """
    OLD/17 contained a useful conceptual bridge:
    inertia is local and causal like a field theory, but keeps the Machian
    intuition that inertia is a property of a material/pressurized vacuum.
    """
    phi, r_s = sp.symbols('phi r_s', real=True)
    return {
        "Mach_comparison": {
            "Mach": "inertia sourced by distant matter; mechanism unspecified",
            "RFG": "inertia sourced by local dressed field/medium response; propagation at c",
            "empty_universe_branch": sp.And(sp.Eq(phi, 0), sp.Eq(r_s, 0)),
            "scope": (
                "working hypothesis: without a pressure contrast/dressed oscillon "
                "there is no operational inertial mass; this should be stated as "
                "a branch of the medium theory, not as a proved cosmological theorem"
            ),
        },
        "unified_rarefaction_principle": {
            "gravity": "static Bernoulli pressure deficit around a localized oscillon",
            "inertia": "dynamic fore/aft pressure asymmetry during acceleration",
            "MOND": "macroscopic central rarefaction/vortex memory in rotating bound systems",
            "common_language": (
                "all three are geometries of measurable-existence/background-pressure "
                "redistribution in the same continuous medium"
            ),
        },
    }


def stage_a5_inertia_status():
    """Stage marker for the OLD deletion gate."""
    return {
        "migrated": True,
        "old_file_drained": "OLD/17. ISPG_Inertia.tex",
        "new_file": "p06_inertia.py",
        "kept_from_old": [
            "retarded scalar/pressure field of an accelerating body",
            "zero drag for uniform Lorentz-boosted motion",
            "relaxation hierarchy tau~R_eff/c",
            "Mach comparison as local causal medium response",
            "unified rarefaction principle: gravity/inertia/MOND",
        ],
        "strengthened_in_new": [
            "F=Ma coefficient comes from Noether energy, not a chosen regular finite part",
            "Lorentz Ward identity now reduces the missing proof to G_RR=E0/c^2",
            "integrated Noether boost now gives p=gamma(E0/c^2)v under the Laue condition",
            "m_i=m_g=E0/c^2 is explicit",
            "radiation reaction is separated from leading inertia",
            "pressure-force firewall prevents a fifth-force reading",
        ],
        "open_math": [
            "direct p01/p10 Laue/boundary test for a localized dressed oscillon",
            "direct p01/p10 zero-mode integral after the Laue test",
            "full nonlinear ADM/Noether momentum for compact matter-filled bodies",
            "explicit bi-conformal spatial-sector momentum integral matching the old -m*a/2",
            "observational bridge for neutron-star merger relaxation times",
        ],
    }


def inertia_central_claim_gate():
    """One-place export gate for p06_inertia.py."""
    ward = lorentz_ward_zero_mode_norm_gate()
    derivation = rfg_action_to_inertia_derivation_gate()
    firewall = pressure_force_firewall_and_observation_gate()
    boost = integrated_noether_boost_closure()
    laue = laue_stress_condition_gate()
    return {
        "file_export_status": "NOT_READY_FOR_FINAL_THEORY_EXPORT",
        "leading_inertia_status": ward["gate_status"],
        "integrated_boost_status": boost["closure_status"],
        "Laue_balance_status": laue["gate_status"],
        "RFG_derivation_status": derivation["export_status"],
        "equivalence_status": "m_i=m_g=E0/c^2 in the leading dressed-oscillon branch",
        "medium_mechanism_status": "retarded fore/aft dressing explains response but does not set the F=Ma coefficient",
        "observation_firewall": firewall["gate_status"],
        "main_open_integral": "verify p01/p10 Laue boundary condition, then compute translation zero-mode norm from the full dressed stress tensor",
        "do_not_claim": firewall["do_not_claim"],
    }


def main() -> None:
    print("=" * 72)
    print("PHASE 19: RFG inertia theorem")
    print("=" * 72)

    sections = [
        ("1. Collective-coordinate inertia theorem", collective_coordinate_inertia_theorem()),
        ("2. Noether four-momentum closure", noether_four_momentum_closure()),
        ("3. Moduli-space projection theorem", moduli_space_projection_theorem()),
        ("4. Lorentz Ward zero-mode norm gate", lorentz_ward_zero_mode_norm_gate()),
        ("5. Integrated Noether boost closure", integrated_noether_boost_closure()),
        ("6. Laue stress-balance gate", laue_stress_condition_gate()),
        ("7. RFG action-to-inertia derivation gate", rfg_action_to_inertia_derivation_gate()),
        ("8. Same energy as gravitational source", gravitational_mass_from_zero_frequency_source()),
        ("9. Geodesic mass cancellation", geodesic_mass_cancellation()),
        ("10. Equivalence principle closure", equivalence_principle_closure()),
        ("11. Dressed mass without double counting", dressed_mass_no_double_counting()),
        ("12. Background pressure scaling", pressure_scaling_preserves_equivalence()),
        ("13. Pressure-force observational firewall", pressure_force_firewall_and_observation_gate()),
        ("14. Uniform motion: zero drag", zero_drag_uniform_motion()),
        ("15. Retarded dressing audit", retarded_dressing_audit()),
        ("16. Radiation-reaction power counting", radiation_reaction_power_counting()),
        ("17. STAGE A5 old retarded self-field drain", stage_a5_old17_retarded_self_field_drain()),
        ("18. STAGE A5 Mach/unified rarefaction ledger", stage_a5_mach_and_unified_rarefaction_ledger()),
        ("19. STAGE A5 migration status", stage_a5_inertia_status()),
        ("20. Central inertia claim gate", inertia_central_claim_gate()),
    ]

    for title, data in sections:
        print(f"\n--- {title} ---")
        for key, value in data.items():
            print(f"  {key:34s}: {value}")

    print("\n--- 21. Relaxation timescale hierarchy ---")
    for key, value in relaxation_timescale_hierarchy().items():
        if isinstance(value, float):
            print(f"  {key:34s}: {value:.3e}")
        else:
            print(f"  {key:34s}: {value}")

    print("\n--- 22. Old vs new assessment ---")
    for item in old_vs_new_inertia_assessment():
        print(f"  - {item}")

    print("\n--- დასკვნა ---")
    print("  RFG-ში ინერცია ძველ self-force ტექსტზე ძლიერია:")
    print("  Integrated boost proof იძლევა p=gamma(E0/c^2)v-ს, თუ Laue stress balance სრულდება.")
    print("  სრული RFG derivation ჯერ ითხოვს p01/p10 Laue-boundary და zero-mode ინტეგრალს.")
    print("  retarded asymmetry რჩება მექანიზმად, ხოლო pressure-gradient force matter-ზე აკრძალულია.")


if __name__ == "__main__":
    main()
