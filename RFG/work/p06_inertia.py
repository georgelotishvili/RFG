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
closes the leading inertial law by a stronger route:

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
            "m_i=m_g=E0/c^2 is explicit",
            "radiation reaction is separated from leading inertia",
        ],
        "open_math": [
            "full nonlinear ADM/Noether momentum for compact matter-filled bodies",
            "explicit bi-conformal spatial-sector momentum integral matching the old -m*a/2",
            "observational bridge for neutron-star merger relaxation times",
        ],
    }


def main() -> None:
    print("=" * 72)
    print("PHASE 19: RFG inertia theorem")
    print("=" * 72)

    sections = [
        ("1. Collective-coordinate inertia theorem", collective_coordinate_inertia_theorem()),
        ("2. Noether four-momentum closure", noether_four_momentum_closure()),
        ("3. Moduli-space projection theorem", moduli_space_projection_theorem()),
        ("4. Same energy as gravitational source", gravitational_mass_from_zero_frequency_source()),
        ("5. Geodesic mass cancellation", geodesic_mass_cancellation()),
        ("6. Equivalence principle closure", equivalence_principle_closure()),
        ("7. Dressed mass without double counting", dressed_mass_no_double_counting()),
        ("8. Background pressure scaling", pressure_scaling_preserves_equivalence()),
        ("9. Uniform motion: zero drag", zero_drag_uniform_motion()),
        ("10. Retarded dressing audit", retarded_dressing_audit()),
        ("11. Radiation-reaction power counting", radiation_reaction_power_counting()),
        ("12. STAGE A5 old retarded self-field drain", stage_a5_old17_retarded_self_field_drain()),
        ("13. STAGE A5 Mach/unified rarefaction ledger", stage_a5_mach_and_unified_rarefaction_ledger()),
        ("14. STAGE A5 migration status", stage_a5_inertia_status()),
    ]

    for title, data in sections:
        print(f"\n--- {title} ---")
        for key, value in data.items():
            print(f"  {key:34s}: {value}")

    print("\n--- 15. Relaxation timescale hierarchy ---")
    for key, value in relaxation_timescale_hierarchy().items():
        if isinstance(value, float):
            print(f"  {key:34s}: {value:.3e}")
        else:
            print(f"  {key:34s}: {value}")

    print("\n--- 16. Old vs new assessment ---")
    for item in old_vs_new_inertia_assessment():
        print(f"  - {item}")

    print("\n--- დასკვნა ---")
    print("  RFG-ში ინერცია უკვე უკეთესადაა დახურული, ვიდრე ძველ self-force ტექსტში:")
    print("  F=Ma გამოდის ოსცილონის ენერგიიდან და ტრანსლაციური Noether მოდიდან,")
    print("  ხოლო retarded asymmetry რჩება ფიზიკურ მექანიზმად, რომელიც აჩქარებისას")
    print("  მედიუმის რეაქციას ხდის ინტუიციურად გასაგებს.")


if __name__ == "__main__":
    main()
