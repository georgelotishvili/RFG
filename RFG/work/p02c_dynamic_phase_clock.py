# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# Active coefficient scheme: Y-scheme; c_Y means c_Y^(Y), not X-scheme c_X.

"""
p02c: dynamic phase-clock branch ledger.

Scope:
- derive the normalized and canonical Phi-current on FLRW;
- record the dynamic Phi_dot(a) branch as a candidate, not a metric claim;
- check the late nonzero zero-current algebraic w=-1 result;
- show how dynamic Phi_dot(a) reshuffles early scaling terms;
- block process-time identification until channel separation, perturbations, and fit.

This file is not the primary FLRW/CMB metric branch. Do not import its clock
factor into H(z), CMB, BBN, BAO, or SN calculations without a separate fit.
"""

import sympy as sp


def phase_current_self_check():
    """
    Check the Phi-current directly from the FLRW phase/mixing Lagrangian.

    The canonical Noether density has an overall factor 2 relative to the
    normalized current used for the algebraic branch. The factor is absorbable
    into the integration constant, but the distinction is recorded explicitly.
    """
    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    u = sp.Function("u")(t)  # u(t) = Phi_dot(t)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    I1 = 3 / a**2
    L_phase_mix = c_Y * u**2 + c_Y2 * u**4 + c_YI1 * u**2 * I1

    canonical_density = sp.simplify(sp.diff(L_phase_mix, u))
    normalized_density = sp.simplify(canonical_density / 2)
    expected_normalized = sp.simplify(
        u * (c_Y + 2 * c_Y2 * u**2 + 3 * c_YI1 / a**2)
    )

    canonical_eom = sp.simplify(sp.diff(a**3 * canonical_density, t) / a**3)
    normalized_eom = sp.simplify(2 * sp.diff(a**3 * expected_normalized, t) / a**3)

    return {
        "status": "PASS" if sp.simplify(normalized_density - expected_normalized) == 0
        and sp.simplify(canonical_eom - normalized_eom) == 0
        else "CHECK",
        "canonical_density": canonical_density,
        "normalized_density": normalized_density,
        "normalized_density_residual": sp.simplify(
            normalized_density - expected_normalized
        ),
        "canonical_EOM": canonical_eom,
        "canonical_minus_normalized_EOM_residual": sp.simplify(
            canonical_eom - normalized_eom
        ),
        "normalization_note": (
            "canonical current = 2 * normalized current; Q_Phi changes by the "
            "same constant factor"
        ),
    }


def dynamic_phase_clock_branch():
    """Algebraic candidate branch for u(a)=Phi_dot(a)."""
    a = sp.Symbol("a", positive=True)
    u = sp.Symbol("u", real=True)
    Q_norm, Q_can = sp.symbols("Q_norm Q_can", real=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    normalized_density = sp.simplify(
        u * (c_Y + 2 * c_Y2 * u**2 + 3 * c_YI1 / a**2)
    )
    canonical_density = sp.simplify(2 * normalized_density)
    normalized_charge = sp.simplify(a**3 * normalized_density)
    canonical_charge = sp.simplify(a**3 * canonical_density)

    dynamic_equation_normalized = sp.Eq(
        sp.simplify(2 * c_Y2 * u**3 + (c_Y + 3 * c_YI1 / a**2) * u),
        Q_norm / a**3,
    )
    zero_current_nonzero_u2 = sp.simplify(
        -(c_Y + 3 * c_YI1 / a**2) / (2 * c_Y2)
    )

    return {
        "status": "ALGEBRAIC_DYNAMIC_BRANCH_CANDIDATE",
        "normalized_current": sp.Eq(normalized_charge, Q_norm),
        "canonical_current": sp.Eq(canonical_charge, Q_can),
        "dynamic_equation_normalized": dynamic_equation_normalized,
        "zero_current_roots": [
            sp.Eq(u, 0),
            sp.Eq(u**2, zero_current_nonzero_u2),
        ],
        "nonzero_zero_current_late_u2": sp.Eq(
            sp.Symbol("u_late", real=True) ** 2,
            sp.simplify(sp.limit(zero_current_nonzero_u2, a, sp.oo)),
        ),
        "status_note": (
            "candidate branch only; not a process-time identity and not a "
            "dark-energy solution"
        ),
    }


def late_zero_current_candidate():
    """Late nonzero zero-current algebra, checked by direct substitution."""
    u2 = sp.Symbol("u2", real=True)
    c_Y, c_Y2 = sp.symbols("c_Y c_Y2", real=True)

    late_u2 = sp.simplify(-c_Y / (2 * c_Y2))
    rho_late = sp.simplify(c_Y * u2 + 3 * c_Y2 * u2**2)
    p_late = sp.simplify(c_Y * u2 + c_Y2 * u2**2)

    rho_sub = sp.simplify(rho_late.subs(u2, late_u2))
    p_sub = sp.simplify(p_late.subs(u2, late_u2))
    w_sub = sp.simplify(p_sub / rho_sub)

    return {
        "status": "ALGEBRAIC_LATE_BRANCH_CANDIDATE",
        "branch": "Q_norm=0, u != 0, a -> infinity",
        "late_u2": sp.Eq(u2, late_u2),
        "rho_after_substitution": rho_sub,
        "p_after_substitution": p_sub,
        "w_after_substitution": w_sub,
        "theorem_pass_background_only": sp.simplify(w_sub + 1) == 0,
        "late_kinetic_prefactor": sp.simplify(4 * c_Y2 * late_u2),
        "viable_sign_window_hint": [
            sp.Gt(c_Y2, 0),
            sp.Lt(c_Y, 0),
        ],
        "sign_note": (
            "c_Y < 0 here is not a pure c_Y no-ghost rule; the effective "
            "background kinetic prefactor includes the Y^2 branch"
        ),
        "blocked_until": [
            "u^2 positivity over the epoch range",
            "no-ghost, gradient, sound-speed, and eigenmode stability",
            "full gravity closure",
            "numerical BBN/CMB/Planck/BAO/SN fit",
        ],
    }


def early_scaling_after_zero_current():
    """
    Show that substituting the dynamic branch reshuffles early powers.

    This is the A-level council warning: the unsubstituted p02 early scaling is
    not enough once u(a) is solved dynamically.
    """
    a = sp.Symbol("a", positive=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)
    c_I1, c_I1sq, c_I2, c_I3 = sp.symbols(
        "c_I1 c_I1sq c_I2 c_I3",
        real=True,
    )

    A = c_Y + 3 * c_YI1 / a**2
    rho_phase_mix_zero_current = sp.simplify(A**2 / (4 * c_Y2))
    rho_full_zero_current = sp.expand(
        -3 * c_I1 / a**2
        - (9 * c_I1sq + 3 * c_I2) / a**4
        - c_I3 / a**6
        + rho_phase_mix_zero_current
    )

    return {
        "status": "BLOCKED_UNTIL_NUMERICAL_BBN_CMB_PLANCK_BAO_BOUNDS",
        "rho_phase_mix_zero_current": rho_phase_mix_zero_current,
        "rho_full_zero_current_expanded": rho_full_zero_current,
        "added_effective_a_minus_2": sp.simplify(
            3 * c_Y * c_YI1 / (2 * c_Y2 * a**2)
        ),
        "added_effective_a_minus_4": sp.simplify(
            9 * c_YI1**2 / (4 * c_Y2 * a**4)
        ),
        "interpretation": (
            "dynamic u(a) can move phase-mixing energy into curvature-like and "
            "radiation-like effective pieces"
        ),
        "do_not_claim": (
            "do not use the p02 unsubstituted a^-2/a^-4/a^-6 ledger as a BBN/CMB "
            "pass after imposing this branch"
        ),
    }


def process_time_match_gate():
    """Comparison to p02b process-time, with metric-branch blocking."""
    a = sp.Symbol("a", positive=True)
    T0 = sp.Symbol("T0", positive=True)
    t_age = sp.Function("t_age")
    C_proc = sp.Function("C_proc")
    Q_norm = sp.Symbol("Q_norm", real=True)
    c_Y, c_Y2, c_YI1 = sp.symbols("c_Y c_Y2 c_YI1", real=True)

    return {
        "status": "BLOCKED_UNTIL_CHANNEL_SEPARATION_AND_FIT",
        "p02b_target": sp.Eq(C_proc(a), T0 / t_age(a)),
        "general_match_condition": sp.Eq(
            2 * c_Y2 * C_proc(a) ** 3
            + (c_Y + 3 * c_YI1 / a**2) * C_proc(a),
            Q_norm / a**3,
        ),
        "zero_current_match_condition": sp.Eq(
            C_proc(a) ** 2,
            -(c_Y + 3 * c_YI1 / a**2) / (2 * c_Y2),
        ),
        "main_gap": (
            "constant p02 coefficients have not been shown to reproduce p02b "
            "self-similar C_proc(a)=T0/t_age(a) over the observational range"
        ),
        "metric_branch_guardrail": (
            "this match must not be inserted into H(z), CMB, BBN, BAO, or SN "
            "without a separate channel-separation proof and numerical fit"
        ),
    }


def module_status():
    """Council-facing status for the dynamic phase-clock file."""
    return {
        "scope": "dynamic phase-clock candidate only; not primary FLRW metric branch",
        "current_self_check": phase_current_self_check(),
        "dynamic_branch": dynamic_phase_clock_branch(),
        "late_zero_current_candidate": late_zero_current_candidate(),
        "early_scaling_after_zero_current": early_scaling_after_zero_current(),
        "process_time_match_gate": process_time_match_gate(),
        "export_status": "NOT_READY_FOR_RFG_THEORY_EXPORT",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("p02c: dynamic phase-clock branch ledger")
    print("=" * 72)

    current = phase_current_self_check()
    print("\n1. Current self-check")
    print("status:", current["status"])
    print("canonical density:", current["canonical_density"])
    print("normalized density:", current["normalized_density"])
    print("density residual:", current["normalized_density_residual"])
    print("EOM residual:", current["canonical_minus_normalized_EOM_residual"])
    print("note:", current["normalization_note"])

    branch = dynamic_phase_clock_branch()
    print("\n2. Dynamic branch")
    print("status:", branch["status"])
    print("normalized current:", branch["normalized_current"])
    print("canonical current:", branch["canonical_current"])
    print("dynamic equation:", branch["dynamic_equation_normalized"])
    print("zero-current roots:", branch["zero_current_roots"])

    late = late_zero_current_candidate()
    print("\n3. Late zero-current candidate")
    print("status:", late["status"])
    print("rho:", late["rho_after_substitution"])
    print("p:", late["p_after_substitution"])
    print("w:", late["w_after_substitution"])
    print("blocked until:", late["blocked_until"])

    early = early_scaling_after_zero_current()
    print("\n4. Early-scaling warning")
    print("status:", early["status"])
    print("rho expanded:", early["rho_full_zero_current_expanded"])
    print("added a^-2:", early["added_effective_a_minus_2"])
    print("added a^-4:", early["added_effective_a_minus_4"])

    match = process_time_match_gate()
    print("\n5. Process-time match gate")
    print("status:", match["status"])
    print("target:", match["p02b_target"])
    print("general match:", match["general_match_condition"])
    print("guardrail:", match["metric_branch_guardrail"])
