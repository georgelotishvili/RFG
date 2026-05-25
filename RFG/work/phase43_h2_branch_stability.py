# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 43: Normal-form stability of the h=2 charged-lepton branch

Status:
    Local normal-form stability conditions. This is not yet the full 3D PDE
    oscillon stability proof, but it converts one open theorem item from
    phase42 into explicit algebraic inequalities.

Variables:
    rho   : amplitude of the triaxial strain doublet Q
    beta  : C3 strain phase, with det(Q) = rho^3 cos(3 beta)/4
    theta : reduced framed holonomy, selected by theta = 2/9

Effective potential:
    V = a2 rho^2 + a4 rho^4
        - lambda3 rho^3 cos(3 beta)/4
        + kappa [1 - cos(pi(9 theta - 2))]
        + g rho^2 (theta - 2/9)^2

Meaning:
    - lambda3 term comes from I3=det(B), see phase41.
    - kappa term is the reduced framed holonomy lock, see phase38.
    - g term is the lowest phase-elastic backreaction allowed near theta=2/9.

Main result:
    At beta=0 and theta=2/9, the branch is locally stable if

        lambda3 > 0,
        kappa > 0,
        8 a4 rho0 - 3 lambda3/4 > 0,
        81 pi^2 kappa + 2 g rho0^2 > 0.

    rho0 is the nonzero stationary amplitude solving

        2 a2 + 4 a4 rho0^2 - 3 lambda3 rho0 / 4 = 0.
"""

import sympy as sp


def normal_form_potential():
    rho, beta, theta = sp.symbols("rho beta theta", real=True)
    a2, a4, lambda3, kappa, g = sp.symbols(
        "a2 a4 lambda3 kappa g", real=True
    )
    theta0 = sp.Rational(2, 9)
    V = (
        a2 * rho**2
        + a4 * rho**4
        - lambda3 * rho**3 * sp.cos(3 * beta) / 4
        + kappa * (1 - sp.cos(sp.pi * (9 * theta - 2)))
        + g * rho**2 * (theta - theta0) ** 2
    )
    return rho, beta, theta, a2, a4, lambda3, kappa, g, theta0, V


def stationary_conditions():
    rho, beta, theta, a2, a4, lambda3, kappa, g, theta0, V = normal_form_potential()
    subs_branch = {beta: 0, theta: theta0}
    d_rho = sp.factor(sp.diff(V, rho).subs(subs_branch))
    d_beta = sp.factor(sp.diff(V, beta).subs(subs_branch))
    d_theta = sp.factor(sp.diff(V, theta).subs(subs_branch))
    nonzero_amplitude_eq = sp.factor(d_rho / rho)
    rho_solutions = sp.solve(nonzero_amplitude_eq, rho)
    return {
        "dV_drho_branch": d_rho,
        "dV_dbeta_branch": d_beta,
        "dV_dtheta_branch": d_theta,
        "nonzero_amplitude_eq": nonzero_amplitude_eq,
        "rho_solutions": rho_solutions,
    }


def hessian_on_branch():
    rho, beta, theta, a2, a4, lambda3, kappa, g, theta0, V = normal_form_potential()
    variables = [rho, beta, theta]
    H = sp.Matrix([[sp.diff(V, x, y) for y in variables] for x in variables])
    H_branch = sp.simplify(H.subs({beta: 0, theta: theta0}))

    # Use the nonzero stationarity equation to eliminate a2.
    a2_sub = sp.solve(
        sp.Eq(2 * a2 + 4 * a4 * rho**2 - 3 * lambda3 * rho / 4, 0),
        a2,
    )[0]
    H_stationary = sp.simplify(H_branch.subs(a2, a2_sub))

    return {
        "variables": ["rho", "beta", "theta"],
        "H_branch": H_branch,
        "H_stationary": H_stationary,
        "diagonal_stationary": [sp.factor(H_stationary[i, i]) for i in range(3)],
        "offdiag_stationary": [
            sp.factor(H_stationary[0, 1]),
            sp.factor(H_stationary[0, 2]),
            sp.factor(H_stationary[1, 2]),
        ],
    }


def stability_conditions():
    return {
        "beta_phase_stability": "lambda3 > 0 and rho0 > 0",
        "theta_holonomy_stability": "81*pi^2*kappa + 2*g*rho0^2 > 0",
        "amplitude_stability": "rho0*(8*a4*rho0 - 3*lambda3/4) > 0",
        "large_root_rule": "for a4>0, choose the larger nonzero rho0 root",
        "pde_remaining": "full spatial oscillon stability still requires fluctuation operator spectrum",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 43: h=2 branch normal-form stability")
    print("=" * 72)

    print("\n1. Stationary conditions at beta=0, theta=2/9")
    for key, value in stationary_conditions().items():
        print(f"  {key:26s}: {value}")

    print("\n2. Hessian on the h=2 branch")
    hessian = hessian_on_branch()
    print(f"  variables: {hessian['variables']}")
    print(f"  H_branch:")
    print(hessian["H_branch"])
    print(f"  H_stationary:")
    print(hessian["H_stationary"])
    print(f"  diagonal_stationary: {hessian['diagonal_stationary']}")
    print(f"  offdiag_stationary : {hessian['offdiag_stationary']}")

    print("\n3. Stability conditions")
    for key, value in stability_conditions().items():
        print(f"  {key:26s}: {value}")
