# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 37: C3 Koide operator for charged-lepton frequencies

Status:
    Structural candidate. This file replaces the empirical question

        why N = 5, 72, 295?

    with a sharper operator question:

        why does the lepton oscillon have a C3-cyclic internal spectrum
        with topological phase theta = 2/9?

Core ansatz:
    nu_k = A * [1 + sqrt(2) cos(theta + 2*pi*k/3)],  k = 0,1,2.

For theta = 2/9, sorting the three frequencies and normalizing the
smallest to the electron gives

    1 : 14.37951 : 58.97010

and therefore, with m proportional to nu^2, the charged-lepton masses.

Important:
    phase38_z9_theta_holonomy.py adds a candidate derivation:
    theta = 2/9 from a Z9 reduced framed-holonomy closure
    (C3 axes x C3 phase/braid sectors, spinorial index h=2).
    phase39_action_symmetry_to_z9.py gives action-level support for
    the C3 x C3 closure lattice.
    phase40_projective_spinor_h2.py gives support for h=2 from the
    difference between projective/nematic and oriented framed closure.
    phase41_action_normal_form_theta.py shows that the C3 strain lock is
    already present in the RFG action through I3 = det(B), because
    E^3 + Ebar^3 = 27 det(Q).
    A full theorem still requires deriving that closure from the RFG action.
"""

import cmath
import math


LEPTON_MASSES_MEV = {
    "electron": 0.51099895000,
    "muon": 105.6583755,
    "tau": 1776.86,
}

THETA_TOPOLOGICAL = 2.0 / 9.0
KOIDE_TARGET = 2.0 / 3.0


def measured_frequency_ratios():
    """Observed sqrt(m/m_e) frequency ratios."""
    m_e = LEPTON_MASSES_MEV["electron"]
    return {
        name: math.sqrt(mass / m_e)
        for name, mass in LEPTON_MASSES_MEV.items()
    }


def c3_raw_frequencies(theta=THETA_TOPOLOGICAL):
    """
    Three eigenfrequencies of the C3-cyclic internal operator.

    The ordering is the natural cyclic ordering k=0,1,2. The physical
    charged leptons are assigned after sorting by frequency.
    """
    return [
        1.0 + math.sqrt(2.0) * math.cos(theta + 2.0 * math.pi * k / 3.0)
        for k in range(3)
    ]


def c3_frequency_ratios(theta=THETA_TOPOLOGICAL):
    """Sorted C3 frequencies normalized to the smallest frequency."""
    values = sorted(c3_raw_frequencies(theta))
    base = values[0]
    if base <= 0:
        raise ValueError("C3 ansatz produced a non-positive base frequency.")
    return {
        "electron": values[0] / base,
        "muon": values[1] / base,
        "tau": values[2] / base,
    }


def c3_mass_predictions(theta=THETA_TOPOLOGICAL):
    """Masses predicted from m ~ nu^2, anchored to the electron mass."""
    ratios = c3_frequency_ratios(theta)
    m_e = LEPTON_MASSES_MEV["electron"]
    return {
        name: m_e * ratio * ratio
        for name, ratio in ratios.items()
    }


def koide_ratio_from_frequencies(frequencies):
    """Koide ratio in frequency variables: sum(nu^2) / sum(nu)^2."""
    total = sum(frequencies)
    return sum(nu * nu for nu in frequencies) / (total * total)


def koide_identity(theta=THETA_TOPOLOGICAL):
    """
    The C3 ansatz gives Koide's 2/3 identity independently of theta,
    as long as the three frequencies are used as one C3 triplet.
    """
    frequencies = c3_raw_frequencies(theta)
    return koide_ratio_from_frequencies(frequencies)


def c3_circulant_operator(theta=THETA_TOPOLOGICAL, omega0=1.0):
    """
    Return the 3x3 complex circulant operator

        Omega = omega0 * [I + (e^(i theta) P + e^(-i theta) P^2)/sqrt(2)]

    where P cycles x -> y -> z -> x.
    """
    phase = cmath.exp(1j * theta) / math.sqrt(2.0)
    phase_conj = cmath.exp(-1j * theta) / math.sqrt(2.0)
    return [
        [omega0, omega0 * phase, omega0 * phase_conj],
        [omega0 * phase_conj, omega0, omega0 * phase],
        [omega0 * phase, omega0 * phase_conj, omega0],
    ]


def c3_operator_eigenvalues(theta=THETA_TOPOLOGICAL, omega0=1.0):
    """Analytic eigenvalues of the C3 circulant operator."""
    return [
        omega0 * value
        for value in c3_raw_frequencies(theta)
    ]


def axis_ratios_from_c3(theta=THETA_TOPOLOGICAL):
    """
    If an axis length scales as L ~ 1/nu, return the C3-predicted
    ellipsoid ratios normalized to the electron axis.
    """
    ratios = c3_frequency_ratios(theta)
    return {
        name: 1.0 / ratio
        for name, ratio in ratios.items()
    }


def fit_error(theta):
    """Squared log-error against observed charged-lepton frequency ratios."""
    predicted = c3_frequency_ratios(theta)
    observed = measured_frequency_ratios()
    return sum(
        (math.log(predicted[name]) - math.log(observed[name])) ** 2
        for name in ("electron", "muon", "tau")
    )


def best_fit_theta(grid_size=20000, refine_steps=80):
    """
    Find the best one-parameter theta fit near the theta = 2/9 branch.

    This is an audit: the theory target remains theta = 2/9, not the
    fitted value. The sorted C3 spectrum has equivalent permutation
    branches, so the audit intentionally stays in the branch containing
    theta = 2/9.
    """
    lo = 0.0
    hi = math.pi / 3.0
    best_theta = lo
    best_value = float("inf")

    for i in range(grid_size + 1):
        theta = lo + (hi - lo) * i / grid_size
        try:
            value = fit_error(theta)
        except ValueError:
            continue
        if value < best_value:
            best_theta = theta
            best_value = value

    width = (hi - lo) / grid_size
    a = max(lo, best_theta - 3.0 * width)
    b = min(hi, best_theta + 3.0 * width)

    for _ in range(refine_steps):
        c = a + (b - a) / 3.0
        d = b - (b - a) / 3.0
        try:
            fc = fit_error(c)
        except ValueError:
            fc = float("inf")
        try:
            fd = fit_error(d)
        except ValueError:
            fd = float("inf")
        if fc < fd:
            b = d
        else:
            a = c

    theta = 0.5 * (a + b)
    return {
        "theta_fit": theta,
        "theta_fit_over_pi": theta / math.pi,
        "theta_topological": THETA_TOPOLOGICAL,
        "theta_delta": theta - THETA_TOPOLOGICAL,
        "error_fit": fit_error(theta),
        "error_topological": fit_error(THETA_TOPOLOGICAL),
    }


def prediction_table(theta=THETA_TOPOLOGICAL):
    """Rows comparing observed and C3-predicted masses/frequencies."""
    observed_ratios = measured_frequency_ratios()
    predicted_ratios = c3_frequency_ratios(theta)
    predicted_masses = c3_mass_predictions(theta)
    rows = []
    for name in ("electron", "muon", "tau"):
        observed_mass = LEPTON_MASSES_MEV[name]
        predicted_mass = predicted_masses[name]
        rows.append(
            {
                "particle": name,
                "observed_freq_ratio": observed_ratios[name],
                "predicted_freq_ratio": predicted_ratios[name],
                "observed_mass_MeV": observed_mass,
                "predicted_mass_MeV": predicted_mass,
                "relative_mass_error": (
                    predicted_mass - observed_mass
                ) / observed_mass,
            }
        )
    return rows


def normal_form_target():
    """
    The minimal nonlinear target that could select theta = 2/9.

    In singlet-doublet variables, Koide is |A1| = |E|. A C3-equivariant
    reduced-framing phase lock may then fix theta through
    cos(pi*(9*theta - 2)). See phase38_z9_theta_holonomy.py.
    """
    return {
        "singlet_doublet_balance": "|A1| = |E|  -> Koide K = 2/3",
        "z9_closure": "C3 axes x C3 phase/braid sectors -> 9 slots",
        "spinorial_index": "h = 2 (first non-trivial oriented framed closure)",
        "action_origin": "I3=det(B) contains det(Q)=Re(E^3)/27",
        "phase_locking_term": "V_phase ~ 1 - cos(pi*(9 theta - 2))",
        "stationary_condition": "9 theta = 2",
        "theta": THETA_TOPOLOGICAL,
        "status": "Candidate derivation in phase38; full RFG-action theorem still open.",
    }


def status_assessment():
    return {
        "candidate": "C3 cyclic Koide operator",
        "closed": "Mass ratios follow from theta = 2/9 and m ~ nu^2.",
        "theta_candidate": "phase38 Z9 reduced holonomy: theta = 2/9.",
        "open": "Derive the Z9 framed closure directly from the RFG action.",
        "replaces_question": "Why N=5,72,295?",
        "new_question": "Why C3 plus theta=2/9?",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 37: C3 Koide operator")
    print("=" * 72)

    print("\n1. C3 raw frequencies at theta = 2/9")
    for k, value in enumerate(c3_raw_frequencies()):
        print(f"  k={k}: nu_k/A = {value:.12f}")

    print("\n2. Frequency ratios and mass predictions")
    for row in prediction_table():
        print(
            f"  {row['particle']:8s}: "
            f"nu_obs={row['observed_freq_ratio']:.8f}, "
            f"nu_C3={row['predicted_freq_ratio']:.8f}, "
            f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
            f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
            f"rel_err={row['relative_mass_error']:.3e}"
        )

    print("\n3. Koide identity")
    print(f"  K_C3(theta=2/9) = {koide_identity():.12f}")
    print(f"  K_target        = {KOIDE_TARGET:.12f}")

    print("\n4. Axis ratios if L ~ 1/nu")
    axes = axis_ratios_from_c3()
    print(
        "  "
        f"L_e : L_mu : L_tau = "
        f"{axes['electron']:.12f} : {axes['muon']:.12f} : {axes['tau']:.12f}"
    )

    print("\n5. Fit audit")
    fit = best_fit_theta()
    print(f"  theta_topological = {fit['theta_topological']:.12f}")
    print(f"  theta_fit         = {fit['theta_fit']:.12f}")
    print(f"  delta             = {fit['theta_delta']:.3e}")
    print(f"  error(theta=2/9)  = {fit['error_topological']:.3e}")
    print(f"  error(best fit)   = {fit['error_fit']:.3e}")

    print("\n6. Normal-form target")
    for key, value in normal_form_target().items():
        print(f"  {key:24s}: {value}")

    print("\n7. Status")
    for key, value in status_assessment().items():
        print(f"  {key:18s}: {value}")
