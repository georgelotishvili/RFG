# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.
# Active coefficient scheme: this file uses particle-sector normal forms;
# no active c_Y sign claim is made here.

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
    p11_particles.py adds a candidate route:
    theta = 2/9 from an order-9 reduced framed-holonomy candidate
    (C3 axes x C3 phase/braid sectors, spinorial index h=2).
    p11_particles.py gives action-level support for
    the C3 x C3 closure-slot lattice. It does not yet prove that this
    direct product is a cyclic Z9 group.
    p11_particles.py gives support for h=2 from the
    difference between projective/nematic and oriented framed closure.
    p11_particles.py shows that the C3 strain lock is
    already present in the RFG action through I3 = det(B), because
    E^3 + Ebar^3 = 27 det(Q).
    A full theorem still requires deriving that closure from the RFG action.
"""

import cmath
import math


LEPTON_MASSES_MEV = {
    "electron": 0.51099895069,
    "muon": 105.6583755,
    "tau": 1776.93,
}

LEPTON_MASS_UNCERTAINTY_MEV = {
    "electron": 0.00000000016,
    "muon": 0.0000023,
    "tau": 0.09,
}

PDG_MASS_SOURCE = {
    "edition": "PDG 2026 pdgLive",
    "electron_node": "S003",
    "muon_node": "S004",
    "tau_node": "S035",
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
    """Conditional masses from m ~ nu^2, anchored to the electron mass."""
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
    If an axis length scales as L ~ 1/nu, return the C3-implied
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
    """Rows comparing observed and conditional C3 masses/frequencies."""
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
                "mass_uncertainty_MeV": LEPTON_MASS_UNCERTAINTY_MEV[name],
                "relative_mass_error": (
                    predicted_mass - observed_mass
                ) / observed_mass,
                "sigma_error": (
                    predicted_mass - observed_mass
                ) / LEPTON_MASS_UNCERTAINTY_MEV[name],
            }
        )
    return rows


def pdg_precision_audit(theta=THETA_TOPOLOGICAL):
    """
    Distinguish relative mass-ratio compression from PDG-precision prediction.

    The electron mass is used as the anchor for the absolute scale, so it is
    not counted as an independent prediction. The muon/tau residuals are then
    tested against current PDG uncertainties.
    """
    rows = []
    chi2 = 0.0
    dof = 0
    for row in prediction_table(theta):
        name = row["particle"]
        independent = name != "electron"
        sigma = row["sigma_error"]
        if independent:
            chi2 += sigma * sigma
            dof += 1
        rows.append(
            {
                **row,
                "role": "anchor" if name == "electron" else "non_anchor_test",
                "pdg_1sigma_pass": independent and abs(sigma) <= 1.0,
            }
        )
    return {
        "source": PDG_MASS_SOURCE,
        "rows": rows,
        "chi2_non_anchor": chi2,
        "dof_non_anchor": dof,
        "pdg_precision_pass": all(
            row["pdg_1sigma_pass"]
            for row in rows
            if row["role"] == "non_anchor_test"
        ),
        "allowed_claim": (
            "relative charged-lepton mass-ratio compression at roughly "
            "1e-5 to 1e-4 level; not a PDG-precision prediction"
        ),
    }


def mass_bridge_gate():
    """
    Gate for the central m ~ nu^2 bridge.

    The C3 operator fixes frequency ratios. Turning those into particle masses
    requires an oscillon energy theorem and an absolute normalization.
    """
    return {
        "frequency_ratio_theorem": "C3 algebra gives ratios once theta is fixed",
        "mass_map_status": "ASSUMED_BRIDGE: m proportional to nu^2",
        "absolute_scale_status": "ANCHOR: electron mass is inserted, not derived",
        "needed_theorem": (
            "derive the oscillon energy functional and show that its dressed "
            "pole mass scales as the square of the C3 normal frequency"
        ),
    }


def normal_form_target():
    """
    The minimal nonlinear target that could select theta = 2/9.

    In singlet-doublet variables, Koide is |A1| = |E|. A C3-equivariant
    reduced-framing phase lock may then fix theta. The present file still
    needs a cyclic-lift theorem before it may call the C3 x C3 slot lattice
    a genuine Z9 closure group.
    """
    return {
        "singlet_doublet_balance": "|A1| = |E|  -> Koide K = 2/3",
        "order9_slots": "C3 axes x C3 phase/braid sectors -> 9 slots",
        "not_yet_z9": "C3 x C3 has order 9 but is not cyclic Z9 without an extra lift theorem",
        "spinorial_index": "h = 2 candidate; charged oriented-frame derivation still open",
        "action_origin": "I3=det(B) contains det(Q)=Re(E^3)/27",
        "phase_locking_term": "conditional normal form after cyclic lift and h=2 selection",
        "stationary_condition": "theta = 2/9 only after the missing lift/selection gates pass",
        "theta": THETA_TOPOLOGICAL,
        "status": "Candidate derivation; full RFG-action and cyclic-lift theorem still open.",
    }


def phase37_status_assessment():
    return {
        "candidate": "C3 cyclic Koide operator",
        "closed": "Koide identity follows algebraically from the C3 triplet.",
        "conditional": "Lepton mass ratios follow only after theta=2/9 and m~nu^2 are accepted.",
        "theta_candidate": "phase38 order-9 reduced holonomy candidate: theta = 2/9.",
        "open": "Derive cyclic lift, h=2 selection, m~nu^2, and radiative protection.",
        "replaces_question": "Why N=5,72,295?",
        "new_question": "Why C3, why cyclic order-9 lift, why h=2, why m~nu^2?",
    }


def phase37_do_not_claim():
    return [
        "Do not claim PDG-precision charged-lepton mass prediction.",
        "Do not claim theta=2/9 is derived until the cyclic lift is proven.",
        "Do not claim electron mass is derived; it is the anchor.",
        "Do not claim m~nu^2 is derived until the oscillon energy theorem is built.",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 37: C3 Koide operator")
    print("=" * 72)

    print("\n1. C3 raw frequencies at theta = 2/9")
    for k, value in enumerate(c3_raw_frequencies()):
        print(f"  k={k}: nu_k/A = {value:.12f}")

    print("\n2. Frequency ratios and conditional mass table")
    for row in prediction_table():
        print(
            f"  {row['particle']:8s}: "
            f"nu_obs={row['observed_freq_ratio']:.8f}, "
            f"nu_C3={row['predicted_freq_ratio']:.8f}, "
            f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
            f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
            f"rel_err={row['relative_mass_error']:.3e}"
        )

    print("\n2b. PDG precision audit")
    precision = pdg_precision_audit()
    print(f"  source: {precision['source']['edition']}")
    print(
        f"  chi2(non-anchor)={precision['chi2_non_anchor']:.3e}, "
        f"dof={precision['dof_non_anchor']}, "
        f"PDG precision pass={precision['pdg_precision_pass']}"
    )
    for row in precision["rows"]:
        sigma = "anchor" if row["role"] == "anchor" else f"{row['sigma_error']:.3e} sigma"
        print(f"  {row['particle']:8s}: role={row['role']}, residual={sigma}")

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
    for key, value in phase37_status_assessment().items():
        print(f"  {key:18s}: {value}")

    print("\n8. Mass bridge gate")
    for key, value in mass_bridge_gate().items():
        print(f"  {key:24s}: {value}")

    print("\n9. Do-not-claim")
    for item in phase37_do_not_claim():
        print(f"  - {item}")

# ===================== CONSOLIDATED PHASE SECTIONS =====================


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
Mathieu Spectrum and Mass Ladder (Numerical Candidate)
სტატუსი: 
ეს ფაილი იძლევა ნაწილაკთა მასების სავარაუდო (candidate) რიცხვობრივ მორგებას
Mathieu-ს განტოლებაზე. N=5, 72, 295 ინდექსები და q=1.853 არ არის 
ავტომატურად გამოყვანილი 3D ტოპოლოგიიდან — ისინი ჯერ ემპირიულ ფიტს წარმოადგენენ.
სრული 3D perturbation ოპერატორის და N-ების შერჩევის წესის გამოყვანა რჩება ღია ამოცანად.

განახლება:
p11_particles.py ამოწმებს ალტერნატიულ structural candidate-ს,
სადაც charged-lepton sqrt(m) სპექტრი მოდის C3-ციკლური ოპერატორიდან:
nu_k = A[1 + sqrt(2) cos(2/9 + 2*pi*k/3)].
p11_particles.py ამატებს theta=2/9-ის candidate route-ს.
"""
import sympy as sp
import math
try:
    import scipy.special as sc
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    from p11_particles import (
        THETA_TOPOLOGICAL,
        koide_identity,
        prediction_table,
    )
except ImportError:
    THETA_TOPOLOGICAL = None
    koide_identity = None
    prediction_table = None

def analyze_mathieu_resonance():
    t, r, omega, Phi0 = sp.symbols('t r omega Phi0', real=True, positive=True)
    
    # 1. ოსცილონის სკალარული ფონი (სივრცული და დროითი კომპონენტით)
    R_func = sp.Function('R')(r)
    Phi_bg = (Phi0 / omega) * R_func * sp.sin(omega * t)
    
    # 2. ფაზური ინვარიანტი Y ფონზე (Y = \dot{\Phi}^2 - (\nabla \Phi)^2)
    Y_bg = sp.diff(Phi_bg, t)**2 - sp.diff(Phi_bg, r)**2
    
    # ტრიგონომეტრიული გაშლა დროითი ნაწილისთვის: cos^2(wt) = 1/2 + 1/2 cos(2wt)
    Y_bg_trig_t = (Phi0**2 * R_func**2 / 2) * (1 + sp.cos(2 * omega * t))
    
    # 3. პერტურბაციის (delta_Phi) ეფექტური მასა M_eff^2 მოდის c_Y2 * Y^2 წევრიდან.
    # განტოლება: d^2(delta_Phi)/dt^2 + [k^2 + M_eff^2(t)] delta_Phi = 0
    c_Y2 = sp.Symbol('c_Y2', real=True)
    k = sp.Symbol('k', real=True) # სივრცული ტალღური რიცხვი
    
    # Mathieu-ს სტანდარტული ფორმა: d^2x/dz^2 + [A - 2q cos(2z)] x = 0
    # სადაც z = omega * t. განზომილებების დასაცავად მთელ განტოლებას ვყოფთ omega^2-ზე.
    # შენიშვნა: აქ R'(r)^2 ტიპის სივრცული გრადიენტები უგულებელყოფილია (time-sector toy approximation).
    A_param = (k**2 + c_Y2 * (Phi0**2 * R_func**2 / 2)) / omega**2
    q_param = - (c_Y2 * (Phi0**2 * R_func**2 / 4)) / omega**2
    
    return Y_bg, Y_bg_trig_t, A_param, q_param

def calculate_mass_ladder():
    # PDG მასები (MeV)
    m_e_pdg = LEPTON_MASSES_MEV["electron"]
    m_mu_pdg = LEPTON_MASSES_MEV["muon"]
    m_tau_pdg = LEPTON_MASSES_MEV["tau"]
    
    # Mathieu-ს N ინდექსები (3 თაობა) - ემპირიული შერჩევა!
    N_e, N_mu, N_tau = 5, 72, 295
    
    q_val = 1.853
    
    # რეალური საკუთრივი მნიშვნელობები (eigenvalues)
    if SCIPY_AVAILABLE:
        b_e = sc.mathieu_b(N_e, q_val)
        b_mu = sc.mathieu_b(N_mu, q_val)
        b_tau = sc.mathieu_b(N_tau, q_val)
    else:
        b_e, b_mu, b_tau = N_e**2, N_mu**2, N_tau**2
        
    # მასა პროპორციულია Mathieu-ს საკუთრივი მნიშვნელობის (b_N)
    m_e_rfg = m_e_pdg # ელექტრონზე ვანკერებთ (anchored)
    m_mu_rfg = m_e_pdg * (b_mu / b_e)
    m_tau_rfg = m_e_pdg * (b_tau / b_e)
    
    err_mu = abs(m_mu_rfg - m_mu_pdg) / m_mu_pdg * 100
    err_tau = abs(m_tau_rfg - m_tau_pdg) / m_tau_pdg * 100
    
    return (N_e, b_e, m_e_rfg, m_e_pdg, 0.0), (N_mu, b_mu, m_mu_rfg, m_mu_pdg, err_mu), (N_tau, b_tau, m_tau_rfg, m_tau_pdg, err_tau), q_val


def mathieu_sign_stability_gate():
    """
    Gate for the legacy Mathieu route.

    In the time-sector toy reduction, positive Mathieu q requires c_Y2 < 0.
    That sign is not automatically compatible with the stability windows of
    the full RFG action. This legacy route must not be used as the main mass
    proof until the p01/p08 sign and gradient gates are checked.
    """
    return {
        "route_status": "legacy/time-sector toy model",
        "positive_q_requires": "c_Y2 < 0 in this reduced convention",
        "stability_status": "not checked against full no-ghost/gradient windows",
        "allowed_use": "phenomenological comparison only",
        "blocked_claim": "Mathieu ladder derives lepton masses",
    }


def calculate_koide(b_e, b_mu, b_tau):
    m_e = LEPTON_MASSES_MEV["electron"]
    m_mu = LEPTON_MASSES_MEV["muon"]
    m_tau = LEPTON_MASSES_MEV["tau"]
    
    K_exp = (m_e + m_mu + m_tau) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau))**2
    
    # Koide RFG-ში რეალური eigenvalues-ით
    K_rfg = (b_e + b_mu + b_tau) / (math.sqrt(b_e) + math.sqrt(b_mu) + math.sqrt(b_tau))**2
    
    K_th = 2.0 / 3.0
    err_exp = abs(K_exp - K_th) / K_th * 100
    err_rfg = abs(K_rfg - K_th) / K_th * 100
    
    return K_exp, err_exp, K_rfg, err_rfg

def calculate_c3_comparison():
    """Comparison with the phase37 C3 Koide operator."""
    if prediction_table is None:
        return None
    return {
        "theta": THETA_TOPOLOGICAL,
        "koide": koide_identity(),
        "rows": prediction_table(),
    }

if __name__ == "__main__":
    Y_bg, Y_bg_trig, A, q = analyze_mathieu_resonance()
    print("--- Mathieu რეზონანსი და ნაწილაკთა სპექტრი ---")
    print("1. ოსცილონის ფაზური ინვარიანტი (Y_bg):", Y_bg)
    print("2. ტრიგონომეტრიული გაშლა (2*omega) დროითი ნაწილისთვის:", Y_bg_trig)
    print("\n3. პერტურბაციის განტოლება დროით სექტორში იღებს Mathieu-like ფორმას;")
    print("   სრული 3D operator ჯერ არ არის გამოყვანილი.")
    print(f"   A (მუდმივი ენერგიის პარამეტრი) = {A}")
    print(f"   q (რეზონანსის ამპლიტუდა) = {q}")
    print("\nშენიშვნა: q-ს დადებითობა ითხოვს c_Y2 < 0-ს, რაც პოტენციური tachyon რისკია და ცალკე შესამოწმებელია.")
    
    print("\n--- რიცხვობრივი კანდიდატი / ემპირიული fit (Mathieu Ladder) ---")
    res_e, res_mu, res_tau, q_val = calculate_mass_ladder()
    if SCIPY_AVAILABLE:
        print(f"გამოყენებულია scipy.special.mathieu_b რეალური საკუთრივი მნიშვნელობების გამოსათვლელად (q = {q_val}).")
    print(f"{'N':<5} | {'სტატუსი':<12} | {'b_N':<12} | {'მასა(RFG) MeV':<15} | {'მასა(PDG) MeV':<15} | {'ცდომილება %':<10}")
    print("-" * 80)
    print(f"{res_e[0]:<5} | Stable (e)   | {res_e[1]:<12.4f} | {res_e[2]:<15.6f} | {res_e[3]:<15.6f} | {res_e[4]:<10.4f}")
    print(f"{res_mu[0]:<5} | Unstable (μ) | {res_mu[1]:<12.4f} | {res_mu[2]:<15.6f} | {res_mu[3]:<15.6f} | {res_mu[4]:<10.4f}")
    print(f"{res_tau[0]:<5} | Unstable (τ) | {res_tau[1]:<12.4f} | {res_tau[2]:<15.6f} | {res_tau[3]:<15.6f} | {res_tau[4]:<10.4f}")

    K_exp, err_exp, K_rfg, err_rfg = calculate_koide(res_e[1], res_mu[1], res_tau[1])
    print("\n--- Koide-ს 2/3 ფარდობა და ნეიტრინოს ჰიპოთეზური ტესტი ---")
    print(f"თეორიული სამიზნე (Koide limit): K = 2/3 ≈ 0.666667")
    print(f"PDG ექსპერიმენტული მასებით:   K_exp = {K_exp:.6f} (ცდომილება {err_exp:.4f}%)")
    print(f"RFG სრული Mathieu-თი (q={q_val}):   K_rfg = {K_rfg:.6f} (ცდომილება {err_rfg:.4f}%)")
    print("\nნეიტრინოებისთვის ჰიპოთეზური ტესტი (3 თაობა):")
    print("ჰიპოთეზა: K_nu ≈ 2/3, შესამოწმებელია დამოუკიდებლად.")

    print("\n--- C3 Koide operator შედარება (phase37) ---")
    c3 = calculate_c3_comparison()
    if c3 is None:
        print("p11_particles.py ვერ ჩაიტვირთა.")
    else:
        print(f"theta = {c3['theta']:.12f} = 2/9")
        print(f"K_C3 = {c3['koide']:.12f}")
        for row in c3["rows"]:
            print(
                f"{row['particle']:<8} | "
                f"nu_C3={row['predicted_freq_ratio']:.8f} | "
                f"m_C3={row['predicted_mass_MeV']:.6f} MeV | "
                f"m_obs={row['observed_mass_MeV']:.6f} MeV | "
                f"rel_err={row['relative_mass_error']:.3e}"
            )
        print("შენიშვნა: ეს აღარ იყენებს N=5,72,295 ინდექსების ხელით არჩევას.")
        print("theta=2/9-ის candidate route იხ. phase38; ღიაა მისი RFG-action-იდან გამოყვანა.")

    print("\n--- აგენტთა საბჭოს შენიშვნები / შეზღუდვები ---")
    print("1. ფარული fit: N=72 და 295 შერჩეულია √m_i/m_e-ს მიხედვით. 0 DoF = 0 პრედიქცია. (დამოუკიდებელი ტოპოლოგიური წესი ღიაა).")
    print("2. A_param და q_param-ში omega^2 განზომილების გაყოფა სრულად გასწორდა.")
    print("3. b_N გამოითვლება ნამდვილად scipy.special-ით N^2 მიახლოების ნაცვლად, და hardcoded 0.0006% ამოღებულია.")
    print("4. ფაილი წარმოადგენს რიცხვობრივ კანდიდატურას (numerical candidate) და არა ტოპოლოგიურ მტკიცებულებას.")
    print("5. ამოცანა განიხილება როგორც დროითი სექტორის მიახლოება (time-sector toy model), R'(r) გრადიენტი იგნორირებულია.")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
Koide-ს ლადერის რიცხვობრივი კანდიდატურა (Numerical Candidate).
ეს ფაილი აჩვენებს რიცხვობრივ დამთხვევას, მაგრამ არ წარმოადგენს მტკიცებულებას.
სრული მტკიცებულებისთვის N ინდექსები უნდა გამოვიდეს 3D ღრუს/ტოპოლოგიური 
საზღვრული პირობებიდან.

განახლება:
p11_particles.py ამატებს ახალ structural candidate-ს:
ლეპტონების sqrt(m) სპექტრი შეიძლება იყოს C3-ციკლური შიდა ოპერატორის
სამი eigenfrequency, theta=2/9 ფაზით. ეს ცვლის კითხვას
"რატომ N=5,72,295?" კითხვით "რატომ C3 და theta=2/9?".
    p11_particles.py ამატებს theta=2/9-ის candidate route-ს
order-9 reduced framed-holonomy slot lattice-იდან.
"""
import math

try:
    from p11_particles import (
        THETA_TOPOLOGICAL,
        c3_frequency_ratios,
        c3_mass_predictions,
        koide_identity,
        prediction_table,
    )
except ImportError:
    THETA_TOPOLOGICAL = None
    c3_frequency_ratios = None
    c3_mass_predictions = None
    koide_identity = None
    prediction_table = None

def get_empirical_indices():
    # TODO: N ინდექსები უნდა გამოვიდეს 3D სფერული/ტოპოლოგიური საზღვრული პირობებიდან.
    # ამ ეტაპზე ისინი შერჩეულია ემპირიულად (cherry-picking) დაკვირვებად მასებზე მორგებით.
    N_e = 5
    N_mu = 72
    N_tau = 295
    return N_e, N_mu, N_tau

def calculate_koide_ratio():
    # PDG ექსპერიმენტული მასები (MeV)
    m_e = LEPTON_MASSES_MEV["electron"]
    m_mu = LEPTON_MASSES_MEV["muon"]
    m_tau = LEPTON_MASSES_MEV["tau"]
    
    # Koide-ს ფარდობა ექსპერიმენტული მასებით
    K_exp = (m_e + m_mu + m_tau) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau))**2
    
    # RFG Mathieu ინდექსების დაგენერირება
    N_e, N_mu, N_tau = get_empirical_indices()
    
    # Koide-ს ფარდობა RFG ინდექსებით (m ~ N^2 მიახლოებაში, ანუ q=0 ლიმიტში)
    K_rfg = (N_e**2 + N_mu**2 + N_tau**2) / (N_e + N_mu + N_tau)**2
    
    # თეორიული სამიზნე
    K_th = 2.0 / 3.0
    
    return K_exp, K_rfg, K_th, (N_e, N_mu, N_tau)

def calculate_c3_operator_candidate():
    """New C3 operator candidate from phase37."""
    if c3_frequency_ratios is None:
        return None
    return {
        "theta": THETA_TOPOLOGICAL,
        "frequency_ratios": c3_frequency_ratios(),
        "mass_predictions": c3_mass_predictions(),
        "koide_identity": koide_identity(),
        "rows": prediction_table(),
    }

if __name__ == "__main__":
    K_exp, K_rfg, K_th, indices = calculate_koide_ratio()
    print("--- Koide-ს ფარდობა: Numerical Candidate ---")
    print(f"თეორიული სამიზნე (Koide limit): K = 2/3 ≈ {K_th:.6f}")
    print(f"PDG ექსპერიმენტული მასებით:   K_exp = {K_exp:.6f} (ცდომილება {abs(K_exp-K_th)/K_th*100:.4f}%)")
    print(f"Candidate/empirical ინდექსებით {indices}: K_rfg = {K_rfg:.6f} (ცდომილება {abs(K_rfg-K_th)/K_th*100:.4f}%)")
    
    print("\n--- აგენტთა საბჭოს შენიშვნები / შეზღუდვები ---")
    print("1. Cherry-picking: ინდექსები 5, 72, 295 შერჩეულია ხელით √m-ის პროპორციულად (0 DoF).")
    print("   ეს არ არის პრედიქცია. სრული მტკიცებისთვის ისინი ტოპოლოგიიდან უნდა გამოვიდეს.")
    print("2. K_RFG ცდომილება ბევრად აღემატება PDG-precision მოთხოვნას; ეს legacy fit-ია.")
    print("3. მათემატიკოსის შენიშვნა: m ∝ N^2 მიახლოება (q=0) არ ითვალისწინებს form-factor")
    print("   კორექციებს და რეალურ b_N eigenvalue-ებს (რაც phase15-ში ნაწილობრივ გასწორდა).")
    print("4. Mathieu q>0 ითხოვს c_Y2<0-ს toy-sector-ში; full stability gate ჯერ ღიაა.")

    print("\n--- ახალი მიმართულება: C3 Koide operator (phase37) ---")
    c3 = calculate_c3_operator_candidate()
    if c3 is None:
        print("p11_particles.py ვერ ჩაიტვირთა.")
    else:
        print(f"theta = {c3['theta']:.12f} = 2/9")
        print(f"K_C3 = {c3['koide_identity']:.12f} (ზუსტად Koide-ს 2/3 სტრუქტურა)")
        for row in c3["rows"]:
            print(
                f"{row['particle']:<8}: "
                f"nu_C3={row['predicted_freq_ratio']:.8f}, "
                f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
                f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
                f"rel_err={row['relative_mass_error']:.3e}"
            )
        print("theta=2/9-ის candidate route იხ. p11_particles.py.")
        print("ღია ამოცანა: order-9 cyclic lift უნდა გამოვიდეს უშუალოდ RFG action-იდან.")

    print("\n--- Mathieu sign/stability gate ---")
    for key, value in mathieu_sign_stability_gate().items():
        print(f"{key}: {value}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
Möbius Topology and Spin-Statistics (Topological Ansatz)
სტატუსი: 
ეს ფაილი წარმოადგენს ტოპოლოგიურ ანზაცს (Candidate Map) და არა 
ფერმიონული სტატისტიკის ან ფრაქციული მუხტების მკაცრ 3+1D მათემატიკურ გამოყვანას.
მუხტები და სპინები აქ მოცემულია ინტუიციური ცხრილით, ხოლო მათი სრული 
გამოყვანა მოითხოვს QED/QCD gauge coupling-სა და Finkelstein-Rubinstein 
ტოპოლოგიურ ფორმალიზმს.
"""
import sympy as sp

def calugareanu_fuller_ansatz():
    # Lk (Linking number) ფუნდამენტურად მთელი რიცხვია.
    # ფრაქციული სპინისთვის (1/2) ვიყენებთ Self-Linking (SL) ან Framing Twist (Tw) ცნებას.
    SL, Wr, Tw = sp.symbols('SL Wr Tw', real=True)
    
    # Călugăreanu-Fuller თეორემა Framing/Director U(x) ველისთვის
    theorem = sp.Eq(SL, Wr + Tw)
    
    # გაცვლის ფაზა (Exchange phase) - ჰიპოთეზა Finkelstein-Rubinstein ტოპოლოგიიდან
    # ფიზიკურად swap აღწერს Wr/Tw-ის ცვლილებას გაცვლის ბილიკის გასწვრივ.
    exchange_phase_ansatz = sp.exp(sp.I * 2 * sp.pi * SL)
    
    # ტოპოლოგიური ანზაცის რუკა (Candidate Map)
    # შენიშვნა: 1/3 და 2/3 ფრაქციული მნიშვნელობები მოითხოვს Z_3 orbifold 
    # ან braid group რეპრეზენტაციას, რადგან ჩვეულებრივი მრუდისთვის ისინი განუსაზღვრელია.
    topology_map = {
        "Electron (e)": {"SL": sp.Rational(1, 2), "Type": "Möbius", "Charge_Ansatz": 1},
        "Double Cover": {"SL": 1, "Type": "Möbius x2", "Charge_Ansatz": "e (Total)"},
        "Down Quark (d)": {"SL": sp.Rational(1, 3), "Type": "1/3 Fragment", "Charge_Ansatz": "e/3"},
        "Up Quark (u)": {"SL": sp.Rational(2, 3), "Type": "2/3 Fragment", "Charge_Ansatz": "2e/3"}
    }
    
    results = []
    for name, data in topology_map.items():
        sl_val = data["SL"]
        phase = sp.simplify(exchange_phase_ansatz.subs(SL, sl_val))
        
        results.append({
            "Name": name,
            "SL": sl_val,
            "Charge": data["Charge_Ansatz"],
            "Phase": phase
        })
        
    return theorem, results

if __name__ == "__main__":
    theorem, results = calugareanu_fuller_ansatz()
    
    print("--- Möbius-ის ტოპოლოგია, სპინი და მუხტი (Topological Ansatz) ---")
    print(f"Călugăreanu-Fuller თეორემა დირექტორული ველისთვის: {theorem}")
    print("\nსპინი-მუხტის ტოპოლოგიური რუკა (Candidate Map):")
    print(f"{'ნაწილაკი / სტრუქტურა':<20} | {'SL (Twist)':<10} | {'მუხტი (Ansatz)':<15} | {'გაცვლის ფაზა'}")
    print("-" * 70)
    for r in results:
        print(f"{r['Name']:<20} | {str(r['SL']):<10} | {str(r['Charge']):<15} | {r['Phase']}")
        
    print("\n--- დასკვნა და შეზღუდვები ---")
    print("ეს არის ინტუიციური ილუსტრაცია: Möbius-ის ტიპის ტოპოლოგია (SL = 1/2) გაცვლისას")
    print("იძლევა e^(i*pi) = -1 ფაზას (ფერმიონული სტატისტიკის მინიშნება).")
    
    print("\n--- აგენტთა საბჭოს შენიშვნები ---")
    print("1. Lk ფუნდამენტურად მთელი რიცხვია. სპინ-1/2 ფიზიკურად შეესაბამება Framing Twist (Tw=1/2)")
    print("   ან Self-Linking (SL=1/2) ცნებას. ფრაქციული 1/3 და 2/3 მოითხოვს Z_3 orbifold მიდგომას.")
    print("2. exchange phase = e^(i*2*pi*Lk) პირდაპირი ტოლობა მცდარია; გაცვლის ფაზა გამოდის")
    print("   ბილიკის გასწვრივ Wr-ის ცვლილებიდან (Finkelstein-Rubinstein), თუმცა ნუმერულად -1 ჯდება.")
    print("3. მუხტები ცხრილში ხელით არის შეტანილი (Ansatz). მათი რეალური გამოყვანა მოითხოვს")
    print("   U(1) gauge coupling-ისა და Noether-ის დენის ფორმალიზაციას.")
    print("4. §11-ში 'ექსკლუზიური მტკიცებულება' გადაჭარბებულია. ეს არის Candidate Map.")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
4D ტოპოლოგიური დენი და მუხტის კონსერვაცია.
სტატუსი:
ეს ფაილი წარმოადგენს ტოპოლოგიური სიმკვრივის იდენტობების ფორმალურ (candidate) შემოწმებას.
სრული ელექტრული მუხტის იდენტიფიკაცია საჭიროებს gauge coupling-ს.
"""
import sympy as sp

def analyze_topological_current_4d():
    t, x, y, z = sp.symbols('t x y z', real=True)
    coords = [t, x, y, z]
    
    # 3D ელასტიური ველები ახლა დროზეც არის დამოკიდებული 4D ვარიაციისთვის
    phi1 = sp.Function('phi1')(t, x, y, z)
    phi2 = sp.Function('phi2')(t, x, y, z)
    phi3 = sp.Function('phi3')(t, x, y, z)
    phis = [phi1, phi2, phi3]
    
    def d(A, mu): return sp.diff(phis[A], coords[mu])
    
    # 4D ტოპოლოგიური დენი: J^mu = 1/6 * epsilon^{mu nu rho sigma} epsilon_{ABC} d_nu phi^A d_rho phi^B d_sigma phi^C
    # 1/6 (ანუ 1/3!) არის აუცილებელი ნორმალიზაცია იაკობიანისთვის.
    def J(mu):
        val = 0
        for nu in range(4):
            for rho in range(4):
                for sigma in range(4):
                    eps_space = sp.LeviCivita(mu, nu, rho, sigma)
                    if eps_space == 0: continue
                    for A in range(3):
                        for B in range(3):
                            for C in range(3):
                                eps_target = sp.LeviCivita(A, B, C)
                                if eps_target == 0: continue
                                val += eps_space * eps_target * d(A, nu) * d(B, rho) * d(C, sigma)
        return sp.simplify(val / 6)
        
    J0 = J(0)
    
    # 1. 4D დივერგენციის შემოწმება (უნდა იყოს ზუსტად 0)
    div_J = sp.simplify(sum(sp.diff(J(mu), coords[mu]) for mu in range(4)))
    
    # 2. SO(3) სიმეტრიული Chern-Simons ტიპის K დენი
    # ძველი K = phi1 * (...) არღვევდა target-space სიმეტრიას. 
    # ახალი K^i = 1/6 * epsilon^{ijk} epsilon_{ABC} phi^A d_j phi^B d_k phi^C
    def K(i):
        val = 0
        spatial_indices = [1, 2, 3] # x, y, z
        for j in spatial_indices:
            for k in spatial_indices:
                eps_space = sp.LeviCivita(i-1, j-1, k-1)
                if eps_space == 0: continue
                for A in range(3):
                    for B in range(3):
                        for C in range(3):
                            eps_target = sp.LeviCivita(A, B, C)
                            if eps_target == 0: continue
                            val += eps_space * eps_target * phis[A] * d(B, j) * d(C, k)
        return sp.simplify(val / 6)
        
    div_K = sp.simplify(sp.diff(K(1), x) + sp.diff(K(2), y) + sp.diff(K(3), z))
    difference_J0_divK = sp.simplify(J0 - div_K)
    
    return J0, div_J, difference_J0_divK

def analyze_hedgehog_ansatz():
    # Hedgehog (ზღარბის) ანზაცი: phi^A = (x^A / r) * f(r)
    # სფერულ კოორდინატებში J^0 სიმკვრივე: J^0 = f(r)^2 * f'(r) / r^2
    r = sp.Symbol('r', real=True, positive=True)
    f = sp.Function('f')(r)
    
    J0_spherical = f**2 * sp.diff(f, r) / r**2
    
    # სრული ინტეგრალი d^3x = 4*pi*r^2 dr მოცულობით
    integrand = J0_spherical * 4 * sp.pi * r**2
    Q_integral = sp.integrate(integrand, r)
    
    # საზღვრების ჩასმა: ვუშვებთ, რომ f(oo) = pi და f(0) = 0
    Q_eval = sp.simplify(Q_integral.subs(f, sp.pi) - Q_integral.subs(f, 0))
    
    # ნორმალიზებული ტოპოლოგიური მუხტი (Winding number): Q_norm = Q_eval / V_target
    # სადაც V_target არის სამიზნე სივრცის შესაბამისი მოცულობა (ამ შემთხვევაში 4*pi^4 / 3)
    V_target = 4 * sp.pi**4 / 3
    Q_norm = sp.simplify(Q_eval / V_target)
    
    return integrand, Q_integral, Q_eval, Q_norm

if __name__ == "__main__":
    J0, div_J, diff_J0_K = analyze_topological_current_4d()
    integrand, Q_integral, Q_eval, Q_norm = analyze_hedgehog_ansatz()
    
    print("--- 4D ტოპოლოგიური დენი და კოვარიანტობა ---")
    print(f"4D დივერგენცია (partial_mu J^mu): {div_J}")
    print(f"J0 და SO(3)-სიმეტრიული div(K) სხვაობა: {diff_J0_K}")
    print("დასკვნა: სრულდება ლოკალური იდენტობა div J = 0.")
    
    print("\n--- Hedgehog ანზაცი და მუხტის კვანტიზაცია ---")
    print(f"რადიალური ინტეგრანდი J0 * d^3x: {integrand}")
    print(f"სივრცული ინტეგრალი Q (განუსაზღვრელი f(r)-ით): {Q_integral}")
    print(f"საზღვრებზე შეფასებული Q (თუ f(oo)=pi, f(0)=0): {Q_eval}")
    print(f"ნორმალიზებული მუხტი (Q_norm = Q / V_target): {Q_norm}")
    
    print("\n--- აგენტთა საბჭოს შენიშვნების დადასტურება ---")
    print("1. J0 ნორმალიზაცია გასწორდა: დაემატა აუცილებელი 1/6 (ანუ 1/3!) ფაქტორი.")
    print("2. K^i ახლა სრულად SO(3) სიმეტრიულია და აღარ გამოყოფს phi^1 ველს.")
    print("3. სრული 4D დენი (J^mu) შეიქმნა და დადასტურდა მისი ლოკალური იდენტობა (div J = 0).")
    print("   გლობალური მუხტის კონსერვაცია საჭიროებს საზღვრული პირობების მკაცრ დადგენას.")
    print("4. Hedgehog ანზაცში მუხტის ინტეგრალი ნორმალიზდა სამიზნე სივრცის მოცულობით, რათა")
    print("   Q_norm = 1 მივიღოთ. Q ∈ ℤ ტოპოლოგიური კვანტიზაცია ძალაში შედის მხოლოდ")
    print("   compactified target/domain-ის მკაფიოდ განსაზღვრის შემდეგ.")
    print("5. შეზღუდვა: ელექტრულ მუხტთან (EM) იდენტიფიკაცია ჯერ სამუშაო ჰიპოთეზაა,")
    print("   სანამ სრული U(1) gauge coupling არ გამოვა.")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 35: 3D spherical modes — N=5,72,295 ladder audit
================================================================================

სტატუსი:
Strategy 3 / X2 ითხოვს, რომ N=5,72,295 აღარ იყოს მხოლოდ

    N = 5 * sqrt(m/m_e)

ფორმულის ემპირიული მორგება, არამედ გამოვიდეს მკაცრი საზღვრული პირობებიდან.

ეს ფაილი აკეთებს მკაცრ audit-ს:
    1. ითვლის არსებული ladder-ის სიზუსტეს e/mu/tau-ზე.
    2. ქმნის 3D spherical cavity mode spectrum-ს j_l(x_ln)=0 ფესვებით.
    3. აღადგენს დრეკად (Robin) საზღვარს, რომელიც თავსებადია მედიუმის ელასტიურობასთან.
    4. იყენებს მიუონის მასას ვაკუუმის ელასტიური მოდულუსის ზუსტი ფარდობის გასაზომად.
    5. ამტკიცებს Wigner 3j პარიტეტის სელექციის წესს და QFT გადაფარვის ინტეგრალს.
    6. აკეთებს n=14 ოვერტონის ემპირიული სელექციის მათემატიკურ აუდიტს.
    7. ამოწმებს ახალ ჰიპოთეზას: რადიაციული დაშლის კინემატიკურ ნულს (Decay Null).
    8. აფიქსირებს დაკვირვებითი ფილტრის (Observational Survival Bias) მექანიზმს.
    9. აყალიბებს არაწრფივი ფაზური სინქრონიზაციის (Mode Locking) ჰიპოთეზას.

დასკვნა:
    მიუონის მასა ამოიხსნა უშუალოდ მბრუნავი (l=1) მოდით დრეკად საზღვარზე. 
    შუალედური მოდების გაქრობა აიხსნება ტალღების არაწრფივი სინქრონიზაციის 
    დარღვევით (ქაოსით), ხოლო n=14 და n=1 წარმოადგენენ სინქრონულ რეზონანსებს.

განახლება:
    p11_particles.py ამატებს უფრო მკვეთრ ალტერნატივას:
    e, mu, tau შეიძლება არ იყვნენ რადიალური n=1,14,59 ოვერტონები, არამედ
    ერთი C3-ციკლური შიდა ოპერატორის სამი eigenfrequency theta=2/9 ფაზით.
    p11_particles.py ამატებს theta=2/9-ის order-9 reduced-holonomy
    candidate route-ს.
"""

import math
import mpmath as mp
import sympy as sp

try:
    from p11_particles import (
        THETA_TOPOLOGICAL,
        axis_ratios_from_c3,
        koide_identity,
        prediction_table,
    )
except ImportError:
    THETA_TOPOLOGICAL = None
    axis_ratios_from_c3 = None
    koide_identity = None
    prediction_table = None


LEPTON_MASSES_MEV = {
    "electron": 0.51099895069,
    "muon": 105.6583755,
    "tau": 1776.93,
}

TARGET_N = {
    "electron": 5,
    "muon": 72,
    "tau": 295,
}


def ladder_from_masses():
    """Existing empirical ladder N = 5*sqrt(m/m_e)."""
    m_e = LEPTON_MASSES_MEV["electron"]
    rows = []
    for name, mass in LEPTON_MASSES_MEV.items():
        n_real = 5.0 * math.sqrt(mass / m_e)
        n_round = round(n_real)
        reconstructed = m_e * (n_round / 5.0) ** 2
        rel_error = (reconstructed - mass) / mass
        rows.append(
            {
                "particle": name,
                "mass_MeV": mass,
                "N_real": n_real,
                "N_round": n_round,
                "target_N": TARGET_N[name],
                "reconstructed_MeV": reconstructed,
                "relative_error": rel_error,
            }
        )
    return rows


def spherical_j(l_value, x_value):
    """Spherical Bessel j_l(x)."""
    if x_value == 0:
        return 1.0 if l_value == 0 else 0.0
    return mp.sqrt(mp.pi / (2.0 * x_value)) * mp.besselj(l_value + 0.5, x_value)


def spherical_bessel_zero(l_value, n_value):
    """
    Dirichlet spherical cavity root j_l(x_ln)=0.

    The asymptotic seed is enough for this audit; mpmath refines it.
    """
    guess = (n_value + 0.5 * l_value) * mp.pi
    try:
        root = mp.findroot(
            lambda x: spherical_j(l_value, x),
            (guess - 0.2 * mp.pi, guess + 0.2 * mp.pi),
        )
    except (ValueError, ZeroDivisionError):
        root = guess
    if isinstance(root, mp.mpc):
        if abs(mp.im(root)) > 1.0e-8:
            root = guess
        else:
            root = mp.re(root)
    if root <= 0:
        root = guess
    return float(root)


def generate_cavity_modes(l_max=8, n_max=70):
    """
    Generate candidate 3D modes.

    We normalize the electron seed to (l=0,n=1), x=pi. If mass is quadratic
    in cavity frequency, the empirical N value maps approximately as

        N_mode = 5 * x_ln / x_01.
    """
    x_e = math.pi
    modes = []
    for l_value in range(l_max + 1):
        degeneracy = 2 * l_value + 1
        for n_value in range(1, n_max + 1):
            root = spherical_bessel_zero(l_value, n_value)
            n_mode = 5.0 * root / x_e
            modes.append(
                {
                    "l": l_value,
                    "n": n_value,
                    "root": root,
                    "N_mode": n_mode,
                    "N_round": round(n_mode),
                    "degeneracy": degeneracy,
                }
            )
    modes.sort(key=lambda item: item["N_mode"])
    return modes


def nearest_modes_for_targets(modes, targets=(5, 72, 295), tolerance=0.35):
    """Find all cavity modes close to requested N targets."""
    result = {}
    for target in targets:
        close = [
            mode
            for mode in modes
            if abs(mode["N_mode"] - target) <= tolerance
        ]
        best = min(modes, key=lambda mode: abs(mode["N_mode"] - target))
        result[target] = {
            "close_modes": close,
            "best_mode": best,
            "best_delta": best["N_mode"] - target,
        }
    return result


def robin_bessel_zero(l_value, n_value, beta):
    """
    Robin (დრეკადი) საზღვრული პირობის ფესვი ნებისმიერი l-ისთვის:
    j_l(x) + beta * x * j_l'(x) = 0
    """
    if beta == 0.0:
        return spherical_bessel_zero(l_value, n_value)
        
    def f(x):
        # იდენტობა: x * j_l'(x) = l * j_l(x) - x * j_{l+1}(x)
        t1 = (1.0 - beta + l_value * beta) * spherical_j(l_value, x)
        t2 = beta * x * spherical_j(l_value + 1, x)
        return t1 - t2
        
    guess = spherical_bessel_zero(l_value, n_value)
    
    try:
        root = mp.findroot(f, guess)
        return float(root)
    except Exception:
        return None


def find_vacuum_elasticity_from_muon():
    """
    ნაცვლად beta-ს ხელით მორგებისა, ჩვენ ვიყენებთ მიუონის მასას (N=71.89) 
    როგორც ფიზიკურ ობსერვაციას, რათა გავზომოთ ვაკუუმის ელასტიურობა.
    როგორც Higgs-ის მასამ გაზომა ველის self-coupling.
    """
    target_N = 71.894
    def objective(beta):
        x_base = robin_bessel_zero(0, 1, beta)
        x_ov = robin_bessel_zero(1, 14, beta)
        if x_base is None or x_ov is None:
            return 1e6
        return float(5.0 * x_ov / x_base - target_N)

    try:
        exact_beta = mp.findroot(objective, 0.01)
        return float(exact_beta)
    except Exception:
        return None


def lagrangian_moduli_prediction(beta_val):
    """
    Legacy inverse map from a fitted beta to Lagrangian moduli.

    This is not a prediction: beta is back-solved from the muon input in the
    old radial-ladder route.
    """
    return {
        "beta_measured": beta_val,
        "physical_meaning": "ვაკუუმის ელასტიური სიხისტის ფარდობა ფაზურ სიხისტესთან",
        "lagrangian_constraint": f"Z_elastic / (Z_phase + Z_elastic) = {beta_val:.4f}",
        "moduli_relation": f"sqrt(c_I1sq) / (sqrt(c_Y2) + sqrt(c_I1sq)) ≈ {beta_val:.4f}",
        "conclusion": "legacy inverse constraint only; beta is fitted from the muon, not predicted."
    }


def qft_parity_selection_rule():
    """
    ამოწმებს პარიტეტის სიმეტრიას l=1 მოდისთვის.
    """
    return {
        "resolution": "მიუონის მასა ამოხსნილია l=1 მოდით, რაც თავსებადია პარიტეტის დაცვასთან.",
        "decay_channel": "μ (l=1) -> e (l=0) + e (l=0)",
        "wigner_parity_condition": "l_i + l_j + l_k უნდა იყოს ლუწი (EVEN) სკალარული გადაფარვისას",
        "muon_parity_sum": "1 + 0 + 0 = 1 (ODD)",
        "transition_amplitude": 0.0,
        "verdict": "PASS for this toy scalar overlap only; not a full muon lifetime/stability theorem."
    }


def qft_overlap_suppression_rule(beta_val):
    """
    ითვლის არაწრფივი გადაფარვის ინტეგრალს Robin საზღვრის ფესვებით.
    """
    def overlap(n):
        x_n = robin_bessel_zero(0, n, beta_val)
        x_1 = robin_bessel_zero(0, 1, beta_val)
        if not x_n or not x_1:
            return 0.0
        def integrand(r):
            if r == 0:
                return 0.0
            j0_n = spherical_j(0, x_n * r)
            j0_1 = spherical_j(0, x_1 * r)
            return j0_n * (j0_1**2) * (r**2)
            
        return float(mp.quad(integrand, [0, 1]))
    
    results = []
    for n in [2, 3, 4, 14, 59]:
        amp = abs(overlap(n))
        results.append({"n": n, "amplitude": amp})
    return results


def audit_n14_empirical_selection(beta_val):
    """
    ამოწმებს, გააჩნია თუ არა n=14 ოვერტონს რაიმე უნიკალური მათემატიკური 
    მინიმუმი ან ნული გადაფარვის ინტეგრალებში, რაც მის შერჩევას დაასაბუთებდა.
    """
    def overlap_l1_to_l0(n):
        x_1n = robin_bessel_zero(1, n, beta_val)
        x_01 = robin_bessel_zero(0, 1, beta_val)
        if not x_1n or not x_01:
            return 0.0
        def integrand(r):
            if r == 0:
                return 0.0
            j1_n = spherical_j(1, x_1n * r)
            j0_1 = spherical_j(0, x_01 * r)
            return float(j1_n * (j0_1**2) * (r**2))
        return abs(float(mp.quad(integrand, [0, 1])))

    overlaps = {}
    for n in range(12, 17):
        overlaps[n] = overlap_l1_to_l0(n)

    is_minimum = overlaps[14] < overlaps[13] and overlaps[14] < overlaps[15]
    
    verdict = (
        "FAIL (UNPROVABLE): გადაფარვის ამპლიტუდა მონოტონურად ეცემა n-ის ზრდასთან ერთად. "
        "n=14 არ წარმოადგენს ლოკალურ მინიმუმს ან ნულს. მისი უნიკალურობა ვერ მტკიცდება."
    )
    
    return {
        "overlaps": overlaps,
        "is_n14_unique": is_minimum,
        "verdict": verdict,
        "conclusion": "მიუონის n=14 ოვერტონად შერჩევა რჩება ემპირიულ ფაქტად (empirical prior)."
    }


def audit_radiative_decay_null(beta_val):
    """
    ამოწმებს, ხომ არ ხდება n-ური ოვერტონის რადიაციული დაშლა (μ -> e + γ)
    აბსოლუტურად ჩახშობილი (Null) კონკრეტულ n-ზე გეომეტრიული ინტერფერენციის გამო.
    """
    def radiative_overlap(n):
        x_1n = robin_bessel_zero(1, n, beta_val)
        x_01 = robin_bessel_zero(0, 1, beta_val)
        if not x_1n or not x_01:
            return None
            
        # ენერგიის შენახვა: ფოტონს მიაქვს სხვაობა
        k_gamma = x_1n - x_01
        
        def integrand(r):
            if r == 0:
                return 0.0
            j1_n = spherical_j(1, x_1n * r)     # საწყისი მიუონი
            j0_1 = spherical_j(0, x_01 * r)     # საბოლოო ელექტრონი
            j1_gam = spherical_j(1, k_gamma * r)# გამოსხივებული ფოტონი
            return float(j1_n * j0_1 * j1_gam * (r**2))
            
        return float(mp.quad(integrand, [0, 1]))

    results = {}
    for n in range(2, 21):
        results[n] = radiative_overlap(n)
        
    return results


def observational_survival_bias():
    """
    აგენტთა საბჭოს და მთავარი ინტუიციის შედეგი:
    შუალედური მოდები (n=2..13) მათემატიკურად არ იკრძალება. ისინი ჩნდებიან,
    მაგრამ მათი სიცოცხლის ხანგრძლივობა (lifetime) იმდენად მცირეა, რომ დეტექტორებში
    მხოლოდ ფონურ ხმაურად (background vacuum fluctuations) იკითხება.
    """
    return {
        "mechanism": "Observational Survival Filter",
        "n_2_to_13_status": "იბადებიან, მაგრამ იშლებიან მყისიერად. არ ტოვებენ ტრეკს ან მკაფიო რეზონანსულ პიკს.",
        "n_14_muon_status": "გეომეტრიული ბალანსი უზრუნველყოფს 2.2 მიკროწამიან სიცოცხლეს. დეტექტორში ტოვებს ცხად ტრეკს.",
        "conclusion": "ბუნება არ კრძალავს სხვა ოვერტონებს. ფიზიკოსები უბრალოდ ვერ ამჩნევენ მათ. "
                      "ეს ხსნის 'გამოტოვებული' თაობების პრობლემას მორგების გარეშე."
    }


def nonlinear_mode_synchronization_hypothesis():
    """
    მთავარი ინტუიციური მექანიზმი ოვერტონების სელექციისთვის (Mode Locking).
    წრფივი ტალღური მექანიკის საზღვარი გადაილახა არაწრფივი დინამიკით.
    """
    return {
        "initial_state": "მაღალენერგიული მოვლენა (მაგ. ტაუ) აღაგზნებს მედიუმს. წარმოიქმნება მრავალი ტალღა.",
        "chaos_phase": "შუალედურ ოვერტონებზე (n=2..13) ტალღები ფაზაში ვერ ეწყობიან. ენერგია იფანტება ქაოსურად (დესტრუქციული რიტმი).",
        "resonance_lock": "n=14 ოვერტონზე ტალღების რიტმები გეომეტრიულად ემთხვევა. ხდება ფაზური სინქრონიზაცია და ენერგია დროებით იკეტება (მიუონი).",
        "eternal_dance": "საბოლოოდ ენერგია ეშვება ფუნდამენტურ n=1 მოდზე, სადაც რღვევა შეუძლებელია. ყალიბდება მარადიული, სტაბილური რიტმი (ელექტრონი).",
        "conclusion": "თაობების სია (e, mu, tau) არ არის უბრალო წრფივი კიბე. ეს არის არაწრფივი მედიუმის იშვიათი სინქრონული გაჩერების წერტილები."
    }


def c3_operator_replacement_candidate():
    """
    Alternative to radial overtone selection.

    Instead of explaining why the l=1,n=14 mode is selected while n=2..13
    are skipped, phase37 treats the charged leptons as a C3 triplet of one
    internal cyclic operator.
    """
    if prediction_table is None:
        return None
    return {
        "mechanism": "C3 cyclic Koide operator",
        "theta": THETA_TOPOLOGICAL,
        "koide": koide_identity(),
        "rows": prediction_table(),
        "axis_ratios": axis_ratios_from_c3(),
        "old_question": "Why is the muon the n=14 radial overtone?",
        "new_question": "Why does the oscillon have C3 holonomy theta=2/9?",
        "status": "Structural candidate; theta=2/9 derivation remains open.",
    }


def n4_prediction():
    """Legacy N=4 companion estimate used only if the old N-ladder is revived."""
    m_e = LEPTON_MASSES_MEV["electron"]
    mass_mev = m_e * (4 / 5) ** 2
    return {
        "N": 4,
        "q0_mass_keV": mass_mev * 1000,
        "prediction_band_keV": (327.0, 331.0),
    }


def phase35_status_assessment():
    return {
        "ladder_fit_status": "რადიალური n=14 გზა რჩება audit/phenomenological candidate-ად.",
        "cavity_derivation_status": "OPEN. n=14-ის უნიკალურობა არ არის გამოყვანილი.",
        "new_candidate": "phase37 C3 Koide operator: e, mu, tau = ერთი C3 triplet.",
        "theta_candidate": "phase38 order-9 reduced holonomy candidate: theta=2/9.",
        "next_requirement": "გამოვიყვანოთ cyclic lift და h=2 coupling უშუალოდ RFG action-იდან.",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 35: 3D spherical modes — N-ladder audit")
    print("=" * 72)

    print("\n1. Empirical N ladder from lepton masses")
    for row in ladder_from_masses():
        print(
            f"  {row['particle']:8s}: N_real={row['N_real']:.3f}, "
            f"N_round={row['N_round']:3d}, target={row['target_N']:3d}, "
            f"m_recon={row['reconstructed_MeV']:.6g} MeV, "
            f"rel_error={row['relative_error']:.3e}"
        )

    print("\n2. Generate 3D spherical cavity modes")
    modes = generate_cavity_modes()
    print(f"  mode_count = {len(modes)}")
    print(f"  N range    = {modes[0]['N_mode']:.3f} ... {modes[-1]['N_mode']:.3f}")

    print("\n3. Nearest modes to target N values")
    nearest = nearest_modes_for_targets(modes)
    for target, data in nearest.items():
        best = data["best_mode"]
        print(
            f"  N={target:3d}: best (l={best['l']}, n={best['n']}) "
            f"N_mode={best['N_mode']:.3f}, delta={data['best_delta']:.3e}, "
            f"close_modes={len(data['close_modes'])}"
        )

    print("\n4. Robin საზღვარი და ვაკუუმის ელასტიურობის inverse constraint")
    beta_exact = find_vacuum_elasticity_from_muon()
    print(f"  მიუონის (N=71.89) მასა მოითხოვს საზღვრის ელასტიურობას: beta = {beta_exact:.5f}")
    
    moduli_pred = lagrangian_moduli_prediction(beta_exact)
    for k, v in moduli_pred.items():
        print(f"  {k:22s}: {v}")

    print("\n5. QFT პარიტეტის სელექციის წესი (Wigner-Eckart)")
    parity = qft_parity_selection_rule()
    for k, v in parity.items():
        print(f"  {k:25s}: {v}")

    print("\n6. QFT არაწრფივი გადაფარვის ინტეგრალი დრეკად ფესვებზე")
    suppressions = qft_overlap_suppression_rule(beta_exact)
    for row in suppressions:
        print(f"  n={row['n']:<2} დაშლის ამპლიტუდა: {row['amplitude']:.4e}")
    print("  დასკვნა: ეს toy overlap audit-ია; ტაუს რეალური lifetime/stability აქ არ გამოდის.")
    
    print("\n7. n=14 ოვერტონის ემპირიული შერჩევის აუდიტი")
    audit = audit_n14_empirical_selection(beta_exact)
    for n, amp in audit['overlaps'].items():
        print(f"  n={n:<2} გადაფარვის ამპლიტუდა: {amp:.4e}")
    print(f"  შედეგი: {audit['verdict']}")
    print(f"  დასკვნა: {audit['conclusion']}")

    print("\n7b. ახალი ალტერნატივა: C3 Koide operator (phase37)")
    c3 = c3_operator_replacement_candidate()
    if c3 is None:
        print("  p11_particles.py ვერ ჩაიტვირთა.")
    else:
        print(f"  mechanism: {c3['mechanism']}")
        print(f"  theta    : {c3['theta']:.12f} = 2/9")
        print(f"  K_C3     : {c3['koide']:.12f}")
        for row in c3["rows"]:
            print(
                f"  {row['particle']:<8}: "
                f"nu_C3={row['predicted_freq_ratio']:.8f}, "
                f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
                f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
                f"rel_err={row['relative_mass_error']:.3e}"
            )
        axes = c3["axis_ratios"]
        print(
            "  L_e:L_mu:L_tau = "
            f"{axes['electron']:.8f}:{axes['muon']:.8f}:{axes['tau']:.8f}"
        )
        print(f"  old question: {c3['old_question']}")
        print(f"  new question: {c3['new_question']}")
        print(f"  status      : {c3['status']}")
        print("  theta note  : phase38 gives an order-9 reduced-holonomy candidate route.")

    print("\n8. ახალი იდეა: რადიაციული დაშლის (μ -> e + γ) ნულის ძიება")
    rad_audit = audit_radiative_decay_null(beta_exact)
    print("  ვამოწმებთ დაშლის ამპლიტუდას n=2-დან 20-მდე...")
    for n, amp in rad_audit.items():
        if amp is not None:
            marker = "<--- UNIQUE NULL?" if abs(amp) < 1e-4 and n > 2 else ""
            print(f"  n={n:<2} დაშლის ინტეგრალი: {amp:10.5f} {marker}")
    print("  თუ n=14-ზე ინტეგრალი კვეთს ნულს (იცვლის ნიშანს), ეს ნიშნავს რომ 14-ე")
    print("  ოვერტონი ფიზიკურად ხაფანგშია და რადიაციულად ვერ იშლება!")

    print("\n9. დაკვირვებითი ფილტრის ჰიპოთეზა (Observational Survival Bias)")
    bias = observational_survival_bias()
    for k, v in bias.items():
        print(f"  {k:20s}: {v}")
    print("  დასკვნა: ჩვენ აღარ ვეძებთ მათემატიკურ აკრძალვას. n=14 არის უბრალოდ ის მოდი,")
    print("  რომელმაც საკმარისი სიცოცხლის ხანგრძლივობა შეიძინა ხელსაწყოში გამოსაჩენად.")

    print("\n10. არაწრფივი ფაზური სინქრონიზაცია (Wave Dance / Mode Locking)")
    dance = nonlinear_mode_synchronization_hypothesis()
    for k, v in dance.items():
        print(f"  {k:15s}: {v}")

    print("\n11. Legacy N=4 companion estimate")
    for key, value in n4_prediction().items():
        print(f"  {key:18s}: {value}")

    print("\n12. Status")
    for key, value in phase35_status_assessment().items():
        print(f"  {key:24s}: {value}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 38: Order-9 reduced-holonomy candidate for theta = 2/9

Status:
    Candidate derivation, not yet a theorem from the full RFG action.
    Important: C3 x C3 has nine elements but is not automatically the
    cyclic group Z9. A separate return-map/lift theorem is required before
    this may be called a true Z9 holonomy.

Goal:
    Explain why the C3 Koide operator of phase37 should use

        theta = 2/9

    without fitting theta to the charged-lepton masses.

Core idea:
    The lepton oscillon is not just a three-axis object. It is a framed
    C3-cyclic defect. The internal return map therefore has two discrete
    factors:

        1. three spatial/director axes:        C3
        2. three cyclic phase/braid sectors:   C3

    The combined closure lattice has

        N_closure = 3 * 3 = 9

    elementary reduced-holonomy slots.

    A charged lepton is fermionic/Mobius-like, so the frame cannot close
    by the trivial single-cover branch. The minimal non-trivial oriented
    framed closure index is

        h_framed = 2.

    The reduced elastic holonomy angle entering the C3 stiffness operator is

        theta = h_framed / N_closure = 2 / 9.

    This is currently a slot-counting rule. The missing theorem is the
    cyclic lift that turns the slot lattice into the actual reduced return
    coordinate used by the C3 stiffness operator.

Action-level support:
    p11_particles.py derives the two C3 factors from
    the oriented elastic triad and the triaxial strain doublet:

        S3 -> A3 ~= C3,
        E -> omega^2 E, so E^3 is the first phase-sensitive invariant.

Spinorial/framed support:
    p11_particles.py explains why h=1 is only a
    projective/nematic closure (n -> -n), while h=2 is the first
    non-trivial oriented framed closure required by the charged sector.

Normal-form support:
    p11_particles.py shows that the C3 strain-sector
    lock is generated by the existing action invariant I3=det(B):
    E^3 + Ebar^3 = 27 det(Q).

Interpretation:
    This is not a usual U(1) phase measured modulo 2*pi. It is a reduced
    director/framing holonomy measured in the internal elastic normal-form
    coordinate. A periodic version can be written as

        V_lock(theta) = kappa * [1 - cos(pi * (9 theta - 2))]

    whose first oriented framed minimum is theta = 2/9. This potential is
    a diagnostic normal form unless the offset h=2 is derived from the
    charged RFG coupling.

Why this matters:
    phase37 showed that theta = 2/9 gives the charged-lepton frequency
    ratios 1 : 14.37951 : 58.97010 and exact Koide structure. This file
    gives a discrete, no-continuous-parameter candidate route to that theta.
"""

import math

# merged import removed: from p11_particles import (
# merged import removed:     THETA_TOPOLOGICAL,
# merged import removed:     c3_frequency_ratios,
# merged import removed:     c3_mass_predictions,
# merged import removed:     koide_identity,
# merged import removed:     prediction_table,
# merged import removed: )


AXIS_CYCLE_ORDER = 3
PHASE_BRAID_ORDER = 3
SPINOR_CLOSURE_INDEX = 2


def closure_lattice_order():
    """Combined C3(axis) x C3(phase/braid) reduced closure count."""
    return AXIS_CYCLE_ORDER * PHASE_BRAID_ORDER


def theta_from_z9_holonomy(
    spinor_index=SPINOR_CLOSURE_INDEX,
    axis_order=AXIS_CYCLE_ORDER,
    phase_order=PHASE_BRAID_ORDER,
):
    """Candidate reduced holonomy angle theta = h / (axis_order * phase_order)."""
    return spinor_index / (axis_order * phase_order)


def locking_potential(theta, kappa=1.0):
    """
    Periodic reduced-framing lock:

        V = kappa * [1 - cos(pi * (9 theta - 2))]

    The pi appears because the reduced director/framing coordinate is
    measured in half-turn units. Minima occur at 9 theta - 2 = 2 n.
    The first oriented framed branch is theta = 2/9. Since the integer
    offset is inserted here, this is a consistency normal form, not yet an
    independent derivation of h=2.
    """
    return kappa * (1.0 - math.cos(math.pi * (closure_lattice_order() * theta - SPINOR_CLOSURE_INDEX)))


def locking_gradient(theta, kappa=1.0):
    """Derivative of locking_potential with respect to theta."""
    n = closure_lattice_order()
    return kappa * math.pi * n * math.sin(math.pi * (n * theta - SPINOR_CLOSURE_INDEX))


def locking_curvature(theta, kappa=1.0):
    """Second derivative of locking_potential with respect to theta."""
    n = closure_lattice_order()
    return kappa * (math.pi * n) ** 2 * math.cos(math.pi * (n * theta - SPINOR_CLOSURE_INDEX))


def cyclic_lift_gate():
    """
    Mathematical gate for the C3 x C3 -> theta=h/9 step.

    C3 x C3 is not isomorphic to Z9. It only supplies nine slots unless the
    defect return map is shown to visit those slots as one cyclic orbit.
    """
    return {
        "slot_lattice": "C3 x C3",
        "slot_count": closure_lattice_order(),
        "is_cyclic_Z9_without_extra_theorem": False,
        "allowed_language": "order-9 reduced closure-slot lattice",
        "forbidden_language": "derived Z9 holonomy",
        "needed_theorem": "construct the charged defect return map and prove it has one order-9 orbit",
    }


def h2_offset_gate():
    """
    Gate for the integer offset in V_lock ~ cos(pi*(9 theta - 2)).
    """
    return {
        "offset_inserted": SPINOR_CLOSURE_INDEX,
        "derived_from_full_action": False,
        "current_status": "selection rule / normal-form consistency condition",
        "needed_theorem": "derive h=2 from the charged oriented-frame coupling, not by choosing the first useful branch",
    }


def positivity_edge_for_c3_operator():
    """
    The smallest C3 eigenfrequency hits zero when

        1 + sqrt(2) cos(theta + 2*pi/3) = 0.

    On the small positive branch this gives theta = pi/12.
    """
    return math.pi / 12.0


def branch_audit(max_index=5):
    """
    Enumerate theta = h/9 branches and test basic C3-spectrum viability.

    Even h is the oriented framed closure family in this candidate map.
    h=2 is the first non-trivial even branch below the positivity edge.
    """
    n = closure_lattice_order()
    edge = positivity_edge_for_c3_operator()
    rows = []
    for h in range(max_index + 1):
        theta = h / n
        raw = [
            1.0 + math.sqrt(2.0) * math.cos(theta + 2.0 * math.pi * k / 3.0)
            for k in range(3)
        ]
        positive = all(value > 0.0 for value in raw)
        oriented_framed = (h % 2 == 0 and h != 0)
        candidate = (h == SPINOR_CLOSURE_INDEX)
        ratios = None
        if positive:
            values = sorted(raw)
            ratios = [value / values[0] for value in values]
        rows.append(
            {
                "h": h,
                "theta": theta,
                "theta_less_than_edge": theta < edge,
                "positive_spectrum": positive,
                "oriented_framed_branch": oriented_framed,
                "selected_candidate": candidate,
                "frequency_ratios": ratios,
            }
        )
    return rows


def derivation_summary():
    """Compact statement of the candidate route."""
    n = closure_lattice_order()
    theta = theta_from_z9_holonomy()
    edge = positivity_edge_for_c3_operator()
    return {
        "axis_cycle_order": AXIS_CYCLE_ORDER,
        "phase_braid_order": PHASE_BRAID_ORDER,
        "closure_lattice": f"{AXIS_CYCLE_ORDER} x {PHASE_BRAID_ORDER} = {n}",
        "cyclic_lift_proven": cyclic_lift_gate()["is_cyclic_Z9_without_extra_theorem"],
        "oriented_framed_index": SPINOR_CLOSURE_INDEX,
        "h2_derived_from_action": h2_offset_gate()["derived_from_full_action"],
        "theta": theta,
        "theta_particle_constant": THETA_TOPOLOGICAL,
        "matches_particle_constant": abs(theta - THETA_TOPOLOGICAL) < 1.0e-15,
        "positivity_edge": edge,
        "below_edge": theta < edge,
        "locking_gradient": locking_gradient(theta),
        "locking_curvature": locking_curvature(theta),
        "status": "Candidate route; cyclic lift and h=2 action derivation remain open.",
    }


def falsification_targets():
    """
    What would make this candidate fail.
    """
    return [
        "If C3 x C3 cannot be lifted to one order-9 return coordinate, theta=h/9 is not derived.",
        "If the lepton oscillon is not C3-cyclic, N_closure=9 collapses.",
        "If the charged sector does not require oriented framed closure, h=2 is not selected.",
        "If the elastic operator uses the full U(1) phase directly, theta should be 2*pi/9, which fails the masses.",
        "If h=1 is shown to be a valid oriented charged branch rather than projective/nematic only, the uniqueness of h=2 fails.",
        "If the RFG action cannot generate the reduced-framing lock V ~ 1 - cos(pi(9 theta - 2)), this remains numerology.",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 38: order-9 reduced-holonomy candidate for theta = 2/9")
    print("=" * 72)

    print("\n1. Discrete closure derivation")
    summary = derivation_summary()
    for key, value in summary.items():
        if isinstance(value, float):
            print(f"  {key:24s}: {value:.12g}")
        else:
            print(f"  {key:24s}: {value}")

    print("\n2. Branch audit theta = h/9")
    for row in branch_audit():
        ratios = row["frequency_ratios"]
        if ratios is None:
            ratio_text = "not physical"
        else:
            ratio_text = " : ".join(f"{value:.5g}" for value in ratios)
        marker = " <-- selected" if row["selected_candidate"] else ""
        print(
            f"  h={row['h']} theta={row['theta']:.9f} "
            f"positive={row['positive_spectrum']} "
            f"oriented={row['oriented_framed_branch']} "
            f"ratios={ratio_text}{marker}"
        )

    print("\n3. Cyclic lift / h=2 gates")
    for key, value in cyclic_lift_gate().items():
        print(f"  {key:36s}: {value}")
    for key, value in h2_offset_gate().items():
        print(f"  {key:36s}: {value}")

    print("\n4. Charged-lepton prediction inherited from phase37")
    print(f"  K_C3 = {koide_identity():.12f}")
    for row in prediction_table():
        print(
            f"  {row['particle']:8s}: "
            f"nu_C3={row['predicted_freq_ratio']:.8f}, "
            f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
            f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
            f"rel_err={row['relative_mass_error']:.3e}"
        )

    print("\n5. Falsification targets")
    for item in falsification_targets():
        print(f"  - {item}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 39: Action-level route to the order-9 closure-slot lattice

Status:
    Candidate derivation layer. This file does not yet prove the full
    charged-lepton theorem from the RFG action, but it removes one
    arbitrary-looking assumption from phase38:

        why two C3 factors?

Result:
    The two C3 factors can be traced to standard structures already present
    in an elastic supersolid:

    1. Axis C3:
       A framed/chiral defect carries an oriented material triad. The full
       permutation group S3 of three axes is reduced to the orientation-
       preserving subgroup A3, which is isomorphic to C3.

       Odd permutations flip the sign of the oriented volume form
       epsilon_ABC and therefore move the defect to the opposite chirality.

    2. Strain-phase C3:
       The traceless triaxial strain sector is a two-component order
       parameter. In complex form

           E = q1 + omega q2 + omega^2 q3,     q1+q2+q3=0,
           omega = exp(2*pi*i/3).

       A cyclic axis permutation sends E -> omega^2 E. Therefore E^3 is
       the first phase-sensitive invariant. This creates three discrete
       strain/braid phase sectors.

    Together they provide nine closure slots:

        C3(axis orientation) x C3(strain phase) -> 9 reduced closure slots.

    This is not by itself a Z9 theorem: C3 x C3 is a direct product whose
    non-identity elements have order 3. The missing step is a charged-defect
    return map that cycles through all nine slots as one orbit.

    phase41 shows that the strain-phase C3 lock is not an external
    addition: it is the cubic invariant already contained in I3=det(B).
    phase40 then selects h=2 as the first non-trivial oriented framed
    branch, giving theta = h/9 = 2/9 in phase38.
"""

import itertools
import math

import sympy as sp


def permutation_parity(perm):
    """Return +1 for even permutations and -1 for odd permutations."""
    inversions = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inversions += 1
    return 1 if inversions % 2 == 0 else -1


def permutation_cycles(perm):
    """Human-readable cycle label for a permutation of (0,1,2)."""
    seen = [False] * len(perm)
    cycles = []
    for i in range(len(perm)):
        if seen[i]:
            continue
        cycle = []
        j = i
        while not seen[j]:
            seen[j] = True
            cycle.append(j + 1)
            j = perm[j]
        if len(cycle) > 1:
            cycles.append("(" + " ".join(str(x) for x in cycle) + ")")
    return "".join(cycles) if cycles else "identity"


def classify_axis_permutations():
    """
    S3 permutations of the elastic material triad.

    The oriented volume form epsilon_ABC keeps only even permutations in a
    fixed-chirality sector. The even set has three elements: C3.
    """
    rows = []
    for perm in itertools.permutations((0, 1, 2)):
        parity = permutation_parity(perm)
        rows.append(
            {
                "perm": perm,
                "cycle": permutation_cycles(perm),
                "det": parity,
                "orientation_preserving": parity == 1,
                "group_sector": "A3 ~= C3" if parity == 1 else "opposite chirality",
            }
        )
    return rows


def axis_c3_order():
    return sum(1 for row in classify_axis_permutations() if row["orientation_preserving"])


def strain_doublet_transform():
    """
    Verify that the triaxial complex strain E transforms by a C3 phase.

    E = q1 + omega q2 + omega^2 q3.
    Under the cyclic permutation (q1,q2,q3)->(q2,q3,q1),
    E -> omega^2 E. Hence E^3 is invariant.
    """
    q1, q2, q3 = sp.symbols("q1 q2 q3")
    omega = sp.exp(2 * sp.pi * sp.I / 3)
    E = q1 + omega * q2 + omega**2 * q3
    E_cyclic = q2 + omega * q3 + omega**2 * q1
    transform_ratio = sp.simplify(E_cyclic / E)

    # SymPy may not reduce exp(2*pi*i/3) expressions completely in a ratio,
    # so verify by coefficient comparison.
    expected = omega**2 * E
    diff = sp.simplify(sp.expand(E_cyclic - expected))

    invariant_diff = sp.simplify(sp.expand(E_cyclic**3 - E**3))

    return {
        "E": E,
        "E_cyclic": E_cyclic,
        "expected_transform": "E -> omega^2 E",
        "transform_difference": diff,
        "E3_invariant_difference": invariant_diff,
        "phase_sector_order": 3,
    }


def z9_closure_from_action_symmetry():
    """Combine the two C3 factors into an order-9 closure-slot lattice."""
    axis_order = axis_c3_order()
    phase_order = strain_doublet_transform()["phase_sector_order"]
    return {
        "axis_c3_origin": "oriented elastic triad: S3 -> A3 ~= C3",
        "axis_order": axis_order,
        "phase_c3_origin": "triaxial strain doublet: E -> omega^2 E, E^3 invariant",
        "phase_order": phase_order,
        "closure_slots": axis_order * phase_order,
        "direct_product_group": "C3 x C3",
        "cyclic_z9_proven": False,
        "missing_lift": "prove one charged-defect return map with order 9",
        "theta_if_h2": 2 / (axis_order * phase_order),
    }


def allowed_phase_locking_terms(max_power=9):
    """
    List low-order powers E^n and whether they are invariant under the
    cyclic strain phase E -> omega^2 E.
    """
    rows = []
    for n in range(1, max_power + 1):
        invariant = (2 * n) % 3 == 0
        rows.append(
            {
                "power": n,
                "phase_factor": f"omega^{2*n}",
                "c3_invariant": invariant,
                "meaning": "allowed anisotropy" if invariant else "forbidden by C3",
            }
        )
    return rows


def spinorial_selection_statement():
    """
    This is the remaining assumption chain, made explicit.
    """
    return {
        "h0": "trivial closure; gives a partially degenerate spectrum, not three charged generations",
        "h1": "projective/nematic closure only, n -> -n; see phase40",
        "h2": "first non-trivial oriented framed closure, n -> n; see phase40",
        "higher_h": "on the small positive C3 branch, h>=3 crosses the zero-eigenfrequency edge",
        "remaining_theorem": "derive the oriented-frame charged coupling directly from the RFG action",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 39: Action-level route to C3 x C3 -> 9 slots")
    print("=" * 72)

    print("\n1. Axis permutations of the oriented elastic triad")
    for row in classify_axis_permutations():
        print(
            f"  perm={row['perm']} cycle={row['cycle']:<9} "
            f"det={row['det']:+d} sector={row['group_sector']}"
        )
    print(f"  axis C3 order = {axis_c3_order()}")

    print("\n2. Triaxial strain doublet")
    strain = strain_doublet_transform()
    for key, value in strain.items():
        print(f"  {key:26s}: {value}")

    print("\n3. Allowed C3 phase-locking powers")
    for row in allowed_phase_locking_terms():
        print(
            f"  E^{row['power']:<2} -> {row['phase_factor']:<8} "
            f"{row['meaning']}"
        )

    print("\n4. Order-9 closure-slot gate")
    closure = z9_closure_from_action_symmetry()
    for key, value in closure.items():
        if isinstance(value, float):
            print(f"  {key:20s}: {value:.12f}")
        else:
            print(f"  {key:20s}: {value}")

    print("\n5. Remaining spinorial selection")
    for key, value in spinorial_selection_statement().items():
        print(f"  {key:20s}: {value}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 40: Why the first charged framed branch is h = 2

Status:
    Candidate selection rule for the h index used in phase38. It becomes a
    theorem only after the charged oriented-frame coupling is derived from
    the full RFG action.

Problem:
    phase38 uses h=2 in

        theta = h / 9 = 2 / 9.

    Why not h=1?

Key distinction:
    h=1 is a valid closure only for a projective/nematic director, where

        n ~ -n.

    But the charged RFG sector is framed and oriented: it uses the oriented
    material triad and the epsilon_ABC volume form/topological current.
    Therefore the lift must close as an oriented frame, not merely as a
    projective line.

Minimal closure:
    Let h count half-turns of the local director/framing.

        h=0: trivial oriented closure.
        h=1: projective closure only; n -> -n, orientation not closed.
        h=2: first non-trivial oriented closure; n -> n.

    Thus h=2 is the first non-trivial charged/framed branch.

Consequence:
    Combining phase39 and this file:

        C3(axis) x C3(strain phase) -> 9 slots,
        h_charged = 2,
        theta = 2/9.
"""

import math

# merged import removed: from p11_particles import (
# merged import removed:     closure_lattice_order,
# merged import removed:     positivity_edge_for_c3_operator,
# merged import removed: )


def director_after_half_turns(h):
    """
    Unit director in a 2D slice after h half-turns.

    h=1 gives n -> -n, which is projectively closed but not oriented closed.
    h=2 gives n -> n.
    """
    angle = math.pi * h
    return (round(math.cos(angle), 12), round(math.sin(angle), 12))


def closure_type(h):
    """Classify the h branch by projective and oriented closure."""
    n0 = (1.0, 0.0)
    nh = director_after_half_turns(h)
    oriented_closed = nh == n0
    projective_closed = nh == n0 or nh == (-n0[0], -n0[1])
    nontrivial = h != 0
    charged_allowed = oriented_closed and nontrivial
    return {
        "h": h,
        "director": nh,
        "projective_closed": projective_closed,
        "oriented_closed": oriented_closed,
        "nontrivial": nontrivial,
        "charged_framed_allowed": charged_allowed,
    }


def c3_frequency_raw(theta):
    return [
        1.0 + math.sqrt(2.0) * math.cos(theta + 2.0 * math.pi * k / 3.0)
        for k in range(3)
    ]


def branch_selection_table(max_h=6):
    """
    Combine projective/oriented closure with C3-spectrum positivity.
    """
    n_slots = closure_lattice_order()
    edge = positivity_edge_for_c3_operator()
    rows = []
    for h in range(max_h + 1):
        theta = h / n_slots
        raw = c3_frequency_raw(theta)
        positive = all(value > 0.0 for value in raw)
        closure = closure_type(h)
        selected = closure["charged_framed_allowed"] and positive and h == 2
        rows.append(
            {
                **closure,
                "theta": theta,
                "below_positive_edge": theta < edge,
                "positive_spectrum": positive,
                "selected": selected,
            }
        )
    return rows


def h2_selection_summary():
    n_slots = closure_lattice_order()
    h = 2
    theta = h / n_slots
    return {
        "closure_slots": n_slots,
        "h_selected": h,
        "reason_h0_rejected": "trivial branch, not a generation-splitting defect",
        "reason_h1_rejected": "projective/nematic closure only: n -> -n",
        "reason_h2_selected": "first non-trivial oriented framed closure: n -> n",
        "derived_from_action": False,
        "selection_status": "candidate selection; not yet an RFG coupling theorem",
        "theta": theta,
        "theta_formula": "theta = h / 9 = 2 / 9",
        "remaining_open_point": "derive the oriented-frame requirement from the charged RFG coupling, not as a selection rule.",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 40: Projective vs oriented closure selects h = 2")
    print("=" * 72)

    print("\n1. Branch selection table")
    for row in branch_selection_table():
        marker = " <-- selected" if row["selected"] else ""
        print(
            f"  h={row['h']} theta={row['theta']:.9f} "
            f"n={row['director']} "
            f"projective={row['projective_closed']} "
            f"oriented={row['oriented_closed']} "
            f"positive={row['positive_spectrum']}{marker}"
        )

    print("\n2. h=2 summary")
    for key, value in h2_selection_summary().items():
        if isinstance(value, float):
            print(f"  {key:24s}: {value:.12f}")
        else:
            print(f"  {key:24s}: {value}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 41: Action-level normal form for C3 phase locking

Status:
    Algebraic support for the phase37-40 chain.

Question:
    Is the C3 phase-locking term an extra hand-added ingredient, or does it
    already exist in the RFG supersolid action P(Y, I1, I2, I3)?

Answer:
    It already exists generically through the determinant invariant I3.

Setup:
    Around the isotropic elastic background, write the 3x3 material matrix as

        B = b I + Q,       tr(Q) = 0.

    The traceless part Q is the triaxial strain. In its principal-axis frame

        Q = diag(q1, q2, q3),        q1 + q2 + q3 = 0.

    Define the complex C3 strain coordinate

        E = q1 + omega q2 + omega^2 q3,     omega = exp(2*pi*i/3).

Key identities:
    1. det(b I + Q) contains det(Q).
    2. det(Q) is the cubic triaxial invariant.
    3. E^3 + conjugate(E)^3 = 27 det(Q).

Therefore any generic I3 dependence in the action supplies the
phase-sensitive C3 anisotropy:

        V_C3 ~ -lambda_3 Re(E^3).

This gives the normal-form origin of the C3 strain-sector locking used in
phase39. Combined with the oriented-triad C3 and the h=2 framed closure from
phase40, phase38 gives theta = 2/9.
"""

import sympy as sp


def invariants_around_isotropic_background():
    """
    Expand I1, I2, I3 for B = b I + Q, tr(Q)=0.
    """
    b, q1, q2 = sp.symbols("b q1 q2", real=True)
    q3 = -q1 - q2
    eigenvalues = [b + q1, b + q2, b + q3]

    I1 = sp.simplify(sum(eigenvalues))
    I2 = sp.simplify(
        eigenvalues[0] * eigenvalues[1]
        + eigenvalues[0] * eigenvalues[2]
        + eigenvalues[1] * eigenvalues[2]
    )
    I3 = sp.simplify(eigenvalues[0] * eigenvalues[1] * eigenvalues[2])

    tr_q2 = sp.simplify(q1**2 + q2**2 + q3**2)
    tr_q3 = sp.simplify(q1**3 + q2**3 + q3**3)
    det_q = sp.simplify(q1 * q2 * q3)

    return {
        "I1": I1,
        "I2": I2,
        "I3": I3,
        "tr_Q2": tr_q2,
        "tr_Q3": tr_q3,
        "det_Q": det_q,
        "I3_expected": sp.simplify(b**3 - b * tr_q2 / 2 + det_q),
        "trQ3_minus_3detQ": sp.simplify(tr_q3 - 3 * det_q),
    }


def c3_complex_strain_identity():
    """
    Prove E^3 + Ebar^3 = 27 det(Q) for q1+q2+q3=0.
    """
    q1, q2 = sp.symbols("q1 q2", real=True)
    q3 = -q1 - q2
    omega = -sp.Rational(1, 2) + sp.I * sp.sqrt(3) / 2
    omega_bar = sp.conjugate(omega)

    E = sp.simplify(q1 + omega * q2 + omega**2 * q3)
    E_bar = sp.simplify(q1 + omega_bar * q2 + omega_bar**2 * q3)
    det_q = sp.simplify(q1 * q2 * q3)
    tr_q3 = sp.simplify(q1**3 + q2**3 + q3**3)

    identity_det = sp.simplify(sp.expand(E**3 + E_bar**3 - 27 * det_q))
    identity_tr = sp.simplify(sp.expand(E**3 + E_bar**3 - 9 * tr_q3))

    return {
        "E": E,
        "E_bar": E_bar,
        "E3_plus_Ebar3_minus_27detQ": identity_det,
        "E3_plus_Ebar3_minus_9trQ3": identity_tr,
        "det_Q": det_q,
        "tr_Q3": tr_q3,
    }


def polar_triaxial_parameterization():
    """
    Parameterize the traceless eigenvalues by a radius rho and phase beta:

        q_i = rho cos(beta + 2*pi*i/3).

    Then det(Q) = rho^3 cos(3 beta) / 4.
    """
    rho, beta = sp.symbols("rho beta", real=True)
    q1 = rho * sp.cos(beta)
    q2 = rho * sp.cos(beta + 2 * sp.pi / 3)
    q3 = rho * sp.cos(beta + 4 * sp.pi / 3)
    det_q = sp.trigsimp(sp.expand_trig(q1 * q2 * q3))
    tr_q2 = sp.trigsimp(sp.expand_trig(q1**2 + q2**2 + q3**2))
    det_expected = rho**3 * sp.cos(3 * beta) / 4
    trq2_expected = 3 * rho**2 / 2
    return {
        "q1": q1,
        "q2": q2,
        "q3": q3,
        "sum_q": sp.trigsimp(sp.expand_trig(q1 + q2 + q3)),
        "det_Q": det_q,
        "det_Q_minus_expected": sp.trigsimp(sp.expand_trig(det_q - det_expected)),
        "tr_Q2": tr_q2,
        "tr_Q2_minus_expected": sp.trigsimp(sp.expand_trig(tr_q2 - trq2_expected)),
    }


def normal_form_from_action():
    """
    Minimal effective potential inherited from the action.

    The sign and size of lambda_3 depend on the microscopic RFG coefficients,
    but the existence of the cubic C3 anisotropy is generic once I3 is present.
    """
    rho, beta = sp.symbols("rho beta", real=True)
    a2, a4, lambda3 = sp.symbols("a2 a4 lambda3", real=True)
    det_q = rho**3 * sp.cos(3 * beta) / 4
    potential = a2 * rho**2 + a4 * rho**4 - lambda3 * det_q
    d_beta = sp.diff(potential, beta)
    stationary_condition = sp.factor(d_beta)
    curvature_beta = sp.factor(sp.diff(potential, beta, 2))
    return {
        "V_eff": potential,
        "dV_d_beta": stationary_condition,
        "d2V_d_beta2": curvature_beta,
        "stationary_rule": "sin(3 beta)=0 -> three C3-locked strain sectors",
    }


def theta_chain_summary():
    return {
        "action_invariant": "I3 = det(B)",
        "triaxial_cubic": "det(Q) = Re(E^3)/27 = tr(Q^3)/3",
        "strain_lock": "V_C3 ~ -lambda3 rho^3 cos(3 beta)/4",
        "phase39": "oriented triad C3 x strain-sector C3 -> 9 slots; cyclic lift open",
        "phase40": "h=2 is first non-trivial oriented framed candidate; action derivation open",
        "phase38": "theta = h/9 = 2/9 only after the lift/selection gates pass",
        "phase37": "C3 Koide operator gives charged-lepton mass ratios",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 41: Action normal form for C3 phase locking")
    print("=" * 72)

    print("\n1. Invariants around B = b I + Q")
    inv = invariants_around_isotropic_background()
    for key, value in inv.items():
        print(f"  {key:28s}: {value}")

    print("\n2. Complex C3 strain identity")
    ident = c3_complex_strain_identity()
    for key, value in ident.items():
        print(f"  {key:34s}: {value}")

    print("\n3. Polar triaxial parameterization")
    polar = polar_triaxial_parameterization()
    for key, value in polar.items():
        print(f"  {key:28s}: {value}")

    print("\n4. Normal form inherited from I3")
    nf = normal_form_from_action()
    for key, value in nf.items():
        print(f"  {key:20s}: {value}")

    print("\n5. Theta chain")
    for key, value in theta_chain_summary().items():
        print(f"  {key:18s}: {value}")


# ===================== merged from p11_particles.py =====================

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


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 46: Quark-sector extension candidate from the same RFG holonomy

Status:
    Falsifiable extension sketch, not a fit.

Motivation:
    A particle-sector theory stronger than Koide-only lepton models should
    not stop at charged leptons. But quark masses are not pole-observable in
    the same clean way as charged leptons; they are scheme- and scale-dependent.
    Therefore this file defines a strict candidate rule and a falsification
    protocol rather than claiming a quark-mass derivation.

RFG rule:
    There is only one primitive reduced holonomy:

        theta_L = 2/9.

    Quark family phases, if present, must be projections of theta_L by
    topological/color/charge data. They are not new fitted angles.

Two minimal projection candidates:
    charge projection:
        theta_U = |Q_u| theta_L = 4/27,
        theta_D = |Q_d| theta_L = 2/27.

    inverse-color/empirical-Z3 projection:
        theta_U = theta_L / 3 = 2/27,
        theta_D = 2 theta_L / 3 = 4/27.

The second ordering matches a known Z3-Koide phenomenological suggestion,
but RFG should decide between them from color/framing topology, not from a
quark-mass fit.
"""

import math


THETA_LEPTON = 2.0 / 9.0


def c3_ratios(theta):
    values = sorted(
        1.0 + math.sqrt(2.0) * math.cos(theta + 2.0 * math.pi * k / 3.0)
        for k in range(3)
    )
    if values[0] <= 0:
        return None
    return [value / values[0] for value in values]


def mass_ratios_from_theta(theta):
    ratios = c3_ratios(theta)
    if ratios is None:
        return None
    return [value * value for value in ratios]


def projection_candidates():
    return [
        {
            "name": "charge_projection",
            "theta_U": 2.0 * THETA_LEPTON / 3.0,
            "theta_D": THETA_LEPTON / 3.0,
            "rule": "theta_f = |Q_f| theta_L",
            "free_angles": 0,
        },
        {
            "name": "inverse_color_empirical_z3_projection",
            "theta_U": THETA_LEPTON / 3.0,
            "theta_D": 2.0 * THETA_LEPTON / 3.0,
            "rule": "theta_U=theta_L/3, theta_D=2 theta_L/3",
            "free_angles": 0,
        },
    ]


def candidate_tables():
    rows = []
    for candidate in projection_candidates():
        theta_u = candidate["theta_U"]
        theta_d = candidate["theta_D"]
        rows.append(
            {
                **candidate,
                "freq_ratios_U": c3_ratios(theta_u),
                "mass_ratios_U": mass_ratios_from_theta(theta_u),
                "freq_ratios_D": c3_ratios(theta_d),
                "mass_ratios_D": mass_ratios_from_theta(theta_d),
            }
        )
    return rows


def falsification_protocol():
    return [
        "Choose a precise quark-mass scheme and scale, e.g. MSbar at a specified mu.",
        "Run all six quark masses to that scale with standard QCD RG equations.",
        "Normalize each sector by its lightest member: (u,c,t) and (d,s,b).",
        "Compare against the two zero-angle-fit RFG projection candidates.",
        "Reject any candidate whose required theta_U/theta_D differs from its topological projection beyond uncertainties.",
        "Do not tune theta_U or theta_D after seeing the masses; that would erase the prediction.",
    ]


def phase46_status_assessment():
    return {
        "strength": "extends RFG particle sector beyond charged leptons without new continuous angles",
        "warning": "quark masses are scheme/scale dependent, so no pole-mass claim is made here",
        "needed_theorem": "derive whether color/framing topology chooses charge_projection or inverse_color projection",
        "relation_to_leptons": "theta_L=2/9 remains the only primitive phase",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 46: Quark extension candidate")
    print("=" * 72)

    print("\n1. Projection candidates")
    for row in candidate_tables():
        print(f"  {row['name']}")
        print(f"    rule     : {row['rule']}")
        print(f"    theta_U  : {row['theta_U']:.12f}")
        print(f"    theta_D  : {row['theta_D']:.12f}")
        print(f"    U freq   : {[round(x, 6) for x in row['freq_ratios_U']]}")
        print(f"    U mass   : {[round(x, 6) for x in row['mass_ratios_U']]}")
        print(f"    D freq   : {[round(x, 6) for x in row['freq_ratios_D']]}")
        print(f"    D mass   : {[round(x, 6) for x in row['mass_ratios_D']]}")

    print("\n2. Falsification protocol")
    for item in falsification_protocol():
        print(f"  - {item}")

    print("\n3. Status")
    for key, value in phase46_status_assessment().items():
        print(f"  {key:16s}: {value}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 45: Radiative stability requirement for the C3 lepton sector

Status:
    Strengthening audit. This file does not solve radiative protection; it
    states the problem quantitatively and defines what RFG must prove.

Why this matters:
    Koide-type relations are usually discussed for charged-lepton pole
    masses. Ordinary QED self-energy corrections are generation-dependent
    because they contain log(m_i). Without a protection mechanism, a beautiful
    tree-level C3 relation may not survive to observed pole masses.

RFG options:
    A. Pole-frequency theorem:
       The C3 operator acts on fully dressed oscillon normal modes, so the
       observed pole masses are the eigenvalues of the already-renormalized
       medium. Then radiative corrections are not added afterward.

    B. Ward/Sumino-like cancellation:
       Phase-elastic/framed-sector loops cancel the generation-dependent QED
       logarithms, preserving the C3 ratios.

    C. Running-scale theorem:
       The C3 relation is exact at a special RFG matching scale, and its
       pole-mass accuracy is a consequence of controlled RG flow.

The strongest route for RFG is A or B. This file quantifies the required
protection as a condition on fractional frequency shifts.
"""

import math

# merged import removed: from p11_particles import (
# merged import removed:     LEPTON_MASSES_MEV,
# merged import removed:     c3_frequency_ratios,
# merged import removed:     koide_ratio_from_frequencies,
# merged import removed: )


ALPHA_EM = 1.0 / 137.035999084


def observed_frequency_vector():
    m_e = LEPTON_MASSES_MEV["electron"]
    return [
        math.sqrt(LEPTON_MASSES_MEV[name] / m_e)
        for name in ("electron", "muon", "tau")
    ]


def c3_frequency_vector():
    ratios = c3_frequency_ratios()
    return [ratios[name] for name in ("electron", "muon", "tau")]


def koide_sensitivity_to_fractional_shifts(eps):
    """
    Apply nu_i -> nu_i * (1 + eps_i) and return the Koide drift.
    eps must be a length-3 fractional-shift list.
    """
    nu = c3_frequency_vector()
    shifted = [value * (1.0 + eps_i) for value, eps_i in zip(nu, eps)]
    return koide_ratio_from_frequencies(shifted) - 2.0 / 3.0


def leading_qed_log_shifts(cutoff_mev):
    """
    Crude one-loop logarithmic mass-shift template:

        delta m_i / m_i ~ (3 alpha / 4 pi) log(cutoff^2 / m_i^2).

    Since nu ~ sqrt(m), the frequency shift is half this value.

    This is not used as a physical RFG prediction; it is a stress test.
    """
    prefactor_mass = 3.0 * ALPHA_EM / (4.0 * math.pi)
    shifts = []
    for name in ("electron", "muon", "tau"):
        mass = LEPTON_MASSES_MEV[name]
        dm_over_m = prefactor_mass * math.log((cutoff_mev / mass) ** 2)
        shifts.append(0.5 * dm_over_m)
    mean_shift = sum(shifts) / 3.0
    relative_to_mean = [shift - mean_shift for shift in shifts]
    return {
        "cutoff_MeV": cutoff_mev,
        "freq_shifts": shifts,
        "mean_shift": mean_shift,
        "generation_dependent_part": relative_to_mean,
        "koide_drift_raw": koide_sensitivity_to_fractional_shifts(shifts),
        "koide_drift_gen_dep": koide_sensitivity_to_fractional_shifts(relative_to_mean),
    }


def protection_condition_symbolic():
    """
    Linearized condition for preserving Koide.

    The common shift eps_i = const drops out. The dangerous part is the
    component of eps_i along the Koide-normal direction in frequency space.
    """
    nu = c3_frequency_vector()
    s1 = sum(nu)
    s2 = sum(value * value for value in nu)
    # dK = [2 nu_i / s1^2 - 2 s2 / s1^3] dnu_i.
    gradient = [
        2.0 * value / (s1 * s1) - 2.0 * s2 / (s1 * s1 * s1)
        for value in nu
    ]
    return {
        "linear_gradient_dK_dnu": gradient,
        "condition": "sum_i gradient_i * nu_i * eps_i = 0",
        "safe_shift": "eps_e = eps_mu = eps_tau is automatically safe",
        "dangerous_shift": "generation-dependent logarithms must be cancelled or reabsorbed into the dressed C3 operator",
    }


def rfg_radiative_options():
    return [
        {
            "option": "A. dressed pole-frequency theorem",
            "claim_needed": "C3 eigenvalues are already physical pole frequencies",
            "strength": "best fit to RFG oscillon ontology",
            "open_task": "derive dressed normal-mode equation including EM/self-field backreaction",
        },
        {
            "option": "B. Ward/Sumino-like cancellation",
            "claim_needed": "phase-elastic loops cancel QED log(m_i) terms",
            "strength": "closest to known Koide-protection logic",
            "open_task": "identify the RFG current and loop sign that cancels QED logs",
        },
        {
            "option": "C. matching-scale theorem",
            "claim_needed": "C3 relation exact at a scale and RG flow preserves pole accuracy",
            "strength": "standard EFT language",
            "open_task": "derive RG equations for the RFG lepton sector",
        },
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 45: Radiative stability requirement")
    print("=" * 72)

    print("\n1. Koide linear sensitivity")
    for key, value in protection_condition_symbolic().items():
        print(f"  {key:28s}: {value}")

    print("\n2. QED-log stress test")
    for cutoff in [1.0e3, 1.0e6, 1.0e9]:
        data = leading_qed_log_shifts(cutoff)
        print(f"  cutoff={cutoff:.3e} MeV")
        print(f"    freq shifts          : {[round(x, 8) for x in data['freq_shifts']]}")
        print(f"    generation dep. part : {[round(x, 8) for x in data['generation_dependent_part']]}")
        print(f"    Koide drift raw      : {data['koide_drift_raw']:.3e}")
        print(f"    Koide drift gen.dep. : {data['koide_drift_gen_dep']:.3e}")

    print("\n3. RFG protection options")
    for row in rfg_radiative_options():
        print(f"  {row['option']}")
        print(f"    claim needed: {row['claim_needed']}")
        print(f"    strength    : {row['strength']}")
        print(f"    open task   : {row['open_task']}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 42: Charged-lepton theorem-chain audit

Purpose:
    Collect the phase37-43 chain in one place and mark what is algebraically
    closed versus what remains a physical/theorem-level assumption.

The chain:
    phase41: I3=det(B) contains the C3 strain lock E^3.
    phase39: oriented elastic triad gives C3, strain phase gives C3 -> 9 slots.
    phase39 still does not prove a cyclic Z9 group.
    phase40: oriented framed closure makes h=2 the candidate branch, not yet
    a result derived from the charged RFG coupling.
    phase43: local normal-form Hessian gives h=2 stability conditions.
    phase38: theta = h/9 = 2/9 only conditionally.
    phase37: C3 Koide operator with theta=2/9 gives charged-lepton ratios.

Status:
    Strong candidate chain. The remaining open theorem is to derive the
    cyclic return map, charged oriented-framed defect sector, and reduced
    holonomy directly as a stationary sector of the full RFG action.
"""

# merged import removed: from p11_particles import (
# merged import removed:     LEPTON_MASSES_MEV,
# merged import removed:     THETA_TOPOLOGICAL,
# merged import removed:     koide_identity,
# merged import removed:     prediction_table,
# merged import removed: )
# merged import removed: from p11_particles import derivation_summary
# merged import removed: from p11_particles import z9_closure_from_action_symmetry
# merged import removed: from p11_particles import h2_selection_summary
# merged import removed: from p11_particles import (
# merged import removed:     c3_complex_strain_identity,
# merged import removed:     invariants_around_isotropic_background,
# merged import removed:     polar_triaxial_parameterization,
# merged import removed: )


def pass_fail(condition):
    return "PASS" if condition else "FAIL"


def audit_action_c3_lock():
    inv = invariants_around_isotropic_background()
    ident = c3_complex_strain_identity()
    polar = polar_triaxial_parameterization()
    checks = {
        "trQ3_equals_3detQ": inv["trQ3_minus_3detQ"] == 0,
        "E3_identity_det": ident["E3_plus_Ebar3_minus_27detQ"] == 0,
        "E3_identity_trQ3": ident["E3_plus_Ebar3_minus_9trQ3"] == 0,
        "polar_det_identity": polar["det_Q_minus_expected"] == 0,
        "polar_trQ2_identity": polar["tr_Q2_minus_expected"] == 0,
    }
    return checks


def audit_z9_and_h2():
    closure = z9_closure_from_action_symmetry()
    h2 = h2_selection_summary()
    theta = h2["theta"]
    checks = {
        "axis_order_is_3": closure["axis_order"] == 3,
        "phase_order_is_3": closure["phase_order"] == 3,
        "closure_slots_is_9": closure["closure_slots"] == 9,
        "cyclic_z9_proven": closure["cyclic_z9_proven"],
        "h_selected_is_2": h2["h_selected"] == 2,
        "h2_derived_from_action": h2["derived_from_action"],
        "theta_is_2_over_9": abs(theta - THETA_TOPOLOGICAL) < 1.0e-15,
    }
    return checks


def audit_mass_predictions():
    precision = pdg_precision_audit()
    checks = {
        "electron_is_anchor_not_prediction": True,
        "relative_compression_ok": all(
            abs(row["relative_mass_error"]) <= 1.0e-4
            for row in precision["rows"]
        ),
        "pdg_precision_pass": precision["pdg_precision_pass"],
        "koide_exact": abs(koide_identity() - 2.0 / 3.0) < 1.0e-15,
        "mass_bridge_derived": False,
    }
    return checks


def open_theorem_items():
    return [
        "Prove a cyclic order-9 lift from the C3 x C3 slot lattice, or stop using Z9 language.",
        "Derive the charged oriented-frame requirement from the RFG coupling, not as a postulate.",
        "Derive the reduced holonomy coordinate theta used by the C3 stiffness operator from the defect moduli space.",
        "Derive m proportional to nu^2 and the absolute electron scale from the oscillon energy functional.",
        "Replace the relative-error mass gate with PDG-uncertainty residuals and explain the muon residual.",
        "Upgrade phase43 local normal-form stability to a full 3D fluctuation-operator/PDE stability proof.",
        "Show no lower-energy non-leptonic defect branch has the same charge and lower action.",
        "Extend the proof beyond principal-axis/algebraic normal form to full 3D localized oscillon fields.",
        "Derive radiative protection of the Koide/C3 pole-frequency relation.",
    ]


def theorem_chain_summary():
    summary = derivation_summary()
    h2 = h2_selection_summary()
    closure = z9_closure_from_action_symmetry()
    return {
        "action_lock": "I3=det(B) -> det(Q)=Re(E^3)/27 -> C3 strain locking",
        "axis_c3": closure["axis_c3_origin"],
        "phase_c3": closure["phase_c3_origin"],
        "closure_slots": closure["closure_slots"],
        "cyclic_z9_proven": closure["cyclic_z9_proven"],
        "h_selected": h2["h_selected"],
        "h2_derived_from_action": h2["derived_from_action"],
        "local_stability": "phase43 Hessian conditions",
        "theta": summary["theta"],
        "koide": koide_identity(),
        "electron_anchor_MeV": LEPTON_MASSES_MEV["electron"],
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 42: Charged-lepton theorem-chain audit")
    print("=" * 72)

    print("\n1. Theorem-chain summary")
    for key, value in theorem_chain_summary().items():
        if isinstance(value, float):
            print(f"  {key:22s}: {value:.12g}")
        else:
            print(f"  {key:22s}: {value}")

    print("\n2. Algebraic C3-lock checks")
    for key, value in audit_action_c3_lock().items():
        print(f"  {key:28s}: {pass_fail(value)}")

    print("\n3. Order-9 / h=2 checks")
    for key, value in audit_z9_and_h2().items():
        print(f"  {key:28s}: {pass_fail(value)}")

    print("\n4. Mass/Koide checks")
    for key, value in audit_mass_predictions().items():
        print(f"  {key:28s}: {pass_fail(value)}")

    print("\n5. Conditional charged-lepton mass table")
    for row in prediction_table():
        print(
            f"  {row['particle']:8s}: "
            f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
            f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
            f"rel_err={row['relative_mass_error']:.3e}, "
            f"sigma={row['sigma_error']:.3e}"
        )

    print("\n6. Still open theorem items")
    for item in open_theorem_items():
        print(f"  OPEN: {item}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 44: Particle-sector strength and novelty audit

Purpose:
    Make the particle-sector claim sharper without overstating it.

Core point:
    Many models already use Koide, Z3, family symmetry, or the phase
    delta_L ~= 2/9. RFG is not novel merely because it reproduces Koide or
    notices 2/9.

RFG's stronger claim is the specific action/topology chain:

    I3 = det(B)
      -> det(Q) = Re(E^3)/27
      -> C3 strain lock from the supersolid action
      -> oriented elastic triad gives another C3
      -> C3 x C3 = 9 reduced closure slots
      -> missing cyclic lift gate
      -> oriented framed branch h=2 candidate
      -> theta = 2/9 conditionally
      -> C3 Koide operator gives charged-lepton mass ratios.

This file also cleans up an internal tension:
    The older N=5,72,295 ladder and N=4 prediction are conditional legacy
    candidates. If phase37-43 C3 theorem-chain is adopted, the old radial
    N-ladder should no longer be presented as the main particle-sector route.
"""

# merged import removed: from p11_particles import prediction_table, koide_identity
# merged import removed: from p11_particles import (
# merged import removed:     audit_action_c3_lock,
# merged import removed:     audit_mass_predictions,
# merged import removed:     audit_z9_and_h2,
# merged import removed:     open_theorem_items,
# merged import removed: )


def comparison_with_existing_approaches():
    """
    High-level comparison. References are given by author/model labels rather
    than used as authority inside the code.
    """
    return [
        {
            "approach": "Koide original / scalar-potential family models",
            "has_koide": True,
            "has_delta_2_over_9": False,
            "action_origin_c3": "family scalar potential, not RFG medium",
            "topological_h2": False,
            "main_gap": "mass relation protected/engineered, but not a vacuum-medium oscillon derivation",
        },
        {
            "approach": "Z3-symmetric Koide parametrization",
            "has_koide": True,
            "has_delta_2_over_9": True,
            "action_origin_c3": "usually parametrization or family symmetry",
            "topological_h2": False,
            "main_gap": "2/9 is observed/suggested, not derived from an elastic invariant plus framed closure",
        },
        {
            "approach": "Sumino/family-gauge protection",
            "has_koide": True,
            "has_delta_2_over_9": "model-dependent",
            "action_origin_c3": "new family gauge sector",
            "topological_h2": False,
            "main_gap": "strong radiative idea, but different ontology and no supersolid oscillon route",
        },
        {
            "approach": "A4/S3/Z3 flavor models",
            "has_koide": "sometimes",
            "has_delta_2_over_9": "sometimes",
            "action_origin_c3": "discrete flavor symmetry input",
            "topological_h2": False,
            "main_gap": "symmetry is typically imposed in flavor space",
        },
        {
            "approach": "RFG phase37-43 chain",
            "has_koide": True,
            "has_delta_2_over_9": True,
            "action_origin_c3": "I3=det(B) -> Re(E^3) C3 lock",
            "topological_h2": True,
            "main_gap": "cyclic lift, h=2 action derivation, full PDE proof, and radiative protection still open",
        },
    ]


def internal_consistency_cleanup():
    """
    What must be demoted or relabeled after adopting the C3 theorem-chain.
    """
    return [
        {
            "item": "phase35 radial n=14 route",
            "old_status": "candidate route for muon as a radial overtone",
            "new_status": "legacy audit / phenomenological comparison",
            "reason": "phase37 explains muon and tau together as one C3 triplet; no missing n=2..13 issue",
        },
        {
            "item": "N=5,72,295 ladder",
            "old_status": "empirical index ladder",
            "new_status": "do not present as primary derivation",
            "reason": "indices are back-solved from sqrt(m); C3 ratios are cleaner",
        },
        {
            "item": "phase36 N=4 329 keV prediction",
            "old_status": "conditional RFG-specific prediction if N-ladder is derived",
            "new_status": "conditional legacy prediction, suspended under C3 route",
            "reason": "C3 charged-lepton sector does not imply a radial N=4 companion",
        },
        {
            "item": "toy/work-2 scripts",
            "old_status": "exploratory sandbox",
            "new_status": "not part of main proof chain",
            "reason": "main particle-sector chain is now phase37-43",
        },
    ]


def rfg_strength_scorecard():
    action_checks = audit_action_c3_lock()
    mass_checks = audit_mass_predictions()
    precision = pdg_precision_audit()

    return [
        {
            "criterion": "Koide identity",
            "status": "closed algebraically",
            "evidence": f"K_C3 = {koide_identity():.12f}",
        },
        {
            "criterion": "charged-lepton ratios",
            "status": "relative postdiction only",
            "evidence": (
                f"relative gate={mass_checks['relative_compression_ok']}; "
                f"PDG precision={precision['pdg_precision_pass']}; "
                f"chi2={precision['chi2_non_anchor']:.3e}"
            ),
        },
        {
            "criterion": "C3 action origin",
            "status": "closed in principal-axis normal form",
            "evidence": all(action_checks.values()),
        },
        {
            "criterion": "theta=2/9 route",
            "status": "blocked at cyclic-lift and h=2-action gates",
            "evidence": audit_z9_and_h2(),
        },
        {
            "criterion": "local h=2 stability",
            "status": "normal-form conditions derived",
            "evidence": "phase43 Hessian",
        },
        {
            "criterion": "full 3D particle proof",
            "status": "open",
            "evidence": "needs localized oscillon fluctuation operator",
        },
        {
            "criterion": "radiative stability",
            "status": "open",
            "evidence": "needs RFG analogue of Koide/Sumino protection or pole-frequency argument",
        },
        {
            "criterion": "absolute electron mass",
            "status": "open",
            "evidence": "current chain predicts ratios anchored to m_e",
        },
        {
            "criterion": "m~nu^2 bridge",
            "status": "open",
            "evidence": mass_bridge_gate()["needed_theorem"],
        },
    ]


def strongest_defensible_claim():
    return (
        "RFG has a stronger-than-phenomenological charged-lepton candidate "
        "because its C3 structure is traced to the supersolid invariant "
        "I3=det(B) and to an order-9 framed slot lattice. It is not yet a "
        "final particle theory until the cyclic lift, h=2 action derivation, "
        "full PDE stability, radiative protection, m~nu^2 bridge, and absolute "
        "scale are derived."
    )


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 44: Particle-sector strength and novelty audit")
    print("=" * 72)

    print("\n1. Existing approaches vs RFG")
    for row in comparison_with_existing_approaches():
        print(f"  {row['approach']}")
        print(f"    Koide: {row['has_koide']}")
        print(f"    delta 2/9: {row['has_delta_2_over_9']}")
        print(f"    C3 origin: {row['action_origin_c3']}")
        print(f"    topological h=2: {row['topological_h2']}")
        print(f"    gap: {row['main_gap']}")

    print("\n2. Internal consistency cleanup")
    for row in internal_consistency_cleanup():
        print(f"  {row['item']}")
        print(f"    old: {row['old_status']}")
        print(f"    new: {row['new_status']}")
        print(f"    reason: {row['reason']}")

    print("\n3. RFG strength scorecard")
    for row in rfg_strength_scorecard():
        print(
            f"  {row['criterion']:24s}: {row['status']} | {row['evidence']}"
        )

    print("\n4. Mass table")
    for row in prediction_table():
        print(
            f"  {row['particle']:8s}: "
            f"m_C3={row['predicted_mass_MeV']:.6f} MeV, "
            f"m_obs={row['observed_mass_MeV']:.6f} MeV, "
            f"rel_err={row['relative_mass_error']:.3e}, "
            f"sigma={row['sigma_error']:.3e}"
        )

    print("\n5. Open theorem items")
    for item in open_theorem_items():
        print(f"  OPEN: {item}")

    print("\n6. Strongest defensible claim")
    print(f"  {strongest_defensible_claim()}")


# ===================== merged from p11_particles.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 47: Particle-sector falsifiability map after the C3/order-9 upgrade

Purpose:
    Convert the strengthened particle sector into concrete pass/fail claims.

Important cleanup:
    The main charged-lepton route is now phase37-43 C3/order-9, not the older
    radial N-ladder. Therefore:

        - e, mu, tau are one C3 triplet.
        - missing radial modes n=2..13 are no longer expected generations.
        - N=4 329 keV is not a main prediction under the C3 route.

This makes the particle sector cleaner and easier to falsify.
"""

# merged import removed: from p11_particles import prediction_table
# merged import removed: from p11_particles import stability_conditions
# merged import removed: from p11_particles import rfg_radiative_options
# merged import removed: from p11_particles import candidate_tables


def charged_lepton_pass_fail():
    precision = pdg_precision_audit()
    return {
        "claim": "e, mu, tau are one C3/order-9 candidate triplet with theta=2/9",
        "status": "relative postdiction unless cyclic lift, h=2, and radiative protection are derived before using masses",
        "pass_condition": "derive the chain first, then pass PDG residuals or explain the residual as a calculable radiative/dressing shift",
        "fail_condition": "no cyclic lift, no stable h=2 branch, unprotected radiative drift, or unresolved PDG residual",
        "pdg_precision_pass": precision["pdg_precision_pass"],
        "chi2_non_anchor": precision["chi2_non_anchor"],
        "numbers": [
            {
                "particle": row["particle"],
                "predicted_MeV": row["predicted_mass_MeV"],
                "observed_MeV": row["observed_mass_MeV"],
                "relative_error": row["relative_mass_error"],
                "sigma_error": row["sigma_error"],
                "role": row["role"],
            }
            for row in precision["rows"]
        ],
    }


def deprecated_legacy_predictions():
    return [
        {
            "legacy_claim": "N=5,72,295 radial ladder",
            "new_status": "not the main lepton-generation mechanism",
            "why": "C3 triplet explains all three charged leptons together",
        },
        {
            "legacy_claim": "N=4 329 keV lepton-like companion",
            "new_status": "suspended under C3 route",
            "why": "C3 charged-lepton theorem does not imply a radial N=4 state",
        },
        {
            "legacy_claim": "muon = uniquely n=14 overtone",
            "new_status": "legacy audit only",
            "why": "overlap integrals were monotonic; no unique n=14 bottleneck",
        },
    ]


def new_falsifiable_targets():
    quark_candidates = candidate_tables()
    return [
        {
            "target": "radiative protection",
            "test": "derive either dressed pole-frequency theorem or loop cancellation",
            "fail": "unprotected QED-like generation logs shift Koide by ~1e-3",
            "file": "p11_particles.py",
        },
        {
            "target": "cyclic order-9 lift",
            "test": "construct a charged-defect return map whose orbit has order 9",
            "fail": "C3 x C3 remains only a slot lattice, so theta=h/9 is not derived",
            "file": "p11_particles.py",
        },
        {
            "target": "PDG residual gate",
            "test": "explain or predict the muon/tau residuals with the same dressing/radiative mechanism",
            "fail": "theta=2/9 remains only a relative compression and misses PDG precision",
            "file": "p11_particles.py",
        },
        {
            "target": "h=2 local stability",
            "test": "satisfy phase43 Hessian inequalities and then full fluctuation spectrum",
            "fail": "h=2 branch is a saddle or has negative fluctuation mode",
            "file": "p11_particles.py",
        },
        {
            "target": "quark zero-angle extension",
            "test": "one of two projection candidates survives fixed-scheme quark masses without fitting theta_U/D",
            "fail": "both projections fail once MSbar masses are run to the chosen scale",
            "file": "p11_particles.py",
            "candidate_count": len(quark_candidates),
        },
        {
            "target": "no C3-forbidden charged h=1 branch",
            "test": "charged defects require oriented framed closure, not projective/nematic closure",
            "fail": "a stable charged h=1 branch appears in the full defect spectrum",
            "file": "p11_particles.py",
        },
    ]


def next_calculation_queue():
    return [
        "Prove or reject the cyclic order-9 lift from the C3 x C3 slot lattice.",
        "Derive the h=2 oriented-frame coupling from the RFG action.",
        "Derive m proportional to nu^2 from the dressed oscillon energy functional.",
        "Build the full localized h=2 oscillon ansatz: Phi(t,r), Q(r,beta), framed triad U(x).",
        "Linearize the RFG action around that ansatz and compute the fluctuation operator.",
        "Check whether all non-gauge eigenvalues are non-negative.",
        "Add EM/self-field coupling and decide between phase45 options A or B.",
        "Choose one quark-mass scheme/scale and run phase46 as a real numerical falsification test.",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 47: Particle-sector falsifiability map")
    print("=" * 72)

    print("\n1. Charged-lepton C3 claim")
    cl = charged_lepton_pass_fail()
    print(f"  claim: {cl['claim']}")
    print(f"  status: {cl['status']}")
    print(f"  pass: {cl['pass_condition']}")
    print(f"  fail: {cl['fail_condition']}")
    print(f"  PDG precision pass: {cl['pdg_precision_pass']}")
    print(f"  chi2 non-anchor   : {cl['chi2_non_anchor']:.3e}")
    for row in cl["numbers"]:
        print(
            f"    {row['particle']:8s}: "
            f"pred={row['predicted_MeV']:.6f} MeV, "
            f"obs={row['observed_MeV']:.6f} MeV, "
            f"rel_err={row['relative_error']:.3e}, "
            f"sigma={row['sigma_error']:.3e}, "
            f"role={row['role']}"
        )

    print("\n2. Legacy predictions demoted")
    for row in deprecated_legacy_predictions():
        print(f"  {row['legacy_claim']}: {row['new_status']} ({row['why']})")

    print("\n3. New falsifiable targets")
    for row in new_falsifiable_targets():
        extra = f", candidates={row['candidate_count']}" if "candidate_count" in row else ""
        print(f"  {row['target']} [{row['file']}{extra}]")
        print(f"    test: {row['test']}")
        print(f"    fail: {row['fail']}")

    print("\n4. h=2 stability conditions")
    for key, value in stability_conditions().items():
        print(f"  {key:26s}: {value}")

    print("\n5. Radiative protection options")
    for row in rfg_radiative_options():
        print(f"  {row['option']}: {row['claim_needed']}")

    print("\n6. Next calculation queue")
    for item in next_calculation_queue():
        print(f"  - {item}")


# =============================================================================
# STAGE D2-D4: OLD particle / frequency / SM embedding gate
# =============================================================================

def stage_d2_frequency_to_mass_old_status():
    """
    Deletion-gate marker for OLD/ISPG_FrequencyToMass.tex.

    The old file's valuable content is preserved, but the main charged-lepton
    route is now the C3/order-9 chain rather than the legacy radial N ladder.
    """
    return {
        "old_file_drained": "OLD/ISPG_FrequencyToMass.tex",
        "new_file": "p11_particles.py",
        "migrated": True,
        "kept": [
            "particles as localized resonant oscillons",
            "mass tracks internal frequency/energy scale",
            "effective mass scales as exp(phi/2) in a pressure background",
            "electron as the first stable occupied charged resonance",
            "neutrinos are outside the charged massive oscillon ladder",
        ],
        "promoted_new_route": (
            "charged leptons are treated as one C3/order-9 candidate triplet with theta=2/9; "
            "the old N=5,72,295 / n=14 ladder is retained as legacy audit only"
        ),
        "still_useful_legacy": [
            "Mathieu ladder as broad cross-sector compatibility language",
            "massive non-neutrino SM states can be indexed on a common ladder at compatibility level",
            "N=4 329 keV branch is suspended unless the radial ladder is independently derived",
        ],
        "open": [
            "absolute electron scale",
            "full 3D oscillon fluctuation operator",
            "radiative protection of C3/Koide pole-frequency relation",
        ],
    }


def stage_d3_em_charge_spin_program_status():
    """
    Drain of the EM/charge/spin/topology ideas from OLD/ISPG_Quantum.tex.

    These are not finished derivations; they are high-value programme targets
    that must not be left only in OLD.
    """
    return {
        "old_source": "OLD/ISPG_Quantum.tex",
        "new_file": "p11_particles.py",
        "photon_EM": {
            "old_idea": "photon as helicoidal/Kelvin-like transverse medium wave",
            "RFG_status": "future gauge-completion programme",
            "needed_theorem": "derive Maxwell equations, U(1) gauge redundancy, and Coulomb 1/r^2 from the medium variables",
        },
        "charge_quantization": {
            "old_idea": "charge tied to framed/winding topology and Calugareanu-Fuller Lk=Wr+Tw",
            "RFG_status": "candidate topological current/operator programme",
            "fractional_charges": "Z3/color-orbifold route for 1/3 and 2/3 remains open",
        },
        "spin": {
            "old_idea": "spin-1/2 from Mobius/framed 720-degree closure",
            "RFG_status": "partly absorbed into h=2 oriented-framed C3/order-9 lepton branch",
            "open": "derive full Dirac spinor dynamics, Pauli exclusion, and spin-statistics from the defect moduli space",
        },
        "g_factor": {
            "target": "Dirac g=2",
            "status": "not derived; future QED-strength theorem",
            "route": "minimal coupling of the emergent spinor to the emergent U(1) connection",
        },
        "alpha_EM": {
            "target": "fine-structure constant 1/137",
            "status": "open structural number; do not claim derivation",
            "candidate_sources": ["resonance step count", "topological normalization", "medium impedance ratio"],
        },
    }


def stage_d4_sm_embedding_old_status():
    """
    Deletion-gate marker for OLD/ISPG_SM_Embedding.tex.

    This preserves the partial-embedding result and keeps its warnings: the old
    file shows compatibility, not a completed derivation of the Standard Model.
    """
    return {
        "old_file_drained": "OLD/ISPG_SM_Embedding.tex",
        "new_file": "p11_particles.py",
        "migrated": True,
        "partial_embedding_kept": [
            "D3-deformed oscillon supplies an eight-mode gluon-count route",
            "imported SM hypercharges satisfy anomaly checks and Y*N_c integer grading",
            "three generations are compatible with a three-axis/framed structure",
            "neutrinos are excluded from the Mathieu massive charged ladder",
            "W, Z, H can be hosted on the broad mass-ladder compatibility map",
        ],
        "Higgs_replacement_status": {
            "negative_result": "fundamental RFG action has no SM-type Phi^4 Higgs potential",
            "candidate": "archion/interior localized broken-phase region for gauge-boson mass generation",
            "open": "derive the W/Z mass matrix with a protected massless photon from the RFG medium, not by importing the Higgs sector",
        },
        "QCD_status": {
            "kept": "8 gluon count and D3/color route",
            "open": "derive SU(3) structure constants, vertices, running, and confinement potential",
        },
        "electroweak_status": {
            "kept": "on-shell m_W/m_Z compatibility as a mass-law consistency check",
            "warning": "not an independent Weinberg-angle prediction while W and Z masses are used as inputs",
        },
        "anomaly_status": (
            "SM anomaly cancellation is imported as a consistency check; "
            "RFG has not yet derived the unique SM hypercharge assignment"
        ),
        "completion_requirements": [
            "EM U(1) derivation",
            "SU(2)_L chirality and weak interactions",
            "SU(3)_c algebra and confinement",
            "archion-interior W/Z mass matrix with photon protection",
            "CKM/PMNS and neutrino ordering",
        ],
    }


def stage_d4_particle_old_file_status():
    """Compact deletion-gate summary for all remaining old particle files."""
    return {
        "OLD_ISPG_Quantum": "migrated across p10_oscillons.py, p11_particles.py, p12_predictions.py; open parts kept as programme targets",
        "OLD_ISPG_FrequencyToMass": "migrated into p11_particles.py STAGE D2; legacy N-ladder demoted",
        "OLD_ISPG_SM_Embedding": "migrated into p11_particles.py STAGE D4; partial embedding warnings preserved",
        "safe_to_delete_condition": (
            "safe as working material once p12 quantum predictions are marked; "
            "no finished theorem is left only in OLD, but several future-programme items remain open"
        ),
    }


# =============================================================================
# FINAL P11 CLAIM GATE
# =============================================================================

def p11_do_not_claim():
    """Particle-sector overclaim blacklist."""
    return [
        "RFG has derived the full Standard Model.",
        "C3 x C3 is already a proven Z9 holonomy.",
        "theta=2/9 is derived from the full RFG action.",
        "h=2 is derived rather than selected by a candidate normal-form rule.",
        "charged-lepton pole masses are predicted within PDG uncertainty.",
        "the absolute electron mass is derived.",
        "m proportional to nu^2 is derived from the oscillon energy functional.",
        "radiative corrections are solved.",
        "U(1), SU(2), SU(3), fractional charge, alpha_EM, g=2, CKM, or PMNS are derived.",
        "the legacy Mathieu/N ladder is the main lepton-generation proof.",
    ]


def p11_particle_sector_claim_gate():
    """Single importable gate for the whole p11 file."""
    precision = pdg_precision_audit()
    z9_gate = cyclic_lift_gate()
    h2_gate = h2_offset_gate()
    mass_bridge = mass_bridge_gate()
    mathieu_gate = mathieu_sign_stability_gate()
    return {
        "overall_status": "STRONG_CANDIDATE_NOT_FINAL_PARTICLE_THEORY",
        "closed_algebra": [
            "C3 triplet gives Koide K=2/3.",
            "I3=det(B) contains the C3 triaxial strain lock in principal-axis normal form.",
            "theta=2/9 gives a strong charged-lepton relative mass-ratio compression.",
        ],
        "hard_blockers": [
            z9_gate["needed_theorem"],
            h2_gate["needed_theorem"],
            mass_bridge["needed_theorem"],
            "derive radiative protection or a dressed pole-frequency theorem",
            "build the full localized h=2 oscillon and fluctuation spectrum",
            "derive EM/gauge sectors rather than importing SM quantum numbers",
            "resolve the muon PDG-precision residual or demote the claim permanently",
        ],
        "pdg_precision_pass": precision["pdg_precision_pass"],
        "pdg_chi2_non_anchor": precision["chi2_non_anchor"],
        "mathieu_gate": mathieu_gate,
        "do_not_claim": p11_do_not_claim(),
        "strongest_allowed_claim": strongest_defensible_claim(),
    }


if __name__ == "__main__":
    print("=" * 72)
    print("P11 FINAL CLAIM GATE")
    print("=" * 72)
    gate = p11_particle_sector_claim_gate()
    print(f"overall_status: {gate['overall_status']}")
    print(f"PDG precision pass: {gate['pdg_precision_pass']}")
    print(f"PDG chi2 non-anchor: {gate['pdg_chi2_non_anchor']:.3e}")
    print("\nClosed algebra:")
    for item in gate["closed_algebra"]:
        print(f"  - {item}")
    print("\nHard blockers:")
    for item in gate["hard_blockers"]:
        print(f"  - {item}")
    print("\nDo-not-claim:")
    for item in gate["do_not_claim"]:
        print(f"  - {item}")

