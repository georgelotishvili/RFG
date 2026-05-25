# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 26: Inspiral-Merger-Ringdown waveform — RFG vs GR smoke-test
================================================================================

სტატუსი:
Strategy 3 / M4-ის შესრულება. ეს ფაილი აღარ არის მხოლოდ PyCBC-ის ღია
ჩანაწერი: იგი ქმნის lightweight TaylorF2 inspiral waveform-ს, ამატებს RFG
ფაზურ correction-ებს და ითვლის overlap/mismatch-ს GR baseline-თან.

რისი მტკიცება შეიძლება ამ ფაილით:
    - 2.5PN GR TaylorF2 phase runnable არის.
    - RFG correction თუ მცირეა, mismatch მცირეა; თუ დიდია, LIGO-template
      consistency FAIL ხდება.
    - dipole/scalar/QNM ნაწილები პარამეტრიზებულია, მაგრამ მათი coupling ჯერ
      phase9/phase28/phase18 derivation-ზეა დამოკიდებული.

რისი მტკიცება ჯერ არ შეიძლება:
    - full IMRPhenom/SEOBNR waveform;
    - რეალური PyCBC catalog fit;
    - RFG-specific 2PN/3PN coefficient-ის derived მნიშვნელობა.
"""

import math

import numpy as np


MTSUN_SI = 4.92549095e-6  # G*M_sun/c^3 in seconds


def integrate_trapezoid(values, x_values):
    """NumPy compatibility wrapper."""
    if hasattr(np, "trapezoid"):
        return np.trapezoid(values, x_values)
    return np.trapz(values, x_values)


def gw_observations():
    return {
        "GW150914": {
            "type": "BBH",
            "m1_msun": 36.0,
            "m2_msun": 29.0,
            "f_low_Hz": 20.0,
            "f_high_Hz": 250.0,
        },
        "GW170817": {
            "type": "BNS",
            "m1_msun": 1.46,
            "m2_msun": 1.27,
            "f_low_Hz": 20.0,
            "f_high_Hz": 1200.0,
            "ct_bound": "|c_T/c - 1| < O(1e-15)",
        },
        "GW190521": {
            "type": "heavy BBH",
            "m1_msun": 85.0,
            "m2_msun": 66.0,
            "f_low_Hz": 11.0,
            "f_high_Hz": 120.0,
        },
    }


def binary_params(m1_msun, m2_msun):
    total = m1_msun + m2_msun
    eta = (m1_msun * m2_msun) / total**2
    chirp = total * eta ** (3.0 / 5.0)
    return {
        "M_total_msun": total,
        "eta": eta,
        "M_chirp_msun": chirp,
        "M_seconds": total * MTSUN_SI,
    }


def taylorf2_phase_25pn(freq_hz, m1_msun, m2_msun, tc=0.0, phic=0.0):
    """
    Frequency-domain TaylorF2 phase through 2.5PN.

    psi(f) = 2*pi*f*tc - phic - pi/4
             + 3/(128*eta)*v^-5 * PN(v)
    v = (pi*M*f)^(1/3)
    """
    params = binary_params(m1_msun, m2_msun)
    eta = params["eta"]
    M_sec = params["M_seconds"]
    v = np.power(np.pi * M_sec * freq_hz, 1.0 / 3.0)

    pn_0 = 1.0
    pn_2 = (3715.0 / 756.0 + 55.0 * eta / 9.0) * v**2
    pn_3 = -16.0 * np.pi * v**3
    pn_4 = (
        15293365.0 / 508032.0
        + 27145.0 * eta / 504.0
        + 3085.0 * eta**2 / 72.0
    ) * v**4
    pn_5 = np.pi * (38645.0 / 756.0 - 65.0 * eta / 9.0) * v**5

    phase = (
        2.0 * np.pi * freq_hz * tc
        - phic
        - np.pi / 4.0
        + (3.0 / (128.0 * eta)) * v ** (-5.0) * (pn_0 + pn_2 + pn_3 + pn_4 + pn_5)
    )
    return phase


def rfg_phase_correction(freq_hz, m1_msun, m2_msun, params):
    """
    Parametric RFG phase corrections.

    beta_dipole:
        -1PN dipole-like correction. Must be tiny unless phase28 strong-field
        scalar charge derivation proves otherwise.

    beta_2pn / beta_3pn:
        phenomenological higher-PN phase deviations.

    scalar_breathing_amp:
        amplitude channel, not phase; reported separately.
    """
    bin_params = binary_params(m1_msun, m2_msun)
    M_sec = bin_params["M_seconds"]
    v = np.power(np.pi * M_sec * freq_hz, 1.0 / 3.0)

    beta_dipole = params.get("beta_dipole", 0.0)
    beta_2pn = params.get("beta_2pn", 0.0)
    beta_3pn = params.get("beta_3pn", 0.0)

    return beta_dipole * v ** (-7.0) + beta_2pn * v ** (-1.0) + beta_3pn * v


def waveform_frequency_domain(freq_hz, m1_msun, m2_msun, rfg_params=None):
    """
    Restricted-amplitude TaylorF2 waveform h(f) = A f^(-7/6) exp(i psi).
    """
    if rfg_params is None:
        rfg_params = {}
    amp = np.power(freq_hz, -7.0 / 6.0)
    phase = taylorf2_phase_25pn(freq_hz, m1_msun, m2_msun)
    phase = phase + rfg_phase_correction(freq_hz, m1_msun, m2_msun, rfg_params)
    return amp * np.exp(1j * phase)


def toy_psd(freq_hz):
    """
    Smooth analytic PSD-like weight for aLIGO band.
    It is not a substitute for a real detector PSD; it makes overlap runnable.
    """
    x = freq_hz / 215.0
    return x ** (-4.14) - 5.0 * x ** (-2.0) + 111.0 * (1.0 - x**2 + 0.5 * x**4) / (1.0 + 0.5 * x**2)


def inner_product(h1, h2, freq_hz):
    psd = np.maximum(toy_psd(freq_hz), 1.0e-46)
    integrand = np.real(h1 * np.conjugate(h2)) / psd
    return 4.0 * integrate_trapezoid(integrand, freq_hz)


def normalized_overlap(h1, h2, freq_hz):
    norm_11 = inner_product(h1, h1, freq_hz)
    norm_22 = inner_product(h2, h2, freq_hz)
    norm_12 = inner_product(h1, h2, freq_hz)
    if norm_11 <= 0 or norm_22 <= 0:
        return float("nan")
    return norm_12 / math.sqrt(norm_11 * norm_22)


def overlap_smoke_test(event, rfg_params, n_freq=4096):
    freqs = np.linspace(event["f_low_Hz"], event["f_high_Hz"], n_freq)
    h_gr = waveform_frequency_domain(freqs, event["m1_msun"], event["m2_msun"])
    h_rfg = waveform_frequency_domain(freqs, event["m1_msun"], event["m2_msun"], rfg_params)
    overlap = normalized_overlap(h_gr, h_rfg, freqs)
    mismatch = 1.0 - overlap
    return {
        "overlap": overlap,
        "mismatch": mismatch,
        "status": "PASS" if mismatch < 0.03 else "FAIL",
    }


def qnm_ringdown_shift(final_mass_msun, epsilon_core=0.0):
    """
    Schwarzschild l=2 ringdown frequency smoke-test.
    GR: omega_220*M = 0.37367 - 0.08896 i.
    epsilon_core is a phenomenological RFG regular-core fractional shift.
    """
    M_sec = final_mass_msun * MTSUN_SI
    omega_real_gr = 0.37367 / M_sec
    omega_imag_gr = -0.08896 / M_sec
    f_gr_hz = omega_real_gr / (2.0 * np.pi)
    tau_gr_s = -1.0 / omega_imag_gr

    f_rfg_hz = f_gr_hz * (1.0 + epsilon_core)
    tau_rfg_s = tau_gr_s / max(1.0 + epsilon_core, 1.0e-12)

    return {
        "f_220_GR_Hz": f_gr_hz,
        "tau_GR_s": tau_gr_s,
        "epsilon_core": epsilon_core,
        "f_220_RFG_Hz": f_rfg_hz,
        "tau_RFG_s": tau_rfg_s,
        "ringdown_status": "PASS" if abs(epsilon_core) < 0.30 else "FAIL",
    }


def scalar_breathing_channel(rfg_params):
    amp = abs(rfg_params.get("scalar_breathing_amp", 0.0))
    return {
        "A_breathing_over_A_TT": amp,
        "current_status": "parameterized only; phase9/phase28 strong-field source derivation needed",
        "ligo_smoke_bound": "PASS" if amp < 0.10 else "FAIL",
    }


def pycbc_interface_open():
    return [
        "replace toy PSD with detector PSD from PyCBC",
        "maximize overlap over time/phase analytically or with matched_filter",
        "scan beta_dipole, beta_2pn, beta_3pn against LVK posterior samples",
        "connect beta_dipole to phase28 scalar-charge derivation",
        "connect epsilon_core to phase18 regular-BH metric and QNM calculation",
    ]


def benchmark_rfg_models():
    return {
        "GR_limit": {
            "beta_dipole": 0.0,
            "beta_2pn": 0.0,
            "beta_3pn": 0.0,
            "scalar_breathing_amp": 0.0,
            "epsilon_core": 0.0,
        },
        "small_RFG_deviation": {
            "beta_dipole": 1.0e-6,
            "beta_2pn": 2.0e-3,
            "beta_3pn": 1.0e-3,
            "scalar_breathing_amp": 0.02,
            "epsilon_core": 0.05,
        },
        "excluded_large_deviation": {
            "beta_dipole": 1.0e-3,
            "beta_2pn": 0.20,
            "beta_3pn": 0.10,
            "scalar_breathing_amp": 0.20,
            "epsilon_core": 0.50,
        },
    }


def status_assessment():
    return {
        "closed_now": "TaylorF2 2.5PN phase, parametric RFG phase corrections, overlap and QNM smoke-tests.",
        "still_open": "real PyCBC/LVK catalog fit and RFG-derived beta coefficients.",
        "falsification": "large beta_dipole/beta_PN/scalar/QNM shifts fail waveform overlap or ringdown bounds.",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 26: IMR waveform — RFG vs GR smoke-test")
    print("=" * 72)

    observations = gw_observations()
    print("\n1. Observational anchors")
    for name, event in observations.items():
        params = binary_params(event["m1_msun"], event["m2_msun"])
        print(
            f"  {name:8s}: type={event['type']:10s} "
            f"M={params['M_total_msun']:.2f} Msun, eta={params['eta']:.4f}, "
            f"band={event['f_low_Hz']:.0f}-{event['f_high_Hz']:.0f} Hz"
        )

    print("\n2. TaylorF2 + RFG overlap smoke-test on GW150914-like BBH")
    event = observations["GW150914"]
    for name, model in benchmark_rfg_models().items():
        ov = overlap_smoke_test(event, model)
        breath = scalar_breathing_channel(model)
        qnm = qnm_ringdown_shift(final_mass_msun=62.0, epsilon_core=model["epsilon_core"])
        print(f"\n  --- {name} ---")
        print(f"    overlap        : {ov['overlap']:.6f}")
        print(f"    mismatch       : {ov['mismatch']:.6e} -> {ov['status']}")
        print(f"    breathing amp  : {breath['A_breathing_over_A_TT']:.3f} -> {breath['ligo_smoke_bound']}")
        print(f"    ringdown f_GR  : {qnm['f_220_GR_Hz']:.2f} Hz")
        print(f"    ringdown f_RFG : {qnm['f_220_RFG_Hz']:.2f} Hz -> {qnm['ringdown_status']}")

    print("\n3. GW170817 speed/dipole guard")
    bns = observations["GW170817"]
    print(f"  c_T filter: {bns['ct_bound']}")
    print("  dipole filter: beta_dipole must remain near zero unless phase28 derives tiny scalar charge.")

    print("\n4. PyCBC/LVK open interface")
    for i, task in enumerate(pycbc_interface_open(), 1):
        print(f"  {i}. {task}")

    print("\n5. Status")
    for key, value in status_assessment().items():
        print(f"  {key:14s}: {value}")
