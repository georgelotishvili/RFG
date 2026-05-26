# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.
# Active coefficient scheme: Y-scheme; c_Y means c_Y^(Y), not X-scheme c_X.

"""
p02b: process-time ledger.

ეს ფაილი არის ცალკე ledger process-time იდეებისთვის. ის არ არის primary
FLRW/CMB metric branch და არ ამტკიცებს დამოუკიდებელ metric mode-ს.
"""

import sympy as sp


def get_process_time_relation():
    """
    ისტორიული process-time ხაზობრივი ესკიზის დაკონსოლიდებული ფორმა.

    აქტიური normalization ემთხვევა p10-ის ბი-კონფორმულ სკალირებას:
    d tau / dt = exp(phi/2). ეს ledger არ შედის primary FLRW/CMB metric
    branch-ში.
    """
    t = sp.Symbol("t", real=True)
    tau = sp.Function("tau")(t)
    phi = sp.Symbol("phi", real=True)

    relation = sp.Eq(sp.diff(tau, t), sp.exp(phi / 2))
    weak_field = sp.series(sp.exp(phi / 2), phi, 0, 2).removeO()

    return relation, weak_field, tau, t, phi


def process_time_scaling_ledger():
    """
    წნევის-ისტორიის process-time სკალირება.

    C(z) ნორმირებულია C(0)=1-ით და გამოიყენება მხოლოდ წინასწარ მონიშნულ
    intrinsic process-rate სიდიდეებზე.
    """
    z = sp.Symbol("z", real=True, nonnegative=True)
    phi_bg = sp.Function("phi_bg")
    C = sp.exp((phi_bg(z) - phi_bg(0)) / 2)

    return {
        "C_of_z": sp.Eq(sp.Function("C")(z), C),
        "process_clock": "d tau_proc / dt = C(z)",
        "rod_scale": "d ell(z) / d ell(0) = C(z)^(-1)",
        "mass_scale": "m_eff(z) / m_eff(0) = C(z)",
        "tail_power_perturbative": "P_tail(z) / P_tail(0) = C(z)^4",
        "normalization": "C(0)=1 by present-epoch calibration",
        "status": "bookkeeping only; not an independent metric mode",
    }


def flrw_vs_process_time_separation():
    """
    წესი, რომელიც process-time ledger-ს primary FLRW/CMB branch-ისგან ყოფს.
    """
    return {
        "metric_FLRW_branch": "phi_0(t)=0 matter-clock cosmic-time gauge-ში",
        "process_diagnostic": "phi_bg(z)=coarse-grained pressure-history diagnostic",
        "not_a_second_metric_mode": True,
        "CMB_rule": (
            "phi_bg(z) არ შეიტანო Einstein-Boltzmann/CMB სექტორში, როგორც "
            "დამატებითი ჰომოგენური metric perturbation."
        ),
        "where_it_acts": (
            "C(z) გამოიყენე მხოლოდ მონიშნულ intrinsic matter/resonance "
            "process-rate კანონებში, formation-history ledger-ში და tail-emission "
            "bookkeeping-ში."
        ),
    }


def process_time_integrals():
    """Coordinate lookback და pressure-weighted process-time ინტეგრალები."""
    z, zp = sp.symbols("z z_prime", real=True, nonnegative=True)
    H = sp.Function("H")
    C = sp.Function("C")

    lookback = sp.Integral(1 / ((1 + zp) * H(zp)), (zp, 0, z))
    process = sp.Integral(C(zp) / ((1 + zp) * H(zp)), (zp, 0, z))
    enhancement = sp.Symbol("A_proc")

    return {
        "coordinate_lookback_time": sp.Eq(sp.Function("t_lb")(z), lookback),
        "cumulative_process_time": sp.Eq(sp.Function("tau_proc")(z), process),
        "enhancement_factor": sp.Eq(enhancement, process / lookback),
        "past_reading": "თუ C(z)>1 წარსულში, intrinsic process history იზრდება",
        "future_reading": "თუ C(t)->0 საკმარისად სწრაფად, process-time ბიუჯეტი შეიძლება სასრული გახდეს",
    }


def rate_conversion_no_double_counting():
    """რომელ სიჩქარეში შედის C(z) და სად არ უნდა შევიდეს."""
    return [
        {
            "rate_tag": "intrinsic_process_time_rate",
            "symbolic_rule": "dX/dt = C(z) * dX/dtau_proc",
            "use_for": "კანონები, რომლებიც პირდაპირ tau_proc-ზეა განსაზღვრული",
            "double_counting_check": "გამოიყენე ზუსტად ერთი C(z) ფაქტორი",
        },
        {
            "rate_tag": "local_proper_time_rate",
            "symbolic_rule": "გადაყვანა მხოლოდ local clock calibration-ის მითითების შემდეგ",
            "use_for": "ადგილობრივი decay/oscillon/tail emission ერთ ეპოქაში",
            "double_counting_check": "არ დაუმატო process-time ფაქტორი, თუ კანონი tau_proc-ზე არ გადაიწერა",
        },
        {
            "rate_tag": "coordinate_cosmic_time_rate",
            "symbolic_rule": "დამატებითი C(z) არ შედის",
            "use_for": "FLRW background, Boltzmann evolution, merger crossing times",
            "double_counting_check": "dot უკვე დროის ცვლადს აფიქსირებს",
        },
        {
            "rate_tag": "observed_time_rate",
            "symbolic_rule": "რჩება სტანდარტული metric redshift/time dilation",
            "use_for": "დაკვირვებული სპექტრები, light curves, inferred rates",
            "double_counting_check": "observed redshift stretch არ გადაითარგმნოს process-time-ად",
        },
    ]


def process_time_application_map():
    """process-time ledger-ის გამოყენების რუკა სხვა კოსმოლოგიურ ბლოკებთან."""
    return {
        "MOND_vortex_maturation": (
            "coordinate-time H_eff რჩება ძირითად კალიბრაციად; process weighting "
            "მხოლოდ მონიშნულ intrinsic formation-rate კანონებში შედის"
        ),
        "JWST_high_z_galaxies": (
            "შეიძლება შეამციროს საჭირო intrinsic development time მხოლოდ მაშინ, "
            "თუ formation rate process-controlled-ად არის მონიშნული"
        ),
        "resonant_tail_dark_energy": (
            "tail-power C^4 არის bookkeeping input; P_tail-ის დროის ცვლადი "
            "ინტეგრაციამდე უნდა გამოცხადდეს"
        ),
        "cluster_mergers": (
            "merger crossing და redistribution times ადგილობრივი ეპოქის სიდიდეებია; "
            "არ გაამრავლო Bullet/cluster crossing times C(z)-ზე"
        ),
        "CMB": (
            "primary CMB იყენებს locked FLRW branch-ს; C(z) არ არის დამატებითი "
            "background metric mode"
        ),
    }


def stage_a2_process_time_status():
    """p02b-ში დარჩენილი process-time ledger-ის სრული სტატუსი."""
    return {
        "source_provenance": "OLD/18-derived ledger; not active authority",
        "scaling": process_time_scaling_ledger(),
        "flrw_separation": flrw_vs_process_time_separation(),
        "integrals": process_time_integrals(),
        "rate_table": rate_conversion_no_double_counting(),
        "applications": process_time_application_map(),
        "integration_status": "separate ledger in p02b; not part of p02 FLRW proof",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("p02b: process-time ledger")
    print("=" * 72)

    relation, weak_field, *_ = get_process_time_relation()
    print("process-time relation:", relation)
    print("weak-field limit:", weak_field)
    print("FLRW/CMB separation:", flrw_vs_process_time_separation()["not_a_second_metric_mode"])
    print("integration status:", stage_a2_process_time_status()["integration_status"])
