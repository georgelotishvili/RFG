# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 9: gravitational waves - old ISPG predictions in RFG language.

Recovered old-theory predictions:
1. Tensor propagation is exactly luminal: alpha_T=0 -> c_g=c.
2. Two GR tensor modes remain.
3. A scalar/breathing channel is allowed by the scalar medium, but its
   source amplitude is parametrically suppressed.
4. Leading dipole radiation cancels when the compact-body sensitivity is
   universal, s=1/2; residuals require a full PN source calculation.

This file keeps the distinction sharp:
    c_g=c is closed structurally;
    breathing amplitude is a controlled working estimate;
    scalar-dipole nulling is leading-order, not a final waveform theorem.
"""

import sympy as sp


def analyze_gw_full():
    """
    Solid-sector TT check already used in earlier RFG phases.

    The RFG solid invariants do not generate h_dot^2 or h_z^2 corrections
    for a TT perturbation on FLRW.  Therefore the solid sector does not
    shift c_T.  A mass term constraint remains separate.
    """
    t, z = sp.symbols('t z', real=True)
    h = sp.Function('h')(t, z)
    h_dot = sp.diff(h, t)
    h_z = sp.diff(h, z)
    eps = sp.Symbol('eps', real=True)
    a = sp.Symbol('a', real=True, positive=True)

    g_11 = -a**2 * (1 - eps * h)
    g_22 = -a**2 * (1 + eps * h)
    g_33 = -a**2

    det_g = a**6 * (1 - eps**2 * h**2)
    sqrt_g = sp.sqrt(det_g)

    g_inv_11 = sp.series(1 / g_11, eps, 0, 3).removeO()
    g_inv_22 = sp.series(1 / g_22, eps, 0, 3).removeO()
    g_inv_33 = 1 / g_33

    b11 = -g_inv_11
    b22 = -g_inv_22
    b33 = -g_inv_33

    i1_pert = sp.simplify(b11 + b22 + b33)
    i2_pert = sp.simplify(sp.Rational(1, 2) * (i1_pert**2 - (b11**2 + b22**2 + b33**2)))
    i3_pert = sp.simplify(b11 * b22 * b33)
    y_pert = 1

    from p01_core import get_polynomial_lagrangian

    y_s, i1_s, i2_s, i3_s = sp.symbols('Y I1 I2 I3', real=True)
    l_poly = get_polynomial_lagrangian(y_s, i1_s, i2_s, i3_s)

    l_eval = l_poly.subs({y_s: y_pert, i1_s: i1_pert, i2_s: i2_pert, i3_s: i3_pert})
    l_density = sp.series(sqrt_g * l_eval, eps, 0, 3).removeO()
    l_o2 = sp.simplify(l_density.coeff(eps, 2) / a**3)

    coeff_h_dot2 = l_o2.coeff(h_dot**2)
    coeff_h_z2 = l_o2.coeff(h_z**2)
    mass_term = sp.simplify(l_o2.coeff(h**2))

    return coeff_h_dot2, coeff_h_z2, mass_term


def analyze_horndeski_luminal_speed():
    """
    Old Appendix 10 result in RFG notation.

    Horndeski tensor-speed excess alpha_T receives contributions from
    G_{4,X} and G_5.  RFG's Einstein-Hilbert backbone has:
        G4 = const, G4_X = 0, G5 = 0.
    """
    c, c_g, alpha_T, G4_X, G5 = sp.symbols('c c_g alpha_T G4_X G5', real=True)
    alpha_t_value = sp.Integer(0)
    c_g_value = c * sp.sqrt(1 + alpha_t_value)

    return {
        "theorem": "tensor gravitational waves are exactly luminal",
        "Horndeski_conditions": "G4=const, G4_X=0, G5=0",
        "alpha_T_definition": sp.Eq(alpha_T, c_g**2 / c**2 - 1),
        "G4_X": sp.Eq(G4_X, 0),
        "G5": sp.Eq(G5, 0),
        "alpha_T": sp.Eq(alpha_T, alpha_t_value),
        "c_g": sp.Eq(c_g, c_g_value),
        "GW170817_status": "satisfied structurally, not by parameter tuning",
    }


def analyze_scalar_breathing_estimate():
    """
    Breathing-mode working estimate inherited from the old theory.

    If scalar charge per mass is universal at leading order, s=1/2, the
    monopole is stationary and the leading dipole cancels.  The first
    candidate radiative scalar channel is the trace quadrupole, estimated
    as A_b/A_t ~ v^2/c^2 ~ r_s/r.
    """
    r, r_s, v, c, s_A, s_B = sp.symbols(
        'r r_s v c s_A s_B',
        real=True,
        positive=True,
    )
    sensitivity_universal = sp.Eq(s_A, sp.Rational(1, 2))
    sensitivity_match = sp.Eq(s_A - s_B, 0)
    amplitude_ratio = sp.simplify(v**2 / c**2)
    virial_ratio = r_s / r

    return {
        "scalar_channel": "breathing polarization h_b is allowed",
        "universal_sensitivity": sensitivity_universal,
        "dipole_charge_difference": sensitivity_match,
        "dipole_status": "leading dipole cancels when s_A=s_B=1/2",
        "amplitude_ratio_working": sp.Eq(sp.Symbol('A_b/A_t'), amplitude_ratio),
        "virial_estimate": sp.Eq(sp.Symbol('A_b/A_t'), virial_ratio),
        "weak_field_example_r_10rs": sp.N(virial_ratio.subs(r, 10 * r_s), 8),
        "status": "parametric working estimate; full PN scalar quadrupole coefficient remains a waveform task",
    }


def gw_prediction_ledger():
    return [
        "Closed: alpha_T=0 -> c_g=c exactly.",
        "Closed: the solid sector adds no h_dot^2 or h_z^2 TT kinetic-gradient correction.",
        "Constraint: the TT mass term must be tuned/constrained to avoid massive graviton dispersion.",
        "Recovered old estimate: scalar breathing amplitude A_b/A_t ~ r_s/r.",
        "Leading scalar dipole is suppressed if compact-body sensitivity is universal, s=1/2.",
        "Open waveform task: exact scalar quadrupole coefficient and comparable-mass IMR templates.",
        "ISCO proxy from phase18: f_ISCO=0.931 f_ISCO_GR is a strong-field timing target, not a full waveform by itself.",
    ]


if __name__ == "__main__":
    coeff_h_dot2, coeff_h_z2, mass_term_flrw = analyze_gw_full()
    a = sp.Symbol('a', real=True, positive=True)
    mass_term_mink = sp.simplify(mass_term_flrw.subs(a, 1))

    print("--- solid-sector TT check: c_T არ იცვლება ---")
    print(f"L_solid-ის კინეტიკური წევრი (h_dot^2): {coeff_h_dot2}")
    print(f"L_solid-ის გრადიენტული წევრი (h_z^2): {coeff_h_z2}")
    print("დასკვნა: L_solid არ შეიცავს h_dot^2 ან h_z^2 წევრებს.")
    print("c_T^2 = c^2 (1 + delta), სადაც delta = 0 ზუსტად სრულდება.")
    print("\nთუმცა L_solid წარმოქმნის გრავიტონის ეფექტურ მასას (m_g^2 * h^2):")
    print(f"FLRW ფონზე მასის კოეფიციენტი: {mass_term_flrw}")
    print(f"Minkowski ფონზე (a=1): {mass_term_mink}")
    print("GW170817-ის და მასიური დისპერსიის ასარიდებლად მოითხოვება ეს mass-term constraint.")

    print("\n--- Horndeski/EFT luminal theorem ---")
    for key, value in analyze_horndeski_luminal_speed().items():
        print(f"{key:28s}: {value}")

    print("\n--- scalar breathing / dipole ledger ---")
    for key, value in analyze_scalar_breathing_estimate().items():
        print(f"{key:28s}: {value}")

    print("\n--- prediction ledger ---")
    for item in gw_prediction_ledger():
        print(f"  - {item}")

# ===================== CONSOLIDATED PHASE SECTIONS =====================


# ===================== merged from p04_gw.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 18: LIGO-ის დიფერენციალური პასუხი როგორც spatial metric-sector consistency bridge
================================================================================

ცენტრალური დაკავშირება:
    ორივე შემთხვევაში spatial metric sector მნიშვნელოვანია, მაგრამ LIGO 
    პირდაპირ არ ამტკიცებს static light-deflection factor 2-ს.

ფიზიკური ჯაჭვი:
    1. ბი-კონფორმობა: g_tt + g_ii ერთად ცვლის, ფაქტორი 2
    2. სტატიკურად: m_eff ∝ e^(φ/2), L_oper ∝ e^(φ/2), მაგრამ c_coord ∝ e^φ
    3. სინათლის გავლის დრო ფიქს. კოორდინატულ მანძილზე: T ∝ e^(-φ)
       (ექსპონენტი 1 — *ორმაგი* L-ის 1/2-ის მიმართ)
    4. TT GW = წმინდა *სივრცული* პერტურბაცია (h_00 = 0)
    5. LIGO ხედავს: ΔT/T₀ = h_+ (დიფერენციალური სიგრძე-სხვაობა ცალმხრივად)

სტატუსი და ფიზიკური დასკვნა:
    ფაქტი, რომ LIGO ხედავს GW, წარმოადგენს ინტუიციურ/გეომეტრიულ 
    consistency bridge-ს თეორიისთვის. LIGO-ს სიგნალი *თავსებადია* 
    ბი-კონფორმულ სტრუქტურასთან (აუცილებელი, მაგრამ არა საკმარისი პირობა, 
    რადგან იგივეს აკეთებს GR). 

    RFG თეორიის GR-სგან გასარჩევად LIGO/Virgo-ს მონაცემებში საჭიროა 
    სკალარული/გრძივი დამატებითი მოდების ძიება ან ულტრა-ზუსტი დისპერსიის 
    გაზომვა. GW170817-მა (c_T = c) უკვე დაადგინა მკაცრი კონსტრეინტი; 
    phase9 მხოლოდ აჩვენებს, რომ solid sector TT kinetic/gradient-ს არ ცვლის; 
    სრული შედარება (extra modes, damping, polarizations, dispersion) ჯერ ღიაა.

References:
    - Intuitive_Theory.md §4 (c_coord ∝ L_oper²)
    - Intuitive_Theory.md §9.1 (LIGO-ის სიგნალის ფიზიკური წყარო)
    - OLD/0. MAIN.tex §sec:solar (light deflection 2× factor)
    - p10_oscillons.py (factor 2 emergence)
"""

import sympy as sp
from sympy import exp, sqrt, simplify, series, symbols, Symbol, Rational, pi


# ==============================================================================
# Setup
# ==============================================================================

def setup():
    """ბი-კონფორმული ფონი + TT გრავიტაციული ტალღის პერტურბაცია"""
    phi = Symbol('phi', real=True)          # ფონური წნევითი პოტენციალი
    h_p = Symbol('h_+', real=True)           # TT GW ამპლიტუდა (plus polarization)
    L_0 = Symbol('L_0', positive=True)       # კოორდინატული მკლავის სიგრძე
    c = Symbol('c', positive=True)           # სინათლის სიჩქარე
    return phi, h_p, L_0, c


# ==============================================================================
# ნაბიჯი 1: სტატიკური ბი-კონფორმობის შეჯამება (phase17-დან)
# ==============================================================================

def step1_static_scaling_recap():
    """
    ბი-კონფორმული მეტრიკიდან გამოდის:
        m_eff/m₀ = e^(φ/2)       (ექსპონენტი 1/2)
        L_oper/L₀ = e^(φ/2)      (ექსპონენტი 1/2)
        T_period/T₀ = e^(-φ/2)   (ექსპონენტი -1/2)
        c_coord/c = e^φ          (ექსპონენტი 1 ← *ორმაგი*)
    """
    phi = Symbol('phi', real=True)
    return {
        'm_eff':    exp(phi/2),
        'L_oper':   exp(phi/2),
        'T_period': exp(-phi/2),
        'c_coord':  exp(phi),                       # ← ფაქტორი 2
    }


# ==============================================================================
# ნაბიჯი 2: სინათლის გავლის დრო ფიქს. კოორდინატულ მანძილზე
# ==============================================================================

def step2_coordinate_light_travel_time():
    """
    ფიქსირებული კოორდინატული მკლავის სიგრძე L_0.
    სინათლის ერთმხრივი გავლის დრო:
        T = L_0 / c_coord = (L_0/c) · e^(-φ)

    ექსპონენტი -1 → "ორმაგი" L_oper-ის (-1/2)-თან შედარებით.

    LIGO რეალურად ზომავს phase/round-trip differential response-ს; აქ კი მოცემულია idealized coordinate-time toy derivation.
    """
    phi, h_p, L_0, c = setup()
    c_coord = c * exp(phi)
    T_light = L_0 / c_coord
    return simplify(T_light)


# ==============================================================================
# ნაბიჯი 3: TT GW-ის გავლენა მკლავის გავლის დროზე
# ==============================================================================

def step3_tt_gw_perturbation():
    """
    TT GW (plus polarization) ცვლის სივრცულ მეტრიკას:
        g_xx → e^(-φ)(1 + h_+)
        g_yy → e^(-φ)(1 - h_+)
        g_tt → -e^φ           ← UNCHANGED (h_00 = 0 TT-კალიბრში)

    სინათლის სიჩქარე x-ში (ds² = 0):
        c_coord_x² = -g_tt/g_xx = e^(2φ)/(1 + h_+)
        c_coord_x  = c · e^φ · 1/√(1 + h_+)
                   ≈ c · e^φ · (1 - h_+/2)         (სუსტ h_+ ლიმიტში)

    გავლის დრო L_0-ზე:
        T_x ≈ (L_0/c) · e^(-φ) · (1 + h_+/2)

    და y-ში:
        T_y ≈ (L_0/c) · e^(-φ) · (1 - h_+/2)
    """
    phi, h_p, L_0, c = setup()

    # c_coord_x ბი-კონფორმულ ფონზე + TT GW
    c_coord_x_sq = exp(2*phi) / (1 + h_p)
    c_coord_x = c * sp.sqrt(c_coord_x_sq)
    c_coord_x_lead = series(c_coord_x, h_p, 0, 2).removeO()

    c_coord_y_sq = exp(2*phi) / (1 - h_p)
    c_coord_y = c * sp.sqrt(c_coord_y_sq)
    c_coord_y_lead = series(c_coord_y, h_p, 0, 2).removeO()

    # Travel times
    T_x = L_0 / c_coord_x
    T_x_lead = simplify(series(T_x, h_p, 0, 2).removeO())

    T_y = L_0 / c_coord_y
    T_y_lead = simplify(series(T_y, h_p, 0, 2).removeO())

    return c_coord_x_lead, c_coord_y_lead, T_x_lead, T_y_lead


# ==============================================================================
# ნაბიჯი 4: LIGO-ის დიფერენციალური სიგნალი
# ==============================================================================

def step4_ligo_differential_signal():
    """
    LIGO ზომავს მკლავების სხვაობას (იდეალიზებული ცალმხრივი გავლისას):
        ΔT_oneway = T_x - T_y ≈ (L_0/c) · e^(-φ) · h_+
        ΔT/T_0 = h_+

    ეს არის LIGO-ის ფუნდამენტური ფიზიკური სიგნალი — დიფერენციალური სიგრძე-სხვაობა,
    რომელიც პროპორციულია h_+ GW ამპლიტუდის.
    """
    phi, h_p, L_0, c = setup()
    c_x, c_y, T_x, T_y = step3_tt_gw_perturbation()

    Delta_T = simplify(T_x - T_y)
    Delta_T_lead = simplify(series(Delta_T, h_p, 0, 2).removeO())

    # ფონური T_0 (no GW):
    T_0 = (L_0/c) * exp(-phi)

    # Relative strain
    relative_strain = simplify(Delta_T_lead / T_0)

    return Delta_T_lead, T_0, relative_strain


# ==============================================================================
# ნაბიჯი 5: ჰიპოთეტური "1911-სტილის" უნივერსი
# ==============================================================================

def step5_hypothetical_1911_only():
    """
    "1911-სტილის" უნივერსი არის მხოლოდ კონცეპტუალური ილუსტრაცია (straw-man).
    ნებისმიერი ცოცხალი მეტრიკული თეორია (მათ შორის GR) ეყრდნობა სივრცული 
    სექტორის აქტიურობას.

    ეს ილუსტრაცია უბრალოდ აჩვენებს, რომ გრავიტაციული ტალღის დასაფიქსირებლად
    აუცილებელია ტენზორული/სივრცული მეტრიკის ვარიაცია. ბი-კონფორმული სკალირება 
    ამ მოთხოვნას ბუნებრივად აკმაყოფილებს.
    """
    h_p = Symbol('h_+', real=True)
    L_0, c = symbols('L_0 c', positive=True)

    # "1911-სტილში" g_xx = 1 უცვლელი TT GW-ის ქვეშ
    # c_coord_x = c (უცვლელი)
    # T_x = L_0/c (უცვლელი)
    # T_y = L_0/c (უცვლელი)
    Delta_T_1911 = sp.Integer(0)  # ცხადადაა ნული

    return Delta_T_1911


# ==============================================================================
# ნაბიჯი 6: Einstein 1911 → 1915 ფაქტორი 2-თან კავშირი
# ==============================================================================

def step6_einstein_1911_1915_connection():
    """
    Einstein-ის ისტორიული ფაქტი:
        1911: მხოლოდ დროითი g_tt → სინათლის გადახრა = 0.87''
        1915: დროითი + სივრცული g_ii → სინათლის გადახრა = 1.75''
        ფაქტორი 2-ის წყარო = სივრცული სექტორი (g_ii)

    LIGO-ში იგივე სივრცული სექტორი მუშაობს:
        TT GW = წმინდა სივრცული პერტურბაცია (h_xx, h_yy)
        თუ მხოლოდ 1911-ის დროითი ნაწილი არსებობდა — LIGO ვერ ნახავდა
        ფაქტი, რომ ხედავს → სივრცული სექტორი აქტიურია

    გადახრის რიცხვობრივი ფაქტი (phase17-დან):
        δ_temporal (1911) = r_s/b
        δ_spatial (1915 add) = r_s/b
        δ_total = 2 r_s/b   ← ფაქტორი 2
    """
    return {
        '1911_only': '0.87 arcsec (temporal sector)',
        '1915_full': '1.75 arcsec (temporal + spatial)',
        'factor_2_source': 'spatial sector (g_ii)',
        'LIGO_consequence': 'TT GW = spatial perturbation → LIGO sees because g_ii is active'
    }


# ==============================================================================
# ნაბიჯი 7: GW150914 order-of-magnitude illustration
# ==============================================================================

def step7_gw150914_illustration():
    """
    GW150914 (LIGO პირველი დეტექცია, 2015):
        h_+ ≈ 1 × 10⁻²¹  (peak strain)
        L_arm = 4 km     (LIGO Hanford/Livingston)
        c = 3 × 10⁸ m/s

    პროგნოზი:
        ΔL = L_arm · h_+ ≈ 4 × 10⁻¹⁸ m   (= 4 attometer, << 0.001 of proton radius)
        
    რეალური ინტერფერომეტრი (Fabry-Pérot) ზომავს ორმხრივ (round-trip) დროს:
        ΔT_roundtrip = 2 · L_arm · h_+ / c ≈ 2.7 × 10⁻²⁶ s

    LIGO-მ ეს ცხადადაა გაზომა.
    """
    h_plus = 1e-21          # GW150914 peak strain
    L_arm = 4e3             # 4 km
    c_si = 2.998e8

    delta_L = L_arm * h_plus
    delta_T = 2 * L_arm * h_plus / c_si

    # შედარების მასშტაბები
    proton_radius = 0.84e-15  # m
    fraction_proton = delta_L / proton_radius

    return delta_L, delta_T, fraction_proton


# ==============================================================================
# დამატებითი: ფაქტორი 2 → LIGO ეპისტემოლოგიური დასკვნა
# ==============================================================================

def step8_epistemic_summary():
    """
    ცხადი ლოგიკური ჯაჭვი:

    (1) ბი-კონფორმობა: g_tt · g_ii = -c²
            ↓
    (2) m_eff, L_oper სკალირდება ექსპონენტით 1/2; c_coord — ექსპონენტით 1
            ↓
    (3) c_coord/c = (L_oper/L₀)²   [phase17 ნაბიჯი 6]
            ↓
    (4) სტატიკურად: სინათლის გადახრის ფაქტორი 2 (1.75''), Pound-Rebka
            ↓
    (5) დინამიკურად: TT GW = სივრცული პერტურბაცია, რომელიც აქტივობს
        c_coord-ის ცვლის გავლით
            ↓
    (6) LIGO ხედავს ΔT/T₀ = h_+ ≠ 0

    შესაბამისობით:
        LIGO სიგნალი თავსებადია ბი-კონფორმობასთან (აუცილებელი პირობა).
        თუმცა, იგივე სიგნალს პროგნოზირებს GR-იც.
        
    RFG-სა და GR-ის ემპირიული გამიჯვნა GW დეტექტორებში მოითხოვს:
        - GW170817-ით დადგენილი c_T = c (მოითხოვს phase9 კონსტრეინტს).
        - სკალარული ან გრძივი მოდების პოვნას.
    """
    pass


# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 18: LIGO-ის დიფერენციალური პასუხი როგორც consistency bridge")
    print("=" * 72)

    # ნაბიჯი 1
    print("\n--- ნაბიჯი 1: ბი-კონფორმობის სტატიკური სკალირება ---")
    scaling = step1_static_scaling_recap()
    print("  სიდიდე       სკალირება        ექსპონენტი")
    print("  ────────     ──────────       ──────────")
    for name, expr in scaling.items():
        print(f"  {name:<10} {str(expr):<15}")
    print("  → c_coord-ის ექსპონენტი (1) *ორმაგია* L_oper-ის (1/2)-თან")
    print("     ეს არის ფაქტორი 2 ემერჯენტული.")

    # ნაბიჯი 2
    print("\n--- ნაბიჯი 2: სინათლის გავლის დრო ფიქს. კოორდინატულ მანძილზე ---")
    T = step2_coordinate_light_travel_time()
    print(f"  T_light = L_0/c_coord = {T}")
    print(f"  ექსპონენტი e^(-φ) — *ორმაგი* L_oper-ის e^(φ/2)-თან")
    print(f"  შენიშვნა: აქ მოცემულია idealized coordinate-time toy derivation.")

    # ნაბიჯი 3
    print("\n--- ნაბიჯი 3: TT GW-ის გავლენა ---")
    c_x, c_y, T_x, T_y = step3_tt_gw_perturbation()
    print(f"  c_coord_x ≈ {c_x}")
    print(f"  c_coord_y ≈ {c_y}")
    print(f"  T_x ≈ {T_x}")
    print(f"  T_y ≈ {T_y}")

    # ნაბიჯი 4
    print("\n--- ნაბიჯი 4: LIGO-ის დიფერენციალური სიგნალი ---")
    delta_T, T_0, rel = step4_ligo_differential_signal()
    print(f"  ΔT_oneway = T_x - T_y = {delta_T}")
    print(f"  T_0 (no GW)    = {T_0}")
    print(f"  ΔT/T_0         = {rel}")
    print(f"  → დიფერენციალური სიგრძე-სხვაობა = h_+   ✓")
    print(f"  ეს არის LIGO-ის სიგნალის ცხადი მათემატიკური წყარო.")

    # ნაბიჯი 5
    print("\n--- ნაბიჯი 5: ჰიპოთეტური 1911-სტილის უნივერსი ---")
    delta_T_1911 = step5_hypothetical_1911_only()
    print(f"  მხოლოდ დროითი ეფექტი (g_ii უცვლელი):")
    print(f"    ΔT_oneway = {delta_T_1911}   ← LIGO ვერ ნახავდა GW-ს")
    print(f"  შენიშვნა: ეს არის straw-man ილუსტრაცია. რეალური თეორიები სივრცულ სექტორს შეიცავენ.")

    # ნაბიჯი 6
    print("\n--- ნაბიჯი 6: Einstein 1911 → 1915 კავშირი ---")
    einstein = step6_einstein_1911_1915_connection()
    for k, v in einstein.items():
        print(f"  {k:<22}: {v}")

    # ნაბიჯი 7
    print("\n--- ნაბიჯი 7: GW150914 order-of-magnitude illustration ---")
    dL, dT, frac = step7_gw150914_illustration()
    print(f"  GW150914: h_+ ≈ 10⁻²¹, L_arm = 4 km")
    print(f"  ΔL (პროგნოზი) = {dL:.2e} m")
    print(f"  ΔT_roundtrip (პროგნოზი) = {dT:.2e} s")
    print(f"  შედარება: {frac:.2e} × პროტონის რადიუსი (10⁻¹⁵ m)  [One-way vs Round-trip გამიჯნულია]")
    print(f"  → LIGO-მ რეალურად დააფიქსირა ეს მცირე დიფერენციალური ფაზური სხვაობა ✓")

    # შემაჯამებელი
    print("\n" + "=" * 72)
    print("ცენტრალური დასკვნა")
    print("=" * 72)
    print("""
    ფაქტი:  LIGO ხედავს GW-ს (GW150914, GW170817...)
            ↓
    შედეგი: სივრცული მეტრიკის სექტორი (g_ii) რეაგირებს გრავიტაციაზე
            ↓
    შედეგი: ბი-კონფორმული თეორიის (c_coord ∝ L_oper²) სტრუქტურა თავსებადია
            ამ დაკვირვებასთან (აუცილებელი პირობა დაცულია).

    LIGO-ის სიგნალი წარმოადგენს *consistency bridge*-ს RFG-სთვის და არა მის
    ექსკლუზიურ მტკიცებულებას (ვინაიდან GR-იც იგივე სიგნალს იძლევა).
    მთავარი ფილტრი — c_T = c (GW170817-დან) მოითხოვს phase9 კონსტრეინტის გათვალისწინებას; სრული ანალიზი ჯერ ღიაა.
    """)
    print("=" * 72)
    print("სტატიისთვის გამოსატანი ბლოკი: LIGO და ფაქტორი 2")
    print("=" * 72)


# ===================== merged from p04_gw.py =====================

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


# ===================== merged from p04_gw.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 27: Double Pulsar PSR J0737-3039A/B — PPK ცდები
================================================================================

რეფერენცია: p04_gw.py double-pulsar/PPK სამუშაო ბლოკი

დაკვირვება — PSR J0737-3039A/B (Kramer et al. 2021, PRX 11:041050):
- ორმაგი პულსარი (ერთადერთი ცნობილი)
- 6 post-Keplerian (PPK) პარამეტრი 0.05% precision
- GR ცდის ერთ-ერთი მკაცრი ფარგლი

PPK პარამეტრები:
1. ω̇ (periastron advance) — RFG: PPN γ, β-დან
2. γ (Einstein delay) — gravitational redshift + time dilation
3. P_b_dot (orbital period decay) — quadrupole + dipole radiation
4. r (Shapiro range) — phase13 1PN-დან
5. s (Shapiro shape) — sin(inclination)
6. Ω̇ (geodetic precession) — phase30 Lense-Thirring related

RFG-ის პროგნოზი:
- 1PN: GR-ის იდენტური (PPN γ=β=1)
- 2.5PN orbital decay: ცდის სქელეტი (phase26)
- Scalar dipole: phase9 ღია (s_A = 0 postulated)

ცდის ფარგლი GR vs RFG:
- 1PN ემთხვევა
- 2.5PN ცდის სქელეტი
- Scalar dipole — RFG-ის ღია ცდა
"""

import math


PSR_J0737_DATA = {
    "reference": "Kramer et al. 2021, PRX 11:041050",
    "P_b_orbital_period_days": 0.10225156248,
    "eccentricity": 0.0877775,
    "M_A_total_solar": 1.338185,
    "M_B_total_solar": 1.248868,
    "periastron_advance_deg_yr": 16.899323,  # ω̇
    "Einstein_delay_ms": 0.384045,  # γ
    "P_b_dot_obs": -1.247920e-12,  # observed orbital decay
    "P_b_dot_GR_pred": -1.247843e-12,  # GR quadrupole prediction
    "Shapiro_r_M_sun_s": 6.162e-6,  # r in solar mass · seconds
    "Shapiro_s": 0.999936,  # sin(i)
    "geodetic_omega_deg_yr": 4.78,  # B's geodetic precession
}


def gr_predictions():
    """GR-ის ცხადი 1PN პრედიქცია PSR J0737-სთვის."""
    G = 6.674e-11
    c = 2.998e8
    M_sun = 1.989e30
    M_A = PSR_J0737_DATA["M_A_total_solar"] * M_sun
    M_B = PSR_J0737_DATA["M_B_total_solar"] * M_sun
    M_tot = M_A + M_B
    P_b = PSR_J0737_DATA["P_b_orbital_period_days"] * 86400
    e = PSR_J0737_DATA["eccentricity"]

    # ω̇_GR = 3 · (GM/c²a)^(5/3) · (2π/P_b)^(5/3) / (1-e²)
    # სიმარტივისთვის — სიმბოლურად
    n = 2 * math.pi / P_b
    omega_dot_GR_rad = 3 * n ** (5 / 3) * (G * M_tot / c**3) ** (2 / 3) / (1 - e**2)
    omega_dot_GR_deg_yr = omega_dot_GR_rad * (180 / math.pi) * 86400 * 365.25

    return {
        "M_total_kg": M_tot,
        "P_b_seconds": P_b,
        "omega_dot_GR_deg_per_yr": omega_dot_GR_deg_yr,
        "omega_dot_observed": PSR_J0737_DATA["periastron_advance_deg_yr"],
    }


def rfg_predictions():
    """RFG-ის პრედიქცია — phase8-დან γ=β=1 ⇒ იდენტური GR."""
    return {
        "PPN_gamma": 1.0,  # phase8 (RFG bi-conformal)
        "PPN_beta": 1.0,  # phase8 2PN
        "omega_dot_RFG": "იდენტური GR-ის (γ=β=1)",
        "Einstein_delay_RFG": "იდენტური GR-ის (Pound-Rebka + gravitational time dilation)",
        "P_b_dot_quadrupole": "იდენტური GR-ის leading order (phase9 c_T=c)",
        "P_b_dot_dipole_RFG": "phase9-დან s_A = 0 (postulated) ⇒ no dipole",
        "P_b_dot_dipole_status": "OPEN — strong-field s_A derivation აკლია (Damour-Esposito-Farèse)",
    }


def open_tests():
    """ცდის ღია ნაბიჯები."""
    return [
        "Strong-field s_A derivation Damour-Esposito-Farèse სცენარით",
        "P_b_dot_dipole_RFG რიცხობრივი ცდა |α_0| < 2e-5 (PSR J1738+0333) ფარგლი",
        "Geodetic precession Ω̇ — Lense-Thirring (phase30) გადახედვა",
        "ცდის ცხადი interface PTA კოლაბორაციით (NANOGrav, EPTA)",
    ]


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 27: Double Pulsar PSR J0737-3039A/B")
    print("რეფერენცია: Kramer 2021, PRX 11:041050; p03_solar.py, p04_gw.py")
    print("=" * 72)

    print("\n1. დაკვირვება (PSR J0737-3039, Kramer et al. 2021)")
    for key, val in PSR_J0737_DATA.items():
        print(f"  {key:30s}: {val}")

    print("\n2. GR-ის 1PN პრედიქცია")
    gr = gr_predictions()
    for key, val in gr.items():
        if isinstance(val, float):
            print(f"  {key:30s}: {val:.6e}")
        else:
            print(f"  {key:30s}: {val}")

    print("\n3. RFG-ის პრედიქცია (phase8 + phase9)")
    rfg = rfg_predictions()
    for key, val in rfg.items():
        print(f"  {key:30s}: {val}")

    print("\n4. ღია ცდები")
    for i, task in enumerate(open_tests(), 1):
        print(f"  {i}. {task}")

    print("\n5. სტატუსი")
    print("  - 1PN: RFG=GR ფიქსირდება (phase8 γ=β=1)")
    print("  - 2.5PN orbital decay: GR=RFG leading order, dipole ცდის სქელეტი")
    print("  - Strong-field s_A — ღია (Damour-Esposito-Farèse)")
    print("  - PSR J0737-3039 ცდის სრული χ² fit — ღია")


# ===================== merged from p04_gw.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 28: PSR J1738+0333 — Scalar dipole radiation bound
================================================================================

რეფერენცია: p04_gw.py

დაკვირვება — PSR J1738+0333 (Freire et al. 2012, Antoniadis et al. 2013):
- pulsar-white dwarf binary
- highly asymmetric system (NS vs WD)
- scalar dipole-radiation-ის ცდის მკაცრი ფარგლი
- α_0 (scalar charge coupling) bound: |α_0| < 2 × 10⁻⁵ (95% CL)

scalar-tensor theory predictions:
- Brans-Dicke: dipole ∝ (α_NS - α_WD)²
- Damour-Esposito-Farèse: strong-field "spontaneous scalarization"
- RFG (phase9): s_A = 0 postulated ⇒ no dipole

RFG-ის ცდა:
1. dipole radiation: ΔP_b/P_b ∝ (s_A - s_B)² · (orbital velocity)²
2. RFG-ში s_A = ∂ln(m_A)/∂Φ
3. Minimal coupling ⇒ s_A = 0 (postulated)
4. Strong-field correction (Damour-Esposito-Farèse): s_A ~ Ω/(mc²)
5. PSR J1738 binding energy: Ω/(mc²) ~ 0.1 (NS) vs ~10^(-6) (WD)
6. ცდის ფარგლი: |α_0| ≤ 2 × 10⁻⁵
"""

import math


PSR_J1738_DATA = {
    "type": "NS + WD binary",
    "M_pulsar_solar": 1.46,
    "M_companion_solar": 0.181,
    "orbital_period_days": 0.354790739872,
    "eccentricity": 0.34e-6,  # very circular
    "P_b_dot_obs": -25.9e-15,  # observed
    "P_b_dot_GR": -27.7e-15,  # GR quadrupole prediction
    "alpha0_bound_95CL": 2e-5,  # scalar charge coupling upper bound
    "delta_alpha0_NS_WD_bound": 6e-3,  # asymmetry bound
}


def gr_quadrupole_prediction():
    """GR quadrupole-only orbital decay (no dipole)."""
    return {
        "P_b_dot_GR": PSR_J1738_DATA["P_b_dot_GR"],
        "ratio_obs_to_GR": PSR_J1738_DATA["P_b_dot_obs"] / PSR_J1738_DATA["P_b_dot_GR"],
        "comment": "PSR J1738 P_b_dot_obs/P_b_dot_GR = 0.94 ± 0.10 — GR consistent",
    }


def scalar_dipole_prediction():
    """
    Scalar-tensor dipole correction.

    ΔP_b_dot / P_b_dot_GR = -(4/3) · (G/c²) · (αA - αB)² · M_chirp² / (P_b² · G M_tot)^(2/3)

    Damour-Esposito-Farèse strong-field:
    αA = α_0 + β_0 · (Ω/mc²) + higher order
    """
    return {
        "Brans_Dicke_form": "ΔP_b ∝ (α_A - α_B)²",
        "RFG_postulated_phase9": "s_A = ∂ln(m_A)/∂Φ = 0 (minimal coupling)",
        "Damour_strong_field": "s_A = α_0 + β_0 · (Ω_A/(m_A c²))",
        "NS_binding_energy": "Ω_NS/(m_NS c²) ~ 0.1 (R~10 km, M~1.4 M_sun)",
        "WD_binding_energy": "Ω_WD/(m_WD c²) ~ 10^(-4) (R~10^4 km)",
        "asymmetry_bound": f"|s_NS - s_WD| < {PSR_J1738_DATA['delta_alpha0_NS_WD_bound']:.1e}",
        "alpha0_bound": f"|α_0| < {PSR_J1738_DATA['alpha0_bound_95CL']:.1e}",
    }


def rfg_strong_field_open():
    """RFG-ის strong-field s_A derivation — ცდის სქელეტი."""
    return [
        "Komar-integrand argument (phase9 Appendix 16) — s_A ≈ 1/2 leading order",
        "Bi-conformal weight e^(-φ) · ρ_0 — kinematic redshift + spatial volume",
        "Structural-response correction (open task per phase9)",
        "Strong-field NS: Ω/mc² ~ 0.1 — non-perturbative regime",
        "PSR J1738 χ² fit RFG s_A-დან — ცდა ღია",
        "Future: ngVLA + SKA pulsar timing → 100× precision",
    ]


def falsification_window():
    """RFG-ის ფალსიფიკაციის ფანჯარა PSR J1738-ში."""
    return {
        "current_bound": "|α_0| < 2 × 10⁻⁵ (PSR J1738)",
        "RFG_postulated": "s_A = 0 (consistent with bound)",
        "RFG_derived_value": "OPEN — strong-field calculation needed",
        "if_s_A_nonzero_NS_O(0.01)": "FALSIFIED (RFG αA-αB ≫ 2e-5)",
        "if_s_A_zero_strictly": "RFG ემთხვევა, არ ფალსიფიცირდება",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 28: PSR J1738+0333 — Scalar dipole bound")
    print("რეფერენცია: Freire 2012, Antoniadis 2013, phase9")
    print("=" * 72)

    print("\n1. დაკვირვება (PSR J1738+0333)")
    for key, val in PSR_J1738_DATA.items():
        print(f"  {key:30s}: {val}")

    print("\n2. GR quadrupole prediction")
    gr = gr_quadrupole_prediction()
    for key, val in gr.items():
        print(f"  {key:25s}: {val}")

    print("\n3. Scalar dipole prediction")
    dipole = scalar_dipole_prediction()
    for key, val in dipole.items():
        print(f"  {key:25s}: {val}")

    print("\n4. RFG strong-field s_A — ცდის ღია ნაბიჯები")
    for i, task in enumerate(rfg_strong_field_open(), 1):
        print(f"  {i}. {task}")

    print("\n5. ფალსიფიკაციის ფანჯარა")
    fals = falsification_window()
    for key, val in fals.items():
        print(f"  {key:30s}: {val}")

    print("\n6. სტატუსი")
    print("  - PSR J1738 |α_0| < 2e-5 დაფიქსირებულია")
    print("  - RFG phase9 s_A = 0 postulated (minimal coupling) — consistent")
    print("  - Strong-field s_A derivation — ღია (Damour-Esposito-Farèse-ის ანალოგი)")
    print("  - If derived s_A > 1e-5 in NS — RFG ფალსიფიცირდება")

