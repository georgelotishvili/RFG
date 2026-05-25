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
    - phase17_effective_mass.py (factor 2 emergence)
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
    print("RFG_Theory.md-სთვის: § 13-ის გვერდი — § 14 'LIGO და ფაქტორი 2'")
    print("=" * 72)
