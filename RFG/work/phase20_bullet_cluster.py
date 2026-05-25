# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 20: Bullet Cluster - chi ველის მეხსიერების მექანიკა
================================================================================

სტატუსი:
ეს ფაილი არის კონცეპტუალური Toy Model (საილუსტრაციო მოდელი). ის აჩვენებს
მეხსიერების ეფექტის ლოგიკას, მაგრამ არ არის რეალური κ(x,y) lensing map-ის
სიმულაცია Clowe et al. მონაცემებზე.

ცენტრალური ფაქტი (1E 0657-558, 2006):
    ორი გალაქტიკის გროვა შეჯახდა ~150 მილიონი წლის წინ.
    - გაზი (X-ray, ~90% ბარიონული მასა): შენელდა ram pressure-ით
    - გალაქტიკები (~10%): გაიარეს თავისუფლად (collisionless)
    - გრავიტაციული ლენზირება: პიკი გალაქტიკებთან,
      *არა* გაზთან, სადაც ბარიონული მასის უმეტესობაა

ეს არის *მთავარი წინააღმდეგობა MOND-ისთვის*:
    MOND-ის დამატებითი გრავიტაცია ეფუძნება ლოკალურ ბარიონულ მასას.
    გაზი (უმეტესობა) -> MOND მოელის ლენზირების პიკს გაზთან.
    დაკვირვება: პიკი გალაქტიკებთან -> MOND-ის სტანდარტული ფორმა ვერ ხსნის.

RFG-ის მექანიზმი:
    chi ველი ვორტექსულ-მეხსიერებითი ფაზაა (PLAN.md ფაზა 4).
    მისი თავისებურება: *მეხსიერების შენარჩუნება დიდი დროით*
    (tau_chi >> t_collision).

    chi-ის რეაგირების განტოლება:
        Box chi + (1/tau_chi) * (chi - chi_eq) = J

    სადაც chi_eq წყარო არის კინემატიკური/ვორტიკალური დინება (omega dot v).
    ეს არის სამუშაო ჰიპოთეზა: გაზი შოკის დროს თერმალიზდება და კარგავს ვორტიციტს,
    ხოლო გალაქტიკები ინარჩუნებენ კინემატიკურ სტრუქტურას.
    tau_chi-ის სიდიდე (თავისუფალი პარამეტრი) განსაზღვრავს, რამდენ ხანში
    "მიჰყვება" chi წყაროს.

დინამიკა:
    tau_chi ~ Gyr >> t_collision ~ 100 Myr -> chi ინარჩუნებს თავის
    პოზიციას (გალაქტიკებთან), არ მიჰყვება გაზს.

შესაბამისად:
    Phi_lens(x) = Phi_baryon(x) + Phi_chi(x)
    პიკი:     gas-ში (~90%) + გალაქტიკები-შიგნით (~10%) + chi (გალაქტიკებთან)
    შედეგი:   გალაქტიკები + chi >> გაზი -> პიკი გალაქტიკებთან

References:
    - Clowe et al. 2006 (ApJL 648:L109) - Bullet Cluster discovery
    - Intuitive_Theory.md §6.2 - vortex memory mechanism
    - PLAN.md ფაზა 4 - chi field role
"""

import sympy as sp
from sympy import Heaviside, Symbol, exp, oo, simplify, symbols


# ==============================================================================
# Setup
# ==============================================================================

def setup_collision():
    """შეჯახების ბაზური პარამეტრები."""
    t, x = symbols("t x", real=True)
    v_coll = Symbol("v_coll", positive=True)  # შეჯახების სიჩქარე
    L_cluster = Symbol("L_cluster", positive=True)
    tau_chi = Symbol("tau_chi", positive=True)  # თავისუფალი/საერთო პარამეტრი
    return t, x, v_coll, L_cluster, tau_chi


# ==============================================================================
# ნაბიჯი 1: კომპონენტთა დინამიკა შეჯახებაში
# ==============================================================================

def step1_collision_dynamics():
    """
    შეჯახების სამი კომპონენტი:

    1. გაზი (hydrodynamic):
       v_gas(t) = v_0 * exp(-t/tau_drag)  <- ram pressure-ით შენელება
       tau_drag ~ 10 Myr (გროვის გადაკვეთის დროზე მცირე)
       -> გაზი ცენტრში გროვდება

    2. გალაქტიკები (collisionless):
       v_gal(t) = v_0  <- მუდმივი (ბალისტიკური მოძრაობა)
       -> გალაქტიკები გაივლიან თავიდან ბოლომდე

    3. chi ველი (memory):
       v_chi მიბმულია გალაქტიკების კინემატიკურ სტრუქტურაზე,
       რადგან tau_chi მნიშვნელოვნად დიდია.
       -> chi ჩამორჩება გაზს და მიჰყვება კინემატიკურ სტრუქტურას.
    """
    t, v_0, tau_drag = symbols("t v_0 tau_drag", positive=True)

    # გაზის სიჩქარე (exponential damping)
    v_gas = v_0 * exp(-t / tau_drag)
    x_gas = simplify(sp.integrate(v_gas, (t, 0, t)))

    # გალაქტიკების ბალისტიკური მოძრაობა
    v_gal = v_0
    x_gal = v_0 * t

    # გალაქტიკებსა და გაზს შორის დაშორება
    separation = simplify(x_gal - x_gas)

    return v_gas, x_gas, v_gal, x_gal, separation


# ==============================================================================
# ნაბიჯი 2: chi ველის რეაგირების განტოლება
# ==============================================================================

def step2_chi_field_equation():
    """
    chi ველის რეაგირების განტოლება (toy model):
        d chi / dt = -(1/tau_chi) * (chi - chi_eq(x,t))

    სადაც:
        chi_eq - წონასწორული მნიშვნელობა, რომელიც მოდის ვორტიკალური დინებიდან
                 (სამუშაო ჰიპოთეზა)
        tau_chi - რელაქსაციის/მეხსიერების დრო

    ინტეგრალური ფორმა Green-ის ფუნქციით:
        G(t,t') = exp(-(t - t')/tau_chi) * Heaviside(t - t') / tau_chi
        chi(t) = chi(0) * exp(-t/tau_chi) + integral_0^t G(t,t') * chi_eq(t') dt'

    აგენტთა საბჭოს დაზუსტება:
    vorticity source (omega dot v) არის ჰიპოთეზა, კოდში ჯერ არ არის გამოყვანილი.
    ასევე, ცვლადი წყაროს ინტეგრალური ფორმა აქ არ მოწმდება; კოდი იყენებს
    მხოლოდ constant-source toy solution-ს.
    """
    t, tau_chi = symbols("t tau_chi", positive=True)
    chi_0, omega_v, k_chi = symbols("chi_0 omega_v k_chi", real=True)
    t_p = Symbol("t_prime", positive=True)

    # Causal Green kernel
    green_kernel = exp(-(t - t_p) / tau_chi) * Heaviside(t - t_p) / tau_chi

    # ვორტიკალური წყარო კოდში მონიშნულია როგორც სიმბოლური მუდმივა.
    chi_eq_const = k_chi * omega_v

    # ინტეგრალური ამოხსნა მუდმივი წყაროსთვის (test case)
    chi_t = chi_0 * exp(-t / tau_chi) + chi_eq_const * (1 - exp(-t / tau_chi))
    chi_t = simplify(chi_t)

    # მეხსიერების რეჟიმი: t << tau_chi
    chi_short = simplify(chi_t.series(t, 0, 2).removeO())

    # რელაქსირებული რეჟიმი: t >> tau_chi
    chi_long = sp.limit(chi_t, t, oo)

    return chi_t, chi_short, chi_long, green_kernel


# ==============================================================================
# ნაბიჯი 3: Bullet Cluster-ის რიცხვობრივი მასშტაბები
# ==============================================================================

def step3_bullet_cluster_timescales():
    """
    Bullet Cluster 1E 0657-558 პარამეტრები:
        v_collision ≈ 4700 km/s        (relative velocity)
        L_separation ≈ 720 kpc         (current gas-galaxy separation)
        t_since_collision ≈ 150 Myr

    გროვის გადაკვეთის დრო:
        t_cross = L/v ≈ 720 kpc / 4700 km/s ≈ 150 Myr

    chi-ის რელაქსაციის დრო (უნივერსალური/თავისუფალი პარამეტრი):
        tau_chi ≳ 1 Gyr
        -> t_collision / tau_chi ≈ 0.15
        -> chi თითქმის გაყინულია შეჯახების დროით მასშტაბზე.
    """
    v_coll_km_s = 4700
    L_sep_kpc = 720
    t_since_collision_Myr = 150
    tau_chi_Gyr = 1.0  # საილუსტრაციო მნიშვნელობა; რეალურად თავისუფალი პარამეტრია

    # Cross time
    L_sep_km = L_sep_kpc * 3.086e16
    t_cross_s = L_sep_km / v_coll_km_s
    t_cross_Myr = t_cross_s / 3.156e13

    # Ratio
    ratio = t_since_collision_Myr / (tau_chi_Gyr * 1000)

    return v_coll_km_s, L_sep_kpc, t_cross_Myr, tau_chi_Gyr, ratio


# ==============================================================================
# ნაბიჯი 4: ლენზირების პოტენციალი - Toy Model შედარება
# ==============================================================================

def step4_toy_lensing_peak_comparison():
    """
    ლენზირების პოტენციალის Toy Model შედარება:
        Phi_lens(x) = Phi_baryon(x) + Phi_chi(x)

    Bullet Cluster-ის შემთხვევაში:
        Position           ბარიონული        chi ვორტექსი    ჯამი
        --------------     ----------       ------------    ----------
        გაზის ცენტრში       90%              ~0%             90%
        გალაქტიკის ცენტრში  10%              100% memory     110%

    შედეგი: ლენზირების პიკი ფიქსირდება გალაქტიკის ცენტრში (110% > 90%).
    შენიშვნა: ეს არის მხოლოდ პროცენტული საილუსტრაციო მოდელი. სრული მტკიცებისთვის
    საჭიროა რეალური κ(x,y) რუკის დათვლა (shear map, gas profile,
    galaxy distribution, projected mass) და რელატივისტური პოტენციალების შედარება.
    """
    Phi_baryon_gas, Phi_baryon_gal, Phi_chi = symbols(
        "Phi_baryon_gas Phi_baryon_gal Phi_chi", positive=True
    )

    Phi_at_gas_center = 0.9 * Phi_baryon_gas + 0 * Phi_chi
    Phi_at_gal_center = 0.1 * Phi_baryon_gal + 1.0 * Phi_chi

    return Phi_at_gas_center, Phi_at_gal_center


# ==============================================================================
# ნაბიჯი 5: MOND-ის კონტრასტი
# ==============================================================================

def step5_mond_contrast():
    """
    MOND (Milgrom) vs RFG-chi:
        MOND:  a_extra ∝ sqrt(g_N_baryonic) -> ბარიონების უმეტესობა გაზშია,
               პიკი უნდა იყოს გაზთან.
        RFG:   a_extra ∝ sqrt(g_chi_memory) -> პიკი რჩება გალაქტიკებთან
               chi ვორტექსული მეხსიერების გამო.
    """
    return {
        "MOND_failure": "a_extra ∝ sqrt(g_N_baryonic) -> peak at gas (ეწინააღმდეგება დაკვირვებას)",
        "RFG_success": "a_extra ∝ sqrt(g_chi_memory) -> peak at galaxies (რელატივისტური შედარება ღიაა)",
        "mechanism": "chi memory + vortical source -> decouples from thermalized gas",
    }


# ==============================================================================
# ნაბიჯი 6: ფალსიფიცირებადი პროგნოზები
# ==============================================================================

def step6_falsifiable_predictions():
    """
    RFG-ის ფალსიფიცირებადი პროგნოზები (ტექსტური ესკიზი):
    """
    return {
        "Bullet Cluster": "peak at galaxies",
        "Abell 520 (Train Wreck)": "peak between gas and galaxies (მეხსიერების რელაქსაციის ტესტი)",
        "MACS J0717.5 / El Gordo": "დამატებითი ტესტები tau_chi-ის შესაზღუდად",
        "tau_chi Parameter": "უნდა იყოს უნივერსალური/საერთო პარამეტრი ყველა გროვისთვის",
        "Falsification": "თუ tau_chi მოითხოვს ინდივიდუალურ მორგებას ყოველ გროვაზე, თეორია ფალსიფიცირებულია.",
    }


# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 20: Bullet Cluster - chi ველის მეხსიერების მექანიკა")
    print("=" * 72)

    # ნაბიჯი 1
    print("\n--- ნაბიჯი 1: შეჯახების სამი კომპონენტი ---")
    v_gas, x_gas, v_gal, x_gal, sep = step1_collision_dynamics()
    print(f"  გაზი: v(t) = {v_gas}")
    print(f"  გალაქტიკები: v = {v_gal}, x = {x_gal}")
    print("  -> შედეგი: გაზი ნელდება, გალაქტიკები წინ მიდიან")

    # ნაბიჯი 2
    print("\n--- ნაბიჯი 2: chi ველის რეაგირების განტოლება ---")
    chi_t, chi_short, chi_long, green_k = step2_chi_field_equation()
    print(f"  Green Kernel G(t,t') = {green_k}")
    print(f"  chi(t) = {chi_t}")
    print(f"  t << tau_chi-ში: chi(t) ≈ {chi_short}")
    print("  -> მეხსიერების რეჟიმში chi ინარჩუნებს თავდაპირველ კინემატიკურ კავშირს")

    # ნაბიჯი 3
    print("\n--- ნაბიჯი 3: Bullet Cluster რიცხვობრივი მასშტაბები ---")
    v_c, L_s, t_cross, tau_g, ratio = step3_bullet_cluster_timescales()
    print(f"  v_collision  = {v_c} km/s")
    print(f"  L_separation = {L_s} kpc")
    print(f"  t_cross      = {t_cross:.0f} Myr")
    print(f"  tau_chi (საილუსტრაციო) = {tau_g} Gyr")
    print(f"  t_coll/tau_chi = {ratio:.3f}  (<< 1 -> chi ინარჩუნებს პოზიციას)")

    # ნაბიჯი 4
    print("\n--- ნაბიჯი 4: ლენზირების პოტენციალის Toy Model შედარება ---")
    Phi_gas, Phi_gal = step4_toy_lensing_peak_comparison()
    print(f"  Phi_at_gas_center      = {Phi_gas}")
    print(f"  Phi_at_galaxy_center   = {Phi_gal}")
    print("  -> ლენზირების პიკი ფიქსირდება გალაქტიკის ცენტრში (თავსებადია დაკვირვებასთან)")

    # ნაბიჯი 5
    print("\n--- ნაბიჯი 5: MOND-ის კლასიკური მოლოდინი vs RFG-ის chi-მოდელი ---")
    contrast = step5_mond_contrast()
    for k, v in contrast.items():
        print(f"  {k:18s}: {v}")

    # ნაბიჯი 6
    print("\n--- ნაბიჯი 6: ფალსიფიცირებადი პროგნოზები ---")
    pred = step6_falsifiable_predictions()
    for k, v in pred.items():
        print(f"  {k:30s}: {v}")

    print("\n" + "=" * 72)
    print("სტატუსი და აგენტთა საბჭოს შენიშვნების დადასტურება:")
    print("=" * 72)
    print("""
    1. სტატუსი: Toy Model; სრული validation საჭიროებს κ(x,y), shear map, gas profile,
       galaxy distribution, projected mass და relativistic potentials შედარებას.
    2. Green-ის ფუნქციის ნიშანი: e^(-(t-t')/tau_chi) კოდში დაფიქსირდა,
       თუმცა ცვლადი წყაროს ინტეგრალური ფორმა აქ არ მოწმდება; გამოიყენება
       constant-source toy solution.
    3. chi ველის წყარო: კოდში დაემატა k_chi * omega_v (ვორტიკალური დინება)
       როგორც სიმბოლური ჰიპოთეზა/toy-წყარო.
    4. placeholder-ტექსტები სრულად წაიშალა მთელი ფაილიდან.
    """)
    print("=" * 72)
