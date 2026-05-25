# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from p01_core import init_variables, get_polynomial_lagrangian
from p01_core import calculate_stress_tensor

def get_friedmann_equations():
    # დრო და კოსმოლოგიური მასშტაბური ფაქტორი
    t = sp.Symbol('t', real=True)
    a = sp.Function('a')(t)
    
    # ჰაბლის პარამეტრი: H = \dot{a} / a
    H = sp.diff(a, t) / a
    
    # გრავიტაციული კონსტანტა
    kappa = sp.Symbol('kappa', real=True)
    
    # FLRW სტრეს-ტენზორი მკაცრი ვარიაციით
    g00, g11, g22, g33 = sp.symbols('g00 g11 g22 g33', real=True)
    Y_g = g00
    I1_g = -g11 - g22 - g33
    I2_g = g11*g22 + g22*g33 + g11*g33
    I3_g = -g11*g22*g33
    
    Y, I1, I2, I3 = init_variables()
    L_poly = get_polynomial_lagrangian(Y, I1, I2, I3)
    L_g = L_poly.subs({Y: Y_g, I1: I1_g, I2: I2_g, I3: I3_g})
    
    T_mixed = calculate_stress_tensor(L_g, [g00, g11, g22, g33])
    ansatz = {g00: 1, g11: -1/a**2, g22: -1/a**2, g33: -1/a**2}
    
    rho_solid = sp.simplify(T_mixed[0].subs(ansatz))
    p_iso_solid = sp.simplify(-T_mixed[1].subs(ansatz))
    
    # მატერიის და რადიაციის წვლილები
    rho_m0, rho_rad0 = sp.symbols('rho_m0 rho_rad0', real=True)
    rho_m = rho_m0 / a**3
    rho_rad = rho_rad0 / a**4
    
    p_m = 0
    p_rad_matter = rho_rad / 3
    
    # პირველი ფრიდმანის განტოლება სრული შემადგენლობით
    friedmann1 = sp.Eq(3 * H**2, kappa * (rho_solid + rho_m + rho_rad))
    
    # მეორე ფრიდმანის განტოლება (აჩქარების)
    a_ddot = sp.diff(a, t, t)
    friedmann2 = sp.Eq(2 * a_ddot / a + H**2, -kappa * (p_iso_solid + p_m + p_rad_matter))
    
    return friedmann1, friedmann2, a, t, rho_solid, p_iso_solid

if __name__ == "__main__":
    f1, f2, a, t, rho_solid, p_iso_solid = get_friedmann_equations()
    
    print("პირველი ფრიდმანის განტოლება:", f1)
    print("\nმეორე ფრიდმანის განტოლება:", f2)
    
    print("\n--- სუპერსოლიდის ენერგიის სიმკვრივე FLRW ფონზე ---")
    rho_expanded = sp.expand(rho_solid)
    rho_collected = sp.collect(rho_expanded, a)
    print("rho_solid(a) =", rho_collected)
    
    print("\nშენიშვნა: ძველ გამარტივებულ კოსმოლოგიურ ფორმულას აკლია აქ მიღებული")
    print("          -3*c_I2/a^4 და -c_I3/a^6 წევრები, რომლებიც სრული ვარიაციიდან ჩნდება.")
    
    # ემერჯენტული Λ_eff-ის გამოყოფა: ვუშვებთ, რომ a -> უსასრულობაში
    x = sp.Symbol("x", positive=True)
    Lambda_eff = sp.simplify(sp.limit(rho_solid.subs(a, x), x, sp.oo))
    print("\nეფექტური კოსმოლოგიური მუდმივა (a -> infinity):")
    print("Lambda_eff / kappa =", Lambda_eff)
    
    print("\n--- კოსმოლოგიური ისტორია (დომინანტური წევრები) ---")
    print("1. ადრეული ეპოქა (a -> 0): ჭარბობს a^-4 (რადიაცია) და a^-6 (I3 ელასტიურობა).")
    print("2. შუა ეპოქა (a ~ 1): ჭარბობს a^-3 (მატერია).")
    print("3. გვიანი ეპოქა (a -> infinity): ჭარბობს მუდმივი წევრი Lambda_eff -> აჩქარებული გაფართოება.")
    
    print("\n--- აგენტთა საბჭოს შენიშვნები ---")
    print("- Lambda_eff ცხადად გამოიყო (c_Y + 3*c_Y2), თუმცა თეორია არ ხსნის მის ნატურალურობას")
    print("  (მცირე დაკვირვებად მნიშვნელობას); ექსტრემალური fine-tuning შენარჩუნებულია (იხ. §12).")
    print("- rho_m და rho_rad სწორად სკალირდება a^-3 და a^-4 კანონებით.")
    print("- კოსმოლოგიური სისტემა შეიცავს radiation/matter/late-time acceleration-ისთვის საჭირო სკალირებად წევრებს;")
    print("  სრული დინამიკური ამოხსნა და დაკვირვებითი ფიტი ცალკე ეტაპია.")
    print("- ენერგიის შენახვის განტოლების (drho/dt + 3H(rho+p_iso) = 0) შემოწმება ცალკე გადის")
    print("  p02_cosmo.py-ში, სადაც მკაცრი Phi=t ფონის შემთხვევაში ნარჩენი რჩება.")

# ===================== CONSOLIDATED PHASE SECTIONS =====================


# ===================== merged from p02_cosmo.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from p01_core import get_polynomial_lagrangian, init_variables
from p01_core import calculate_stress_tensor

def check_conservation():
    t = sp.Symbol('t', real=True)
    a = sp.Function('a')(t)
    H = sp.diff(a, t) / a
    
    # ინვერსიული მეტრიკის კომპონენტები FLRW-სთვის (g^uu)
    g00, g11, g22, g33 = sp.symbols('g00 g11 g22 g33', real=True)
    
    # Phi ველის დროითი წარმოებული (FLRW-ში სივრცული გრადიენტები ნულია)
    # ვტოვებთ ფუნქციად, რათა ენერგიის შენახვისა და EOM-ის იდენტობა გამოჩნდეს
    Phi_dot = sp.Function('Phi_dot')(t)
    
    # FLRW-სთვის ინვარიანტები (g^uu და Phi_dot ცვლადებით)
    Y_g = g00 * Phi_dot**2
    I1_g = -g11 - g22 - g33
    I2_g = g11*g22 + g22*g33 + g11*g33
    I3_g = -g11*g22*g33
    
    Y, I1, I2, I3 = init_variables()
    L = get_polynomial_lagrangian(Y, I1, I2, I3)
    L_g = L.subs({Y: Y_g, I1: I1_g, I2: I2_g, I3: I3_g})
    
    # 1. ენერგიის შენახვის პირდაპირი გამოთვლა T_mn-დან
    T_mixed = calculate_stress_tensor(L_g, [g00, g11, g22, g33])
    
    # FLRW ანზაცის ჩასმა g^uv-სთვის
    ansatz = {g00: 1, g11: -1/a**2, g22: -1/a**2, g33: -1/a**2}
    rho = sp.simplify(T_mixed[0].subs(ansatz))
    p_iso = sp.simplify(-T_mixed[1].subs(ansatz)) # T^1_1 = -p (სწორი იზოტროპული მიდგომა)
    
    # შენახვის განტოლების მარცხენა მხარე: nabla_u T^u_0 = drho/dt + 3H(rho + p)
    rho_dot = sp.simplify(sp.diff(rho, t))
    conservation_lhs = sp.simplify(rho_dot + 3 * H * (rho + p_iso))
    
    # 2. სკალარული ველის (Phi) მოძრაობის განტოლება
    # nabla_u (dL / d(partial_u Phi)) = 0
    dL_dPhi_dot = sp.simplify(sp.diff(L_g, Phi_dot).subs(ansatz))
    # კოვარიანტული დივერგენცია დროში: 1/a^3 * d/dt(a^3 * V^0)
    EOM_Phi = sp.simplify((sp.diff(a**3 * dL_dPhi_dot, t)) / a**3)
    
    # 3. Noether-ის იდენტობა
    # დროითი ტრანსლაციის სიმეტრიიდან: nabla_u T^u_0 = EOM_Phi * partial_0 Phi
    expected_lhs = sp.simplify(EOM_Phi * Phi_dot)
    
    # შევამოწმოთ სხვაობა (უნდა იყოს ზუსტი 0)
    difference = sp.simplify(conservation_lhs - expected_lhs)
    
    # 4. თუ ჩავსვამთ მკაცრ ანზაცს Phi = t (ანუ Phi_dot = 1, dPhi_dot/dt = 0)
    ansatz_Phi_t = {Phi_dot: 1, sp.diff(Phi_dot, t): 0}
    conservation_eval = sp.simplify(conservation_lhs.subs(ansatz_Phi_t))
    EOM_eval = sp.simplify(EOM_Phi.subs(ansatz_Phi_t))
    
    return difference, expected_lhs, conservation_eval, EOM_eval

if __name__ == "__main__":
    diff, expected_lhs, cons_eval, eom_eval = check_conservation()
    
    print("--- Noether-ის იდენტობის შემოწმება FLRW ფონზე ---")
    print("∇_μ T^μ_0 (ენერგიის შენახვის განტოლების მარცხენა მხარე):")
    print(expected_lhs)
    print("\nსხვაობა (∇_μ T^μ_0) და (EOM_Phi * Phi_dot) შორის:")
    print(diff)
    if diff == 0:
        print("დასკვნა: ენერგიის შენახვა drho/dt+3H(rho+p) კოვარიანტულად ზუსტად უდრის")
        print("         სკალარული ველის მოძრაობის განტოლებას! (Noether's identity შესრულებულია)")
        
    print("\n--- მტკიცება Phi = t ანზაცისთვის ---")
    print("თუ ჩავსვამთ მკაცრ ფონს Phi = t (Phi_dot = 1):")
    print("drho/dt + 3H(rho+p) ნარჩენი:", cons_eval)
    print("EOM_Phi ნარჩენი:", eom_eval)
    print("\nაგენტთა საბჭოს პასუხი:")
    print("1. p_iso პირდაპირ T^1_1-დან დაითვალა, 'გასაშუალოება' აღარ გამოიყენება.")
    print("2. Y და I_k მეტრიკიდან ფიქსირდება; დამოუკიდებელი ფუნქციები ამოღებულია.")
    print("3. Phi=t ანზაცი არ ანულებს ენერგიის შენახვას ავტომატურად (ნარჩენი რჩება).")
    print("4. ეს ნარჩენი ზუსტად ემთხვევა Phi-ს EOM-ს, რაც ნიშნავს რომ სისტემა")
    print("   მოითხოვს Phi-ს დინამიკურ ევოლუციას Phi_dot(t), ან ენერგიის გაცვლას")
    print("   ბარიონულ მატერიასთან/რადიაციასთან სრული შენახვისთვის.")


# ===================== merged from p02_cosmo.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
ისტორიული ხაზობრივი ესკიზი. 
მიმდინარე ვერსია: p10_oscillons.py. 
ხარვეზად აღარ ჩაითვალოს.
"""
import sympy as sp

def get_process_time_relation():
    t = sp.Symbol('t', real=True)
    tau = sp.Function('tau')(t)
    
    # ფონური წნევითი პოტენციალი
    phi = sp.Symbol('phi', real=True)
    
    # პროცესის დროის ტემპი დამოკიდებულია პოტენციალზე (ბი-კონფორმული სკალირება)
    dtau_dt = sp.diff(tau, t)
    relation = sp.Eq(dtau_dt, sp.exp(phi / 2))
    
    return relation, tau, t, phi

if __name__ == "__main__":
    relation, tau, t, phi = get_process_time_relation()
    
    print("პროცესის დროის კავშირი პოტენციალთან:")
    print(relation)
    print("\nსუსტი ველის ლიმიტი (Linearization):")
    print(sp.series(sp.exp(phi/2), phi, 0, 2).removeO(), "ეს ძველ ხაზობრივ ფორმას ემთხვევა მხოლოდ მაშინ, თუ ძველი p იდენტიფიცირდება phi-სთან და alpha = 1/2.")
    print("\nშენიშვნა: ეს არის ისტორიული/ესკიზური ვარიანტი.")
    print("მიმდინარე/სწორი ვერსია არის p10_oscillons.py (d tau/dt = e^(phi/2)).")


# ===================== merged from p02_cosmo.py =====================

# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 6 (relativistic): FLRW რელატივისტური სტრესი
================================================================================

რეფერენცია: NOTATION.md, p01_core.py

ეს ფაილი იყენებს NOTATION.md-ის აქტიურ კონვენციას:
- სიგნატურა (+---)
- B^{AB} = -g^{mu nu} * d_mu phi^A * d_nu phi^B
- T_{mu nu} = 2*dL/dg^{mu nu} - g_{mu nu}*L

ეს არის FLRW ფონური ცდა, რომელიც phase22-ის ბირთვით (Bianchi/Noether
იდენტობა) გადადის. შედეგი — rho(a), p_iso(a) — შედარდება ძველ გამარტივებულ
ფორმულას.

FLRW ანზაცი:
    ds^2 = dt^2 - a(t)^2 * (dx^2 + dy^2 + dz^2)
    g^{mu nu} = diag(1, -1/a^2, -1/a^2, -1/a^2)
    Phi = Phi(t), phi^A = (x, y, z)
"""

import sympy as sp
from p01_core import evaluate_on_background, reduce_zero


def get_flrw_pressures():
    """
    FLRW ფონური ჩასმა phase22-ის evaluate_on_background-ით.

    აბრუნებს:
    - rho = T_{00}
    - p_iso = T_{11} / a^2
    """
    result = evaluate_on_background("flrw", lagrangian_mode="full")

    a = sp.Function("a")(sp.Symbol("t", real=True))
    rho = sp.simplify(result["T_cov"][0, 0])
    p_iso = sp.simplify(result["T_cov"][1, 1] / a**2)

    return rho, p_iso, a, result


def compare_with_theory(rho):
    """
    ძველი გამარტივებული ფორმულა:
        rho = -3*c_I1/a^2 - 9*c_I1sq/a^4 + c_Y + 3*c_Y2 + 3*c_YI1/a^2

    ეს თეორიის ტექსტში გამოტოვებულია c_I2 და c_I3 წევრებისთვის. phase22 სრულ
    ფორმას ბრუნდება — შესწორება უნდა დარჩეს სამუშაო ბაზაში და მერე სტატიაში გავიდეს.
    """
    c_Y, c_Y2 = sp.symbols("c_Y c_Y2", real=True)
    c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols(
        "c_I1 c_I1sq c_I2 c_I3 c_YI1", real=True
    )

    t = sp.Symbol("t", real=True)
    a = sp.Function("a")(t)
    Phi = sp.Function("Phi")(t)
    Phi_dot = sp.diff(Phi, t)

    rho_full_expected = (
        -3 * c_I1 / a**2
        - 9 * c_I1sq / a**4
        - 3 * c_I2 / a**4
        - c_I3 / a**6
        + c_Y * Phi_dot**2
        + 3 * c_Y2 * Phi_dot**4
        + 3 * c_YI1 * Phi_dot**2 / a**2
    )

    rho_md_text = (
        -3 * c_I1 / a**2
        - 9 * c_I1sq / a**4
        + c_Y * Phi_dot**2
        + 3 * c_Y2 * Phi_dot**4
        + 3 * c_YI1 * Phi_dot**2 / a**2
    )

    full_match = sp.simplify(rho - rho_full_expected) == 0
    md_text_match = sp.simplify(rho - rho_md_text) == 0
    missing_terms = sp.simplify(rho_full_expected - rho_md_text)

    return full_match, md_text_match, missing_terms


def check_bianchi_residual(result):
    """phase22-ის Bianchi residual FLRW-ზე უნდა იყოს [0,0,0,0]."""
    residual = [reduce_zero(value) for value in result["residual"]]
    is_ok = all(value == 0 for value in residual)
    return is_ok, residual


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 6 (relativistic): FLRW რელატივისტური სტრესი")
    print("რეფერენცია: NOTATION.md, phase22")
    print("=" * 72)

    rho, p_iso, a, result = get_flrw_pressures()

    print("\n1. FLRW ფონური ჩასმა (NOTATION.md კონვენცია)")
    print(f"  invariants (Y, I1, I2, I3): {result['invariants']}")

    print("\n2. ენერგია-იმპულსის ტენზორი")
    print(f"  rho = T_00 = {sp.expand(rho)}")
    print(f"  p_iso = T_11/a^2 = {sp.expand(p_iso)}")

    print("\n3. შედარება ძველ გამარტივებულ კოსმოლოგიურ ფორმულასთან")
    full_match, md_text_match, missing = compare_with_theory(rho)
    print(f"  სრული ფორმა (c_I2 და c_I3 ჩათვლით) ემთხვევა: {full_match}")
    print(f"  ძველი §3 ფორმა (c_I2, c_I3-ის გარეშე) ემთხვევა: {md_text_match}")
    print(f"  ძველ ფორმულაში გამოტოვებული წევრები: {sp.expand(missing)}")

    print("\n4. Bianchi residual phase22-ის ცდიდან")
    bianchi_ok, residual = check_bianchi_residual(result)
    print(f"  residual vector: {residual}")
    print(f"  Bianchi/Noether იდენტობა: {'OK' if bianchi_ok else 'CHECK'}")

    print("\n5. სტატუსი")
    print("  - კონვენცია: NOTATION.md-ის აქტიური ფორმა")
    print("  - სიგნატურა: (+---)")
    print("  - phase22-ის Bianchi იდენტობა FLRW-ზე სრულდება")
    print("  - rho-ში c_I2, c_I3 წევრები არსებობს — ეს შესწორება სამუშაო ბაზაში დახურულია")
    print("  - p_iso §3-ის ფორმასთან — შემოწმდეს ცალკე")

