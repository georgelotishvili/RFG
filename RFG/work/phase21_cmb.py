# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 21: CMB — კოსმოლოგიური პერტურბაციები და EFT 
================================================================================

სტატუსი:
ეს ფაილი წარმოადგენს კონცეპტუალურ ხიდს კოსმოლოგიური პერტურბაციების 
ეფექტური ველის თეორიასთან (EFT). სრული CMB სპექტრის (C_l), BAO-სა და 
ლენზირების გამოსათვლელად საჭიროა განტოლებების CLASS/CAMB კოდებში 
ინტეგრირება, რაც ცალკე ამოცანაა.

ცენტრალური დასკვნები:
    1. სტანდარტული Horndeski (Bellini-Sawicki) პარამეტრიზაცია არ არის 
       საკმარისი RFG-სთვის, რადგან I_k ელასტიური სექტორი შეიცავს 
       სივრცულ ტრანსვერსულ მოდებს. საჭიროა "EFT of Solid Inflation" (ESS) ჩარჩო.
    2. c_Y კოეფიციენტი უნდა იყოს უარყოფითი (c_Y < 0), რათა კინეტიკური 
       წევრი (alpha_K) იყოს დადებითი და ავირიდოთ Ghost არასტაბილურობა.
    3. გრავიტაციული ტალღის სიჩქარე c_T = c, მაგრამ მასიური დისპერსიის 
       ასარიდებლად საჭიროა phase9-ის კონსტრეინტი.

References:
    - Bellini & Sawicki 2014, JCAP 07:050
    - Endlich, Nicolis, Wang 2013 (Solid Inflation)
"""

import sympy as sp
from sympy import symbols, Symbol, simplify, diff


def map_rfg_to_horndeski():
    """
    RFG-ის L_solid → Horndeski G_2, G_3, G_4, G_5 (მხოლოდ Y-სექტორით)

    Y = -2X (Bellini-Sawicki convention: X = -½ g^μν ∂_μφ ∂_νφ)
    ⟹ Y-სექტორი ⊂ G_2(X, φ)

    I_k-სექტორი არ ჯდება სტანდარტულ Horndeski-ში — ის მოითხოვს
    EFT of Solid Inflation (ESS) ჩარჩოს.

    Horndeski (მხოლოდ φ) სექტორში:
        G_2 = c_Y·Y + c_Y2·Y² = -2c_Y·X + 4·c_Y2·X²
        G_3 = 0   (no kinetic mixing)
        G_4 = M_Pl²/2
        G_5 = 0
    """
    X = Symbol('X', real=True)
    # აგენტთა საბჭოს შესწორება: c_Y არ უნდა იყოს positive=True, 
    # რადგან alpha_K > 0 მოითხოვს c_Y < 0-ს.
    c_Y = Symbol('c_Y', real=True) 
    c_Y2, M_Pl = symbols('c_Y2 M_Pl', positive=True)

    G_2 = -2 * c_Y * X + 4 * c_Y2 * X**2
    G_3 = sp.Integer(0)
    G_4 = M_Pl**2 / 2
    G_5 = sp.Integer(0)

    return G_2, G_3, G_4, G_5, X


# ==============================================================================
# ნაბიჯი 1: α_T (tensor speed excess)
# ==============================================================================

def compute_alpha_T():
    """
    Bellini-Sawicki:
        α_T = 2X·(G_{4,X} - G_{5,φ}) / M_*²    +  (G_5,X-term)

    RFG-ში:
        G_4 = M_Pl²/2 (X-დამოუკიდებელი) ⟹ G_{4,X} = 0
        G_5 = 0 ⟹ G_{5,φ} = 0

    ⟹ α_T = 0  (c_T = c)

    თუმცა I_k სექტორი წარმოშობს გრავიტონის ეფექტურ მასას (phase9).
    მკაცრი GW170817 თავსებადობისთვის მოითხოვება phase9 კონსტრეინტი:
    -0.5*c_Y - 0.5*c_Y2 + 0.5*c_I1 + 7.5*c_I1sq + 1.5*c_I2 + 0.5*c_I3 + 0.5*c_YI1 = 0
    """
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()
    M_star_sq = Symbol('M_star_sq', positive=True)
    phi = Symbol('phi', real=True)

    # G_{4,X} = 0 since G_4 const in X
    G4_X = diff(G_4, X)
    # G_{5,φ}, G_{5,X} = 0
    G5_phi = diff(G_5, phi)

    # α_T = 2X(G_{4,X} - G_{5,φ})/M_*²
    alpha_T = 2*X*(G4_X - G5_phi) / M_star_sq
    alpha_T = simplify(alpha_T)

    return alpha_T, G4_X, G5_phi


# ==============================================================================
# ნაბიჯი 2: α_M (Planck-mass running)
# ==============================================================================

def compute_alpha_M():
    """
    Bellini-Sawicki:
        M_*² = 2(G_4 - 2X·G_{4,X} + ...)

    α_M = (d ln M_*² / dt) / H
    რადგან G_4 = M_Pl²/2 = const, საბაზისო Horndeski სექტორში α_M = 0.
    """
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()
    M_Pl = Symbol('M_Pl', positive=True)

    M_star_sq = 2 * G_4  # = M_Pl²
    M_star_sq_simplified = simplify(M_star_sq)

    # dM_*²/dt = 0 (M_Pl const)
    t = Symbol('t', real=True)
    M_star_sq_t = M_star_sq_simplified  # no t-dependence
    d_M_star_dt = diff(M_star_sq_t, t)

    alpha_M = sp.Integer(0)

    return alpha_M, M_star_sq_simplified, d_M_star_dt


# ==============================================================================
# ნაბიჯი 3: α_B (braiding)
# ==============================================================================

def compute_alpha_B():
    """
    Bellini-Sawicki:
        α_B = 2(X·G_{3,X}·φ̇/H ·... + G_{4,X}·... ) / M_*²

    α_B = 0, რადგან G_3 = 0 და G_4 = const.
    """
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()

    G3_X = diff(G_3, X)
    G4_X = diff(G_4, X)

    # All zero ⟹ α_B = 0
    alpha_B = sp.Integer(0)

    return alpha_B, G3_X, G4_X


# ==============================================================================
# ნაბიჯი 4: α_K (kineticity)
# ==============================================================================

def compute_alpha_K():
    """
    Bellini-Sawicki:
        α_K = (2X·G_{2,X} + 4X²·G_{2,XX} + ...) / (H²·M_*²)

    G_2 = -2c_Y·X + 4c_Y2·X²
    α_K = (-4c_Y·X + 48c_Y2·X²) / (H²·M_Pl²)
    სტაბილურობა (no-ghost) მოითხოვს α_K > 0.
    ვინაიდან X დადებითია (time-like დერივატივი), c_Y უნდა იყოს უარყოფითი!
    """
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()
    H, M_Pl = symbols('H M_Pl', positive=True)

    G2_X = diff(G_2, X)
    G2_XX = diff(G_2, X, 2)

    alpha_K = (2*X*G2_X + 4*X**2*G2_XX) / (H**2 * M_Pl**2)
    alpha_K = simplify(alpha_K)

    return alpha_K, G2_X, G2_XX


# ==============================================================================
# ნაბიჯი 5: CMB სპექტრის თავსებადობა
# ==============================================================================

def cmb_consistency_check():
    return {
        'alpha_T': '0 (Y-sector; needs phase9 constraint: c_I1/2 + 15*c_I1sq/2 + 3*c_I2/2 + c_I3/2 - c_Y/2 - c_Y2/2 + c_YI1/2 = 0)',
        'alpha_M': '0 (მხოლოდ Y/Horndeski სექტორში; I_k მოითხოვს ESS ანალიზს)',
        'alpha_B': '0 (მხოლოდ Y/Horndeski სექტორში; I_k მოითხოვს ESS ანალიზს)',
        'alpha_K': '~ X·c_Y2/M_Pl² (small if c_Y2 ≪ M_Pl²/X_bg)',
        'CMB_spectrum': 'not tested here; requires CLASS/CAMB C_l computation',
        'ISW_BAO_lensing': 'not tested here; requires BAO/lensing likelihood comparison',
    }


# ==============================================================================
# ნაბიჯი 6: I_k სექტორის წვლილი ფონზე
# ==============================================================================

def i_k_sector_on_flrw():
    """
    I_k სექტორი FLRW ფონზე:
        ρ(I_1)    ∝ 1/a²   ← curvature-like
        ρ(I_1²)   ∝ 1/a⁴   ← radiation-like
        ρ(I_2)    ∝ 1/a⁴   ← radiation-like
        ρ(I_3)    ∝ 1/a⁶   ← stiff fluid

    BBN (დიდი აფეთქების ნუკლეოსინთეზის) შეზღუდვები:
    რადიაციის მსგავსი წევრები (c_I1sq, c_I2) არ უნდა აჭარბებდნენ 
    დასაშვებ ეფექტურ ნეიტრინოთა რაოდენობას (ΔN_eff).
    """
    a = Symbol('a', positive=True)

    rho_I1 = sp.Symbol('c_I1') * 3 / a**2          # 1/a² (curvature-like)
    rho_I1sq = sp.Symbol('c_I1sq') * 9 / a**4      # 1/a⁴ (radiation-like)
    rho_I2 = sp.Symbol('c_I2') * 3 / a**4          # 1/a⁴
    rho_I3 = sp.Symbol('c_I3') * 1 / a**6          # 1/a⁶ (stiff)

    return rho_I1, rho_I1sq, rho_I2, rho_I3


# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 21: CMB — ეფექტური ველის თეორიის (EFT) კავშირები")
    print("=" * 72)

    print("\n--- Horndeski მაპირება ---")
    G_2, G_3, G_4, G_5, X = map_rfg_to_horndeski()
    print(f"  G_2(X) = {G_2}")
    print(f"  G_3    = {G_3}")
    print(f"  G_4    = {G_4}")
    print(f"  G_5    = {G_5}")

    print("\n--- ნაბიჯი 1: α_T (tensor speed excess) ---")
    aT, G4X, G5p = compute_alpha_T()
    print(f"  G_{{4,X}} = {G4X}")
    print(f"  G_{{5,φ}} = {G5p}")
    print(f"  α_T = {aT} Y/Horndeski ქვე-სექტორში; სრული solid-sector GW თავსებადობა დამატებით მოითხოვს phase9-ის მასის კონსტრეინტს.")

    print("\n--- ნაბიჯი 2: α_M (Planck-mass running) ---")
    aM, M_star_sq, dM_dt = compute_alpha_M()
    print(f"  M_*² = 2·G_4 = {M_star_sq}")
    print(f"  dM_*²/dt = {dM_dt}")
    print(f"  α_M = {aM}. ეს მხოლოდ Y/Horndeski ქვე-სექტორშია, ხოლო I_k solid sector მოითხოვს ESS/full perturbation ანალიზს.")

    print("\n--- ნაბიჯი 3: α_B (braiding) ---")
    aB, G3X, G4X2 = compute_alpha_B()
    print(f"  G_{{3,X}} = {G3X}")
    print(f"  G_{{4,X}} = {G4X2}")
    print(f"  α_B = {aB}. ეს მხოლოდ Y/Horndeski ქვე-სექტორშია, ხოლო I_k solid sector მოითხოვს ESS/full perturbation ანალიზს.")

    print("\n--- ნაბიჯი 4: α_K (kineticity) ---")
    aK, G2X, G2XX = compute_alpha_K()
    print(f"  G_{{2,X}}  = {G2X}")
    print(f"  G_{{2,XX}} = {G2XX}")
    print(f"  α_K = {aK}")
    print(f"  Ghost-ის თავიდან ასაცილებლად საჭიროა α_K > 0.")
    print(f"  აქედან გამომდინარეობს კრიტიკული პირობა: c_Y < 0.")

    print("\n--- ნაბიჯი 5: CMB თავსებადობა ---")
    check = cmb_consistency_check()
    for k, v in check.items():
        print(f"  {k:18s}: {v}")

    print("\n--- ნაბიჯი 6: I_k სექტორი და BBN ლიმიტები ---")
    r1, r1sq, r2, r3 = i_k_sector_on_flrw()
    print(f"  ρ(I_1)    = {r1}    (∝ 1/a²)")
    print(f"  ρ(I_1²)   = {r1sq}  (∝ 1/a⁴)")
    print(f"  ρ(I_2)    = {r2}    (∝ 1/a⁴)")
    print(f"  ρ(I_3)    = {r3}    (∝ 1/a⁶)")
    print("  BBN ლიმიტი: c_I1sq და c_I2 ≲ ΔN_eff * ρ_gamma.")

    # შემაჯამებელი
    print("\n" + "=" * 72)
    print("აგენტთა საბჭოს შენიშვნების დადასტურება:")
    print("1. placeholder ტექსტები სრულად გასუფთავდა.")
    print("2. c_Y-ის positive=True დეკლარაცია მოიხსნა. α_K-ს გაანალიზებამ აჩვენა,")
    print("   რომ No-Ghost პირობა მოითხოვს c_Y < 0 (აგენტ-მათემატიკოსის სწორი შენიშვნა).")
    print("3. Bellini-Sawicki პარამეტრიზაციის არასრულფასოვნება RFG-სთვის აღიარებულია.")
    print("   ელასტიური სექტორისთვის აუცილებელია 'EFT of Solid Inflation' (ESS) ჩარჩო.")
    print("4. CMB-ის სპექტრის (C_l), BAO-ს და Lensing-ის ფორმალური რიცხვობრივი მორგება")
    print("   მოითხოვს CLASS/CAMB კოდებში იმპლემენტაციას (Future Work).")
    print("5. GW მასის კონსტრეინტი phase9-დან პირდაპირ იქნა ციტირებული α_T ბლოკში.")
    print("=" * 72)