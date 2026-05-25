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
phase37_c3_koide_operator.py ამოწმებს ალტერნატიულ structural candidate-ს,
სადაც charged-lepton sqrt(m) სპექტრი მოდის C3-ციკლური ოპერატორიდან:
nu_k = A[1 + sqrt(2) cos(2/9 + 2*pi*k/3)].
phase38_z9_theta_holonomy.py ამატებს theta=2/9-ის candidate derivation-ს.
"""
import sympy as sp
import math
try:
    import scipy.special as sc
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    from phase37_c3_koide_operator import (
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
    m_e_pdg = 0.51099895
    m_mu_pdg = 105.658375
    m_tau_pdg = 1776.86
    
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

def calculate_koide(b_e, b_mu, b_tau):
    m_e = 0.51099895
    m_mu = 105.658375
    m_tau = 1776.86
    
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
    print("\n--- Koide-ს 2/3 ფარდობა და ნეიტრინოს პროგნოზი ---")
    print(f"თეორიული სამიზნე (Koide limit): K = 2/3 ≈ 0.666667")
    print(f"PDG ექსპერიმენტული მასებით:   K_exp = {K_exp:.6f} (ცდომილება {err_exp:.4f}%)")
    print(f"RFG სრული Mathieu-თი (q={q_val}):   K_rfg = {K_rfg:.6f} (ცდომილება {err_rfg:.4f}%)")
    print("\nწინასწარმეტყველება ნეიტრინოებისთვის (3 თაობა):")
    print("ჰიპოთეზური პროგნოზი: K_nu ≈ 2/3, შესამოწმებელია დამოუკიდებლად.")

    print("\n--- C3 Koide operator შედარება (phase37) ---")
    c3 = calculate_c3_comparison()
    if c3 is None:
        print("phase37_c3_koide_operator.py ვერ ჩაიტვირთა.")
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
        print("theta=2/9-ის candidate derivation იხ. phase38; ღიაა მისი RFG-action-იდან გამოყვანა.")

    print("\n--- აგენტთა საბჭოს შენიშვნები / შეზღუდვები ---")
    print("1. ფარული fit: N=72 და 295 შერჩეულია √m_i/m_e-ს მიხედვით. 0 DoF = 0 პრედიქცია. (დამოუკიდებელი ტოპოლოგიური წესი ღიაა).")
    print("2. A_param და q_param-ში omega^2 განზომილების გაყოფა სრულად გასწორდა.")
    print("3. b_N გამოითვლება ნამდვილად scipy.special-ით N^2 მიახლოების ნაცვლად, და hardcoded 0.0006% ამოღებულია.")
    print("4. ფაილი წარმოადგენს რიცხვობრივ კანდიდატურას (numerical candidate) და არა ტოპოლოგიურ მტკიცებულებას.")
    print("5. ამოცანა განიხილება როგორც დროითი სექტორის მიახლოება (time-sector toy model), R'(r) გრადიენტი იგნორირებულია.")
