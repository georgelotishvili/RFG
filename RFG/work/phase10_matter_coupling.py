# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp

def analyze_matter_coupling_tov():
    r = sp.Symbol('r', real=True, positive=True)
    
    # ველები
    Phi = sp.Function('Phi')(r) # ნიუტონის პოტენციალი
    rho_solid = sp.Function('rho_solid')(r) # სუპერსოლიდის ენერგიის სიმკვრივე
    
    # სუპერსოლიდის სტრეს-ტენზორის კომპონენტები
    p_rad = sp.Function('p_rad')(r)
    p_tan = sp.Function('p_tan')(r)
    delta_p = p_tan - p_rad # ანიზოტროპია
    
    # სრული ინერციული სიმკვრივე ნიუტონურ ლიმიტში: rho_inert = rho_solid + p_rad
    rho_inert = rho_solid + p_rad
    
    # სრული სისტემის ენერგია-იმპულსის შენახვის კანონი (\nabla_\mu T^{\mu 1} = 0)
    # სუსტ ველში (p << rho) და სტატიკურ სფერულ სიმეტრიაში გვაძლევს ანიზოტროპიულ TOV განტოლებას.
    # ეს აღწერს სითხის შიდა წონასწორობას და არა გარე ტესტ-ნაწილაკის აჩქარებას!
    # ნიუტონური ლიმიტი: d(p_rad)/dr + rho_inert * d(Phi)/dr - 2*delta_p/r = 0
    
    grad_Phi = sp.Symbol('grad_Phi') # პოტენციალის გრადიენტი d(Phi)/dr
    
    tov_eq = sp.Eq(sp.diff(p_rad, r) + rho_inert * grad_Phi - 2 * delta_p / r, 0)
    
    # ამოვხსნათ პოტენციალის გრადიენტისთვის შიდა წონასწორობაში
    sols = sp.solve(tov_eq, grad_Phi)
    if not sols:
        raise ValueError("ვერ მოიძებნა ამოხსნა grad_Phi-სთვის")
        
    g_sol = sols[0]
    return sp.simplify(g_sol)

if __name__ == "__main__":
    g_sol = analyze_matter_coupling_tov()
    print("--- ანიზოტროპიული TOV განტოლება ნიუტონურ ლიმიტში ---")
    print("პოტენციალის გრადიენტი (Phi') შიდა წონასწორობაში:")
    print(g_sol)
    
    print("\n--- კავშირი ტესტ-ნაწილაკის დინამიკასთან ---")
    print("ეს განტოლება აღწერს მედიუმის შიდა წონასწორობას.")
    print("გარე ტესტ-ნაწილაკის გეოდეზიური აჩქარება (a_test = -Phi') და MOND-ის")
    print("სრული ამოხსნა გაგრძელებულია phase11_mond_metric.py-ში.")
    
    print("\n--- აგენტთა საბჭოს შენიშვნები ---")
    print("1. ინერციული სიმკვრივე გასწორდა: rho_inert = rho_solid + p_rad.")
    print("2. განტოლება წარმოადგენს TOV-ის ნიუტონურ ლიმიტს; M-R წირის ასაგებად აკლია EoS და საზღვრის პირობები.")
    print("3. Phi'-ის ფორმულა არის წრიული (იმპლიციტური ODE), რადგან p_rad და rho_solid")
    print("   თვითონ არიან Phi-ს ფუნქციები ლაგრანჟიანიდან. ეს მოითხოვს მკაცრ მიბმას/fine-tuning-ს.")