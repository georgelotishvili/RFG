# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
ისტორიული ხაზობრივი ესკიზი. 
მიმდინარე ვერსია: phase17_effective_mass.py. 
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
    print("მიმდინარე/სწორი ვერსია არის phase17_effective_mass.py (d tau/dt = e^(phi/2)).")