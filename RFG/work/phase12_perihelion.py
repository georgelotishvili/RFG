# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
დაკვირვებითი სანდოობის ტესტი (Sanity-check).
ეს ფაილი დამოუკიდებლად არ ამტკიცებს RFG მეტრიკას. ის ამოწმებს, რომ 
phase8_weak_field.py-დან მიღებული gamma=1 და beta=1 იძლევა მერკურის სწორ პრეცესიას.
"""
import math

def calculate_mercury_precession():
    # გრავიტაციული და ფიზიკური მუდმივები (SI ერთეულებში)
    G_val = 6.67430e-11
    M_sun = 1.98847e30
    c_val = 299792458.0
    
    # მერკურის ორბიტის პარამეტრები
    a_mercury = 57.90905e9   # მეტრი (დიდი ნახევარღერძი)
    e_mercury = 0.205630     # ექსცენტრისიტეტი
    
    # კეპლერის პერიოდი T = 2π√(a³/GM), ფარული ემპირიული შენატანის თავიდან ასაცილებლად
    T_mercury_sec = 2 * math.pi * math.sqrt(a_mercury**3 / (G_val * M_sun))
    T_mercury_days = T_mercury_sec / (24.0 * 3600.0)
    days_per_century = 36525.0
    
    # PPN პარამეტრები გამოყვანილი RFG თეორიიდან (იხ. phase8_weak_field.py)
    gamma_val = 1.0
    beta_val = 1.0
    
    # თეორიული PPN ფაქტორი: (2 + 2*gamma - beta) / 3
    ppn_factor = (2 + 2 * gamma_val - beta_val) / 3.0
    
    # პრეცესია თითო ორბიტაზე (რადიანებში)
    delta_phi_rad = (6 * math.pi * G_val * M_sun) / (c_val**2 * a_mercury * (1 - e_mercury**2)) * ppn_factor
    
    # გადაყვანა არკწამებში თითო საუკუნეზე
    orbits_per_century = days_per_century / T_mercury_days
    rad_to_arcsec = (180.0 / math.pi) * 3600.0
    precession_arcsec_per_century = delta_phi_rad * orbits_per_century * rad_to_arcsec
    
    return ppn_factor, precession_arcsec_per_century

if __name__ == "__main__":
    ppn_factor, precession = calculate_mercury_precession()
    print("--- მერკურის პერიჰელიონის პრეცესია RFG თეორიაში ---")
    print(f"PPN ფაქტორი ((2 + 2*gamma - beta)/3): {ppn_factor}")
    print(f"გამოთვლილი პრეცესია: {precession:.2f} არკწამი/საუკუნეში")
    print("დაკვირვებული მნიშვნელობა (GR): 42.98 არკწამი/საუკუნეში")
    assert abs(precession - 42.98) < 0.1, f"ცდომილება დიდია: პრეცესია = {precession}"
    print("დასკვნა: პრეცესია ზუსტად ემთხვევა დაკვირვებად 42.98″/cy მნიშვნელობას (assert გავლილია).")