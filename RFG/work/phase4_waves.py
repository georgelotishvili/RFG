# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
ლოკალური მოდულუსების შემოწმების ესკიზი.

შენიშვნა:
- ეს ფაილი არ არის ტალღების სრული მტკიცება და არ შეიცავს დიაგონალიზებულ eigenmode-ებს.
- სრული ტალღებისთვის (2x2 კინეტიკური/გრადიენტული მატრიცა) გამოიყენება phase1_action.py.
- გრავიტაციული ტალღებისთვის გამოიყენება phase9_gravitational_waves.py (TT მეტრიკის ექსპანსია).
- H_YY არ უნდა ჩაითვალოს გრავიტაციული ტალღების სიხისტედ; c_T მოდის აინშტაინ-ჰილბერტის სექტორიდან.
"""
import sympy as sp
from phase1_action import init_variables, get_polynomial_lagrangian

def get_local_moduli():
    Y, I1, I2, I3 = init_variables()
    L = get_polynomial_lagrangian(Y, I1, I2, I3)
    
    # ლოკალური მოდულუსებისთვის ვითვლით ლაგრანჟიანის მეორე რიგის წარმოებულებს (ჰესიანს)
    hessian_YY = sp.diff(L, Y, 2)
    hessian_I1I1 = sp.diff(L, I1, 2)
    hessian_YI1 = sp.diff(L, Y, I1)
    
    # ფონზე შეფასება (Minkowski: Y=1, I1=3, I2=3, I3=1)
    bg_subs = {Y: 1, I1: 3, I2: 3, I3: 1}
    
    H_YY_bg = sp.simplify(hessian_YY.subs(bg_subs))
    H_I1I1_bg = sp.simplify(hessian_I1I1.subs(bg_subs))
    H_YI1_bg = sp.simplify(hessian_YI1.subs(bg_subs))
    
    # სრული 2x2 ჰესიანის დეტერმინანტი
    H_det_bg = sp.simplify(H_YY_bg * H_I1I1_bg - H_YI1_bg**2)
    
    return H_YY_bg, H_I1I1_bg, H_YI1_bg, H_det_bg

if __name__ == "__main__":
    h_YY, h_I1I1, h_YI1, h_det = get_local_moduli()
    print("ფაზური ლოკალური მოდულუსი ფონზე (H_YY):", h_YY)
    print("ელასტიური ლოკალური მოდულუსი ფონზე (H_I1I1):", h_I1I1)
    print("შერეული კავშირის მოდულუსი ფონზე (H_YI1):", h_YI1)
    
    print("\n--- ლოკალური სტაბილურობის პირობები (მოდულუსების დადებითობა) ---")
    print(f"დიაგონალური ფაზური სტაბილურობა: H_YY > 0  =>  {h_YY} > 0")
    print(f"დიაგონალური ელასტიური სტაბილურობა: H_I1I1 > 0  =>  {h_I1I1} > 0")
    print(f"სრული 2x2 ჰესიანის დადებითობა (Det > 0): {h_det} > 0")
    print("\nეს აჩვენებს მხოლოდ ლოკალური დიაგონალური და შერეული მოდულუსების დადებითობის პირობებს;")
    print("სრული ტალღური სტაბილურობა და eigenmode-ები მოწმდება phase1_action.py-ში.")