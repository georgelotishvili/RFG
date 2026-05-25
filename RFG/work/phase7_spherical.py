# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

import sympy as sp
from phase1_action import get_polynomial_lagrangian

def solve_static_spherical():
    r = sp.Symbol('r', real=True, positive=True)
    rs = sp.Symbol('rs', real=True, positive=True)
    eps = sp.Symbol('eps', real=True)
    kappa = sp.Symbol('kappa', real=True)
    
    # უცნობი ფუნქციები (სიგნატურით + - - -)
    A = sp.Function('A')(r)
    B = sp.Function('B')(r)
    Psi_p = sp.Function('Psi_p')(r) # ეს არის Psi'(r), სადაც Phi = t + Psi(r)
    
    # აინშტაინის ტენზორი G^t_t და G^r_r მეტრიკისთვის diag(B, -A, -r^2, -r^2 sin^2 theta)
    G_tt = -sp.diff(A, r) / (r * A**2) + (1/A - 1)/r**2
    G_rr = sp.diff(B, r) / (r * A * B) + (1/A - 1)/r**2
    G_thth = sp.diff(B, r, 2)/(2*A*B) - sp.diff(B, r)**2/(4*A*B**2) - sp.diff(A, r)*sp.diff(B, r)/(4*A**2*B) + sp.diff(B, r)/(2*r*A*B) - sp.diff(A, r)/(2*r*A**2)
    
    # ინვარიანტები (Phi = t + Psi(r) და phi^A = x^A comoving ანზაცისთვის)
    Y = 1/B - Psi_p**2 / A
    I1 = 2 + 1/A
    I2 = 1 + 2/A
    I3 = 1/A
    
    Y_s, I1_s, I2_s, I3_s = sp.symbols('Y I1 I2 I3', real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    
    L_eval = L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_Y = sp.diff(L_poly, Y_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I1 = sp.diff(L_poly, I1_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I2 = sp.diff(L_poly, I2_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    L_I3 = sp.diff(L_poly, I3_s).subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3})
    
    # სტრეს-ტენზორი T^mu_nu
    T_tt = 2 * L_Y / B - L_eval
    T_rr = 2 * (L_Y * Psi_p**2 / A + L_I1 / A + 2 * L_I2 / A + L_I3 / A) - L_eval
    T_thth = 2 * (L_I1 + L_I2 * (1 + 1/A) + L_I3 / A) - L_eval
    
    # სკალარული ველის განტოლება: \nabla_\mu ( dL / d(\partial_\mu \Phi) ) = 0
    scalar_eq = sp.Derivative(r**2 * sp.sqrt(B/A) * L_Y * Psi_p, r)
    
    # სუსტი ველის ექსპანსია (ვუშვებთ Psi_p = 0 ცენტრალური მუხტის არარსებობის გამო)
    a1 = sp.Symbol('a1', real=True)
    b2 = sp.Symbol('b2', real=True)
    
    U = eps * rs / r
    A_w = 1 + a1 * U
    B_w = 1 - U + b2 * U**2
    
    # ვანაცვლებთ A და B ცვლადებს სუსტი ველის ფუნქციებით და ვითვლით G და T ტენზორებს
    G_tt_w = -sp.diff(A_w, r) / (r * A_w**2) + (1/A_w - 1)/r**2
    G_rr_w = sp.diff(B_w, r) / (r * A_w * B_w) + (1/A_w - 1)/r**2
    G_thth_w = sp.diff(B_w, r, 2)/(2*A_w*B_w) - sp.diff(B_w, r)**2/(4*A_w*B_w**2) - sp.diff(A_w, r)*sp.diff(B_w, r)/(4*A_w**2*B_w) + sp.diff(B_w, r)/(2*r*A_w*B_w) - sp.diff(A_w, r)/(2*r*A_w**2)
    
    Y_w = 1/B_w
    I1_w = 2 + 1/A_w
    I2_w = 1 + 2/A_w
    I3_w = 1/A_w
    
    L_w = L_poly.subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    L_Y_w = sp.diff(L_poly, Y_s).subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    L_I1_w = sp.diff(L_poly, I1_s).subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    L_I2_w = sp.diff(L_poly, I2_s).subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    L_I3_w = sp.diff(L_poly, I3_s).subs({Y_s: Y_w, I1_s: I1_w, I2_s: I2_w, I3_s: I3_w})
    
    T_tt_w = 2 * L_Y_w / B_w - L_w
    T_rr_w = 2 * (L_I1_w / A_w + 2 * L_I2_w / A_w + L_I3_w / A_w) - L_w
    T_thth_w = 2 * (L_I1_w + L_I2_w * (1 + 1/A_w) + L_I3_w / A_w) - L_w
    
    def get_series(expr):
        return sp.simplify(sp.series(expr, eps, 0, 3).removeO())
        
    Eq_tt_w = get_series(G_tt_w - kappa * T_tt_w)
    Eq_rr_w = get_series(G_rr_w - kappa * T_rr_w)
    Eq_thth_w = get_series(G_thth_w - kappa * T_thth_w)
    
    # გამოვყოთ ნულოვანი (ფონური) და პირველი რიგის განტოლებები
    Eq_tt_O0 = sp.simplify(Eq_tt_w.subs(eps, 0))
    Eq_tt_O1 = sp.simplify(sp.diff(Eq_tt_w, eps).subs(eps, 0))
    
    Eq_rr_O0 = sp.simplify(Eq_rr_w.subs(eps, 0))
    Eq_rr_O1 = sp.simplify(sp.diff(Eq_rr_w, eps).subs(eps, 0))
    
    Eq_thth_O0 = sp.simplify(Eq_thth_w.subs(eps, 0))
    Eq_thth_O1 = sp.simplify(sp.diff(Eq_thth_w, eps).subs(eps, 0))
    
    # ბი-კონფორმობის ანალიზი (წინასწარ a1=1 დაშვების გარეშე)
    Delta_T = sp.simplify(T_tt_w - T_rr_w)
    Delta_T_O1 = sp.simplify(sp.diff(Delta_T, eps).subs(eps, 0))
    
    c_Y, c_Y2, c_I1, c_I1sq, c_I2, c_I3, c_YI1 = sp.symbols('c_Y c_Y2 c_I1 c_I1sq c_I2 c_I3 c_YI1', real=True)
    E_vac = -c_Y - 3*c_Y2 + 3*c_I1 + 9*c_I1sq + 3*c_I2 + c_I3 - 3*c_YI1
    P_vac = c_Y + c_Y2 + 5*c_I1 + 21*c_I1sq + 7*c_I2 + 3*c_I3 + 5*c_YI1
    vac_sols = sp.solve([E_vac, P_vac], [c_Y, c_I1])
    
    clean_Delta_T_O1 = sp.simplify(Delta_T_O1.subs(vac_sols))
    bc_constraint = sp.Eq(clean_Delta_T_O1 / (rs/r), 0)
    
    return G_tt, G_rr, G_thth, T_tt, T_rr, T_thth, scalar_eq, Eq_tt_w, Eq_rr_w, Eq_thth_w, Eq_tt_O0, Eq_tt_O1, Eq_rr_O0, Eq_rr_O1, Eq_thth_O0, Eq_thth_O1, Delta_T_O1, clean_Delta_T_O1, bc_constraint

if __name__ == "__main__":
    res = solve_static_spherical()
    G_tt, G_rr, G_thth, T_tt, T_rr, T_thth, scalar_eq, Eq_tt_w, Eq_rr_w, Eq_thth_w, Eq_tt_O0, Eq_tt_O1, Eq_rr_O0, Eq_rr_O1, Eq_thth_O0, Eq_thth_O1, Delta_T_O1, clean_Delta_T_O1, bc_constraint = res
    print("--- ამოხსნა სფერული ანზაცისთვის ---")
    print("G^t_t =", G_tt)
    print("T^t_t =", T_tt)
    print("G^r_r =", G_rr)
    print("T^r_r =", T_rr)
    print("G^theta_theta =", G_thth)
    print("T^theta_theta =", T_thth)
    print("\nსკალარული ველის განტოლება:")
    print("0 =", scalar_eq)
    print("-> დასკვნა: ცენტრალური სკალარული მუხტის გარეშე Psi'(r) = 0, ანუ Phi(t,r) = t.")
    print("\n--- სუსტი ველის ლიმიტი (O(rs/r)) ---")
    print("Eq_tt (O(eps)):", Eq_tt_w)
    print("Eq_rr (O(eps)):", Eq_rr_w)
    print("Eq_thth (O(eps)):", Eq_thth_w)
    
    print("\n--- ფონური ვაკუუმის განტოლებები (O(1)) ---")
    print("Eq_tt_O0 =", Eq_tt_O0)
    print("Eq_rr_O0 =", Eq_rr_O0)
    print("Eq_thth_O0 =", Eq_thth_O0)
    
    print("\n--- პირველი რიგის განტოლებები (O(eps)) ---")
    print("Eq_tt_O1 =", Eq_tt_O1)
    print("Eq_rr_O1 =", Eq_rr_O1)
    print("Eq_thth_O1 =", Eq_thth_O1)

    print("\n--- ბი-კონფორმობის პირობა (g_tt * g_rr = -1) ---")
    print("Delta_T(O(eps)) =", Delta_T_O1)
    print("Delta_T(O(eps)) სუფთა (Minkowski ვაკუუმის კონსტრეინტებით) =", clean_Delta_T_O1)
    print("კონსტრეინტი ბი-კონფორმულობისთვის:", bc_constraint)
    
    print("\n--- აგენტთა საბჭოს დასკვნები ---")
    print("1. Angular კომპონენტი დაემატა; სრული სფერული სისტემის დახურვა საჭიროებს corrected O(eps) განტოლებების ერთობლივ ამოხსნას.")
    print("2. U=rs/r ანზაცით წარმოებულები სიმბოლურად ითვლება (eps-სერიების აღრევა აღმოიფხვრა).")
    print("3. T_rr-ის ნიშნის შეცდომა გასწორდა (T_rr და T_thth ახლა დადებითი ელასტიური წევრებით იწყება).")
    print("4. ბი-კონფორმობის ანალიზში a1=1 წინასწარ აღარ იდება. Minkowski E_vac=0, P_vac=0 კონსტრეინტებით მიიღება სუფთა შეზღუდვა.")