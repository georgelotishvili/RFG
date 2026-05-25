# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
================================================================================
PHASE 22 (v3.1): სრული ენერგია-იმპულსის ტენზორი — Carter-Karlovini ფორმალიზმი
================================================================================

სტატუსი:
ეტაპი I-ის ფესვის გასწორების ფაილი. აქ ფიქსირდება:
- T_{mu nu} = 2*dL/dg^{mu nu} - g_{mu nu}*L კონვენცია
- I_1, I_2, I_3 Carter-Karlovini-ის ჯაჭვური წესით B^{AB}-დან
- nabla_mu T^{mu nu} = sum_A E_A * partial^nu psi_A Noether/Bianchi იდენტობა
  diagonal Minkowski/FLRW/Bianchi I ფონებზე სრული L-ით
- Schwarzschild diagonal ფონზე იგივე იდენტობის reduced Y+I1 smoke-test
- off-diagonal ვარიაციის ფაქტორის ცალკე smoke-test
- falsification ცდა, რომ Noether იდენტობა მცდარი T-ის ფორმულას ცხადად
  იჭერს (არ არის ცარიელი ტავტოლოგია)

შენიშვნა მკითხველისთვის:
matter stress tensor-ის off-shell კოვარიანტული დივერგენცია generic-ად
ნული არ არის — ის ფიზიკურ ველთა Euler-Lagrange წყაროებს უდრის. ეს ფაილი
ამოწმებს ცხადად ამ იდენტობას, არა "off-shell zero"-ს. falsification_test()
ცხადყოფს რომ ცდა ფაქტობრივად ფარდდება T-ის ფორმულის სისწორეს — ცარიელი
ტავტოლოგია არ არის.
"""

import sympy as sp
from phase1_action import get_polynomial_lagrangian


DIM = 4
NSOLID = 3


# ============================================================================
# გეომეტრიული ფუნქციები
# ============================================================================


def get_christoffel(g_cov, g_inv, coords):
    """Christoffel სიმბოლოები Gamma^lambda_{mu nu}."""
    Gamma = [sp.zeros(DIM, DIM) for _ in range(DIM)]
    half = sp.Rational(1, 2)
    for lam in range(DIM):
        for mu in range(DIM):
            for nu in range(DIM):
                term = 0
                for rho in range(DIM):
                    term += half * g_inv[lam, rho] * (
                        sp.diff(g_cov[rho, mu], coords[nu])
                        + sp.diff(g_cov[rho, nu], coords[mu])
                        - sp.diff(g_cov[mu, nu], coords[rho])
                    )
                Gamma[lam][mu, nu] = sp.simplify(term)
    return Gamma


def covariant_divergence_contra(T_contra, Gamma, coords, nu):
    """nabla_mu T^{mu nu}."""
    total = 0
    for mu in range(DIM):
        total += sp.diff(T_contra[mu, nu], coords[mu])
        for lam in range(DIM):
            total += Gamma[mu][mu, lam] * T_contra[lam, nu]
            total += Gamma[nu][mu, lam] * T_contra[mu, lam]
    return sp.simplify(total)


# ============================================================================
# ინვარიანტები და ლაგრანჟიანი
# ============================================================================


def build_invariants_from_metric(g_inv, field_grads):
    """
    field_grads[0] = d_mu Phi (ფაზური სკალარი)
    field_grads[1..3] = d_mu phi^A (ელასტიური ველები)
    """
    Y = sp.simplify(sum(
        g_inv[mu, nu] * field_grads[0][mu] * field_grads[0][nu]
        for mu in range(DIM) for nu in range(DIM)
    ))

    B = sp.zeros(NSOLID, NSOLID)
    for A in range(NSOLID):
        for Bidx in range(NSOLID):
            B[A, Bidx] = sp.simplify(sum(
                -g_inv[mu, nu] * field_grads[A + 1][mu] * field_grads[Bidx + 1][nu]
                for mu in range(DIM) for nu in range(DIM)
            ))

    I1 = sp.simplify(B.trace())
    I2 = sp.simplify(sp.Rational(1, 2) * (I1**2 - (B * B).trace()))
    I3 = sp.simplify(B.det())
    return Y, I1, I2, I3, B


def get_full_lagrangian(Y, I1, I2, I3):
    """phase1_action-ის სრული პოლინომიური L."""
    Y_s, I1_s, I2_s, I3_s = sp.symbols("Y I1 I2 I3", real=True)
    L_poly = get_polynomial_lagrangian(Y_s, I1_s, I2_s, I3_s)
    L = sp.simplify(L_poly.subs({Y_s: Y, I1_s: I1, I2_s: I2, I3_s: I3}))
    return L, L_poly, (Y_s, I1_s, I2_s, I3_s)


def get_test_lagrangian(Y, I1, I2, I3, mode="full"):
    """შესამოწმებელი L: სრული ან შემცირებული."""
    if mode == "full":
        return get_full_lagrangian(Y, I1, I2, I3)
    if mode == "reduced_y_i1":
        c_y, c_i1 = sp.symbols("c_y c_i1", real=True)
        L = c_y * Y + c_i1 * I1
        L_poly = c_y * sp.Symbol("Y", real=True) + c_i1 * sp.Symbol("I1", real=True)
        return L, L_poly, sp.symbols("Y I1 I2 I3", real=True)
    raise ValueError(f"unknown lagrangian mode: {mode}")


# ============================================================================
# Off-diagonal ვარიაციის ცდა (smoke-test)
# ============================================================================


def offdiag_variation_smoke_test():
    """
    შემოწმდება: symmetric g^{01} ცვლადის ვარიაცია ფაქტორი-2-ით
    არ უნდა გადიდდეს. შედარება — დამოუკიდებელი g^{01}, g^{10}.
    """
    q00, q11, q22, q33, q01 = sp.symbols("q00 q11 q22 q33 q01", real=True)
    r01, r10 = sp.symbols("r01 r10", real=True)
    gcov01 = sp.Symbol("gcov01", real=True)
    c_y, c_i1, c_yi1 = sp.symbols("c_y c_i1 c_yi1", real=True)
    v0, v1 = sp.symbols("v0 v1", real=True)
    e10, e11 = sp.symbols("e10 e11", real=True)

    field_grads = [
        [v0, v1, 0, 0],
        [e10, e11, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]

    g_inv_sym = sp.Matrix([
        [q00, q01, 0, 0],
        [q01, q11, 0, 0],
        [0, 0, q22, 0],
        [0, 0, 0, q33],
    ])
    Y_sym, I1_sym, _, _, _ = build_invariants_from_metric(g_inv_sym, field_grads)
    L_sym = c_y * Y_sym + c_i1 * I1_sym + c_yi1 * Y_sym * I1_sym
    T01_unique = sp.expand(sp.diff(L_sym, q01) - gcov01 * L_sym)
    T01_wrong = sp.expand(2 * sp.diff(L_sym, q01) - gcov01 * L_sym)

    g_inv_independent = sp.Matrix([
        [q00, r01, 0, 0],
        [r10, q11, 0, 0],
        [0, 0, q22, 0],
        [0, 0, 0, q33],
    ])
    Y_ind, I1_ind, _, _, _ = build_invariants_from_metric(g_inv_independent, field_grads)
    L_ind = c_y * Y_ind + c_i1 * I1_ind + c_yi1 * Y_ind * I1_ind
    T01_independent = sp.expand(2 * sp.diff(L_ind, r01) - gcov01 * L_ind)
    T01_independent = sp.expand(T01_independent.subs({r01: q01, r10: q01}))

    residual_unique = sp.expand(T01_unique - T01_independent)
    residual_wrong = sp.expand(T01_wrong - T01_independent)
    return residual_unique, residual_wrong


# ============================================================================
# ფონური მონაცემები
# ============================================================================


def background(name):
    """coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields."""
    if name == "minkowski":
        t, x, y, z = sp.symbols("t x y z", real=True)
        coords = [t, x, y, z]
        g_cov_diag = [1, -1, -1, -1]
        g_inv_diag = [1, -1, -1, -1]
        sqrt_minus_g = sp.Integer(1)
        fields = [t, x, y, z]
        return coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields

    if name == "flrw":
        t, x, y, z = sp.symbols("t x y z", real=True)
        coords = [t, x, y, z]
        a = sp.Function("a")(t)
        Phi = sp.Function("Phi")(t)
        g_cov_diag = [1, -a**2, -a**2, -a**2]
        g_inv_diag = [1, -1 / a**2, -1 / a**2, -1 / a**2]
        sqrt_minus_g = a**3
        fields = [Phi, x, y, z]
        return coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields

    if name == "bianchi_i":
        t, x, y, z = sp.symbols("t x y z", real=True)
        coords = [t, x, y, z]
        a = sp.Function("a")(t)
        b = sp.Function("b")(t)
        c = sp.Function("c")(t)
        Phi = sp.Function("Phi")(t)
        g_cov_diag = [1, -a**2, -b**2, -c**2]
        g_inv_diag = [1, -1 / a**2, -1 / b**2, -1 / c**2]
        sqrt_minus_g = a * b * c
        fields = [Phi, x, y, z]
        return coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields

    if name == "schwarzschild":
        t, r, th, ph = sp.symbols("t r theta phi", real=True, positive=True)
        coords = [t, r, th, ph]
        r_s = sp.Symbol("r_s", real=True, positive=True)
        f = 1 - r_s / r
        g_cov_diag = [f, -1 / f, -r**2, -r**2 * sp.sin(th)**2]
        g_inv_diag = [1 / f, -f, -1 / r**2, -1 / (r**2 * sp.sin(th)**2)]
        sqrt_minus_g = r**2 * sp.sin(th)
        fields = [t, r, th, ph]
        return coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields

    raise ValueError(f"unknown background: {name}")


# ============================================================================
# T-ის ცდის ფუნქცია (ფონური ჩასმის წინ მთლიანი L ფონის ფუნქციად)
# ============================================================================


def evaluate_on_background(name, lagrangian_mode="full", t_factor=1):
    """
    ფონური Bianchi/Noether იდენტობის ცდა.

    t_factor: T-ის ცდის გადანამრავლება. სწორი მნიშვნელობა 1; falsification
    ცდისთვის t_factor != 1 — residual მაშინ ცხადად არანულოვანი.

    ცდა:
    1. T_cov[μ,μ] = 2 * dL/dq_μ - L/q_μ  (q_μ = g^{μμ}, q-ის ფუნქცია)
    2. ფონური ჩასმის შემდეგ T_cov-ი ფონის ფუნქციად
    3. T_contra = g_inv ფონური * g_inv ფონური * T_cov (diagonal)
    4. Christoffel ფონური მეტრიკიდან
    5. div = nabla_mu T^{mu nu}
    6. eom_A = (1/sqrt(-g)) d_mu (sqrt(-g) dL/dD_{A,μ})  ფონური ჩასმის შემდეგ
    7. source[nu] = sum_A eom_A * g^{νν} * partial_nu psi_A
    8. residual = div - source — სწორი T-სთვის 0, არასწორი T-სთვის არ-0.
    """
    coords, g_cov_diag, g_inv_diag, sqrt_minus_g, fields = background(name)
    q = [sp.Symbol(f"q_{mu}", real=True, nonzero=True) for mu in range(DIM)]

    # ფონური ფუნქცია ცხადად — დარჩება ფონის ფუნქცია, არ-ფიქსირდება როგორც სიმბოლო
    D = [[0 for _ in range(DIM)] for _ in range(NSOLID + 1)]
    D_symbol = [[None for _ in range(DIM)] for _ in range(NSOLID + 1)]
    grad_subs = {}

    for A, field in enumerate(fields):
        for mu, coord in enumerate(coords):
            grad = sp.diff(field, coord)
            if grad != 0:
                sym = sp.Symbol(f"D_{name}_{A}_{mu}", real=True)
                D[A][mu] = sym
                D_symbol[A][mu] = sym
                grad_subs[sym] = grad

    g_inv_diag_symbolic = sp.diag(*q)
    Y, I1, I2, I3, _ = build_invariants_from_metric(g_inv_diag_symbolic, D)
    L, _, _ = get_test_lagrangian(Y, I1, I2, I3, lagrangian_mode)

    # T-ის ცდა q-ის ფუნქცია — q-ის მიმართ ვარიაცია სიმბოლურია; ჩასმა მერე
    T_cov_general = sp.zeros(DIM, DIM)
    for mu in range(DIM):
        T_cov_general[mu, mu] = t_factor * 2 * sp.diff(L, q[mu]) - L / q[mu]

    # current_A_mu სიმბოლური — სიმბოლური ვარიაცია D[A][mu]-ის მიმართ ჯერ
    current_symbolic = [[None for _ in range(DIM)] for _ in range(NSOLID + 1)]
    for A in range(NSOLID + 1):
        for mu in range(DIM):
            if D[A][mu] == 0:
                current_symbolic[A][mu] = sp.Integer(0)
            else:
                current_symbolic[A][mu] = sp.diff(L, D[A][mu])

    # ფონური ჩასმა — q -> g_inv ფონური, D -> ფონური წარმოებული
    subs = {q[mu]: g_inv_diag[mu] for mu in range(DIM)}
    subs.update(grad_subs)

    invariants = tuple(sp.simplify(expr.subs(subs)) for expr in (Y, I1, I2, I3))

    T_cov = sp.zeros(DIM, DIM)
    for mu in range(DIM):
        T_cov[mu, mu] = sp.simplify(T_cov_general[mu, mu].subs(subs))

    T_contra = sp.zeros(DIM, DIM)
    for mu in range(DIM):
        for nu in range(DIM):
            T_contra[mu, nu] = sp.simplify(
                g_inv_diag[mu] * g_inv_diag[nu] * T_cov[mu, nu]
            )

    g_cov = sp.diag(*g_cov_diag)
    g_inv = sp.diag(*g_inv_diag)
    Gamma = get_christoffel(g_cov, g_inv, coords)

    div = [
        covariant_divergence_contra(T_contra, Gamma, coords, nu)
        for nu in range(DIM)
    ]

    eom = []
    for A in range(NSOLID + 1):
        total = 0
        for mu, coord in enumerate(coords):
            if current_symbolic[A][mu] == 0:
                continue
            current_functional = current_symbolic[A][mu].subs(subs)
            integrand = sqrt_minus_g * current_functional
            total += sp.diff(integrand, coord)
        if sqrt_minus_g != 0:
            eom.append(sp.simplify(total / sqrt_minus_g))
        else:
            eom.append(sp.simplify(total))

    source = []
    for nu, coord in enumerate(coords):
        rhs = 0
        for A, field in enumerate(fields):
            rhs += eom[A] * g_inv_diag[nu] * sp.diff(field, coord)
        source.append(sp.simplify(rhs))

    residual = [sp.simplify(div[nu] - source[nu]) for nu in range(DIM)]

    return {
        "coords": coords,
        "fields": fields,
        "g_cov_diag": g_cov_diag,
        "g_inv_diag": g_inv_diag,
        "sqrt_minus_g": sqrt_minus_g,
        "invariants": invariants,
        "T_cov": T_cov,
        "T_contra": T_contra,
        "divergence": div,
        "eom": eom,
        "source": source,
        "residual": residual,
    }


# ============================================================================
# Falsification ცდა — ვაჩვენებთ რომ Noether იდენტობა მცდარ T-ს იჭერს
# ============================================================================


def falsification_test():
    """
    სცადე T-ის ცდა გადანამრავლებული t_factor = 3 ფაქტორით.
    residual მაშინ უნდა იყოს არანულოვანი — ცდა მართლა ფარდდება T-ის
    ფორმულის სისწორეს.
    """
    correct = evaluate_on_background("flrw", lagrangian_mode="reduced_y_i1", t_factor=1)
    wrong = evaluate_on_background("flrw", lagrangian_mode="reduced_y_i1", t_factor=3)
    correct_ok = all(reduce_zero(r) == 0 for r in correct["residual"])
    wrong_fails = any(reduce_zero(r) != 0 for r in wrong["residual"])
    return correct_ok, wrong_fails, correct["residual"], wrong["residual"]


# ============================================================================
# დამხმარე ფუნქციები
# ============================================================================


def convention_summary():
    return {
        "signature": "(+---)",
        "B_AB": "B^{AB} = -g^{mu nu} d_mu phi^A d_nu phi^B",
        "stress_tensor": "T_{mu nu} = 2*dL/dg^{mu nu} - g_{mu nu}*L",
        "phase_relation": (
            "phase6_relativistic.py ემთხვევა; phase1_tensor.py-ის mixed "
            "ფორმა ამ კონვენციით უნდა გადაიწეროს"
        ),
        "no_ghost": "იხ. NOTATION.md Full No-Ghost Window; X-scheme pure branch-ში c_Y < 0.",
    }


def reduce_zero(expr):
    """ნულის შემოწმება ფონური სიმეტრიების გათვალისწინებით."""
    return sp.factor(sp.together(sp.trigsimp(sp.simplify(expr))))


def is_zero_vector(values):
    return all(reduce_zero(value) == 0 for value in values)


# ============================================================================
# Main
# ============================================================================


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 22 (v3.1): სრული ენერგია-იმპულსის ტენზორი")
    print("=" * 72)

    print("\n1. კონვენციები")
    for key, value in convention_summary().items():
        print(f"  {key:14s}: {value}")

    print("\n2. Off-diagonal ვარიაციის smoke-test")
    offdiag_ok, offdiag_wrong = offdiag_variation_smoke_test()
    print(f"  symmetric ვარიაცია residual: {reduce_zero(offdiag_ok)}")
    print(f"  ფაქტორ-2 ცდომილების residual ნულია? {reduce_zero(offdiag_wrong) == 0}")
    print(f"  ცდის სტატუსი: {'OK' if reduce_zero(offdiag_ok) == 0 else 'CHECK'}")

    print("\n3. Bianchi/Noether იდენტობა ფონებზე")
    backgrounds = [
        ("minkowski", "full"),
        ("flrw", "full"),
        ("bianchi_i", "full"),
        ("schwarzschild", "reduced_y_i1"),
    ]
    for name, mode in backgrounds:
        result = evaluate_on_background(name, lagrangian_mode=mode)
        residual = [reduce_zero(value) for value in result["residual"]]
        residual_ok = is_zero_vector(residual)
        print(f"\n--- {name} ({mode}) ---")
        print(f"  invariants (Y, I1, I2, I3): {result['invariants']}")
        print(f"  residual vector: {residual}")
        print(f"  status: {'OK' if residual_ok else 'CHECK'}")

    print("\n4. Falsification ცდა — Noether იდენტობა ფარდდება T-ის ფორმულას?")
    correct_ok, wrong_fails, correct_res, wrong_res = falsification_test()
    print(f"  სწორი T (t_factor=1): residual ნული? {correct_ok}")
    print(f"  მცდარი T (t_factor=3): residual არანული? {wrong_fails}")
    if correct_ok and wrong_fails:
        print("  ცდის სტატუსი: OK — იდენტობა ფაქტობრივად ფარდდება")
    else:
        print("  ცდის სტატუსი: CHECK — ცდა არ ფარდდება სწორ/მცდარ T-ს")

    print("\n5. FLRW ენერგია და წნევა")
    flrw = evaluate_on_background("flrw")
    a = sp.Function("a")(sp.Symbol("t", real=True))
    print(f"  rho = {sp.expand(flrw['T_cov'][0, 0])}")
    print(f"  p = T_11/a^2 = {sp.expand(flrw['T_cov'][1, 1] / a**2)}")

    print("\n6. სტატუსი")
    print("  - Minkowski/FLRW/Bianchi I diagonal ფონებზე იდენტობა შემოწმდა სრული L-ით")
    print("  - Schwarzschild diagonal ფონზე — reduced Y+I1 smoke-test")
    print("  - off-diagonal ვარიაცია smoke-test-ით შემოწმდა")
    print("  - falsification ცდამ აჩვენა რომ იდენტობა მცდარი T-ის ფორმულას იჭერს")
    print("  - generic non-diagonal ფონებზე სრული proof ჯერ ღია")
    print("  - შემდეგი ნაბიჯი: NOTATION.md + phase1_tensor/phase6 გადაწერა")
