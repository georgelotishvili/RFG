# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 39: Action-level route to the Z9 closure lattice

Status:
    Candidate derivation layer. This file does not yet prove the full
    charged-lepton theorem from the RFG action, but it removes one
    arbitrary-looking assumption from phase38:

        why C3 x C3?

Result:
    The two C3 factors can be traced to standard structures already present
    in an elastic supersolid:

    1. Axis C3:
       A framed/chiral defect carries an oriented material triad. The full
       permutation group S3 of three axes is reduced to the orientation-
       preserving subgroup A3, which is isomorphic to C3.

       Odd permutations flip the sign of the oriented volume form
       epsilon_ABC and therefore move the defect to the opposite chirality.

    2. Strain-phase C3:
       The traceless triaxial strain sector is a two-component order
       parameter. In complex form

           E = q1 + omega q2 + omega^2 q3,     q1+q2+q3=0,
           omega = exp(2*pi*i/3).

       A cyclic axis permutation sends E -> omega^2 E. Therefore E^3 is
       the first phase-sensitive invariant. This creates three discrete
       strain/braid phase sectors.

    Together:

        C3(axis orientation) x C3(strain phase) -> Z9 reduced closure slots.

    phase41 shows that the strain-phase C3 lock is not an external
    addition: it is the cubic invariant already contained in I3=det(B).
    phase40 then selects h=2 as the first non-trivial oriented framed
    branch, giving theta = h/9 = 2/9 in phase38.
"""

import itertools
import math

import sympy as sp


def permutation_parity(perm):
    """Return +1 for even permutations and -1 for odd permutations."""
    inversions = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inversions += 1
    return 1 if inversions % 2 == 0 else -1


def permutation_cycles(perm):
    """Human-readable cycle label for a permutation of (0,1,2)."""
    seen = [False] * len(perm)
    cycles = []
    for i in range(len(perm)):
        if seen[i]:
            continue
        cycle = []
        j = i
        while not seen[j]:
            seen[j] = True
            cycle.append(j + 1)
            j = perm[j]
        if len(cycle) > 1:
            cycles.append("(" + " ".join(str(x) for x in cycle) + ")")
    return "".join(cycles) if cycles else "identity"


def classify_axis_permutations():
    """
    S3 permutations of the elastic material triad.

    The oriented volume form epsilon_ABC keeps only even permutations in a
    fixed-chirality sector. The even set has three elements: C3.
    """
    rows = []
    for perm in itertools.permutations((0, 1, 2)):
        parity = permutation_parity(perm)
        rows.append(
            {
                "perm": perm,
                "cycle": permutation_cycles(perm),
                "det": parity,
                "orientation_preserving": parity == 1,
                "group_sector": "A3 ~= C3" if parity == 1 else "opposite chirality",
            }
        )
    return rows


def axis_c3_order():
    return sum(1 for row in classify_axis_permutations() if row["orientation_preserving"])


def strain_doublet_transform():
    """
    Verify that the triaxial complex strain E transforms by a C3 phase.

    E = q1 + omega q2 + omega^2 q3.
    Under the cyclic permutation (q1,q2,q3)->(q2,q3,q1),
    E -> omega^2 E. Hence E^3 is invariant.
    """
    q1, q2, q3 = sp.symbols("q1 q2 q3")
    omega = sp.exp(2 * sp.pi * sp.I / 3)
    E = q1 + omega * q2 + omega**2 * q3
    E_cyclic = q2 + omega * q3 + omega**2 * q1
    transform_ratio = sp.simplify(E_cyclic / E)

    # SymPy may not reduce exp(2*pi*i/3) expressions completely in a ratio,
    # so verify by coefficient comparison.
    expected = omega**2 * E
    diff = sp.simplify(sp.expand(E_cyclic - expected))

    invariant_diff = sp.simplify(sp.expand(E_cyclic**3 - E**3))

    return {
        "E": E,
        "E_cyclic": E_cyclic,
        "expected_transform": "E -> omega^2 E",
        "transform_difference": diff,
        "E3_invariant_difference": invariant_diff,
        "phase_sector_order": 3,
    }


def z9_closure_from_action_symmetry():
    """Combine the two C3 factors into the Z9 closure lattice."""
    axis_order = axis_c3_order()
    phase_order = strain_doublet_transform()["phase_sector_order"]
    return {
        "axis_c3_origin": "oriented elastic triad: S3 -> A3 ~= C3",
        "axis_order": axis_order,
        "phase_c3_origin": "triaxial strain doublet: E -> omega^2 E, E^3 invariant",
        "phase_order": phase_order,
        "closure_slots": axis_order * phase_order,
        "theta_if_h2": 2 / (axis_order * phase_order),
    }


def allowed_phase_locking_terms(max_power=9):
    """
    List low-order powers E^n and whether they are invariant under the
    cyclic strain phase E -> omega^2 E.
    """
    rows = []
    for n in range(1, max_power + 1):
        invariant = (2 * n) % 3 == 0
        rows.append(
            {
                "power": n,
                "phase_factor": f"omega^{2*n}",
                "c3_invariant": invariant,
                "meaning": "allowed anisotropy" if invariant else "forbidden by C3",
            }
        )
    return rows


def spinorial_selection_statement():
    """
    This is the remaining assumption chain, made explicit.
    """
    return {
        "h0": "trivial closure; gives a partially degenerate spectrum, not three charged generations",
        "h1": "projective/nematic closure only, n -> -n; see phase40",
        "h2": "first non-trivial oriented framed closure, n -> n; see phase40",
        "higher_h": "on the small positive C3 branch, h>=3 crosses the zero-eigenfrequency edge",
        "remaining_theorem": "derive the oriented-frame charged coupling directly from the RFG action",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 39: Action-level route to C3 x C3 -> Z9")
    print("=" * 72)

    print("\n1. Axis permutations of the oriented elastic triad")
    for row in classify_axis_permutations():
        print(
            f"  perm={row['perm']} cycle={row['cycle']:<9} "
            f"det={row['det']:+d} sector={row['group_sector']}"
        )
    print(f"  axis C3 order = {axis_c3_order()}")

    print("\n2. Triaxial strain doublet")
    strain = strain_doublet_transform()
    for key, value in strain.items():
        print(f"  {key:26s}: {value}")

    print("\n3. Allowed C3 phase-locking powers")
    for row in allowed_phase_locking_terms():
        print(
            f"  E^{row['power']:<2} -> {row['phase_factor']:<8} "
            f"{row['meaning']}"
        )

    print("\n4. Z9 closure")
    closure = z9_closure_from_action_symmetry()
    for key, value in closure.items():
        if isinstance(value, float):
            print(f"  {key:20s}: {value:.12f}")
        else:
            print(f"  {key:20s}: {value}")

    print("\n5. Remaining spinorial selection")
    for key, value in spinorial_selection_statement().items():
        print(f"  {key:20s}: {value}")
