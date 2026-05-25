# Notation header (see NOTATION.md):
# signature (+---); Y = g^mn d_m Phi d_n Phi; B^AB = -g^mn d_m phi^A d_n phi^B.
# T_mn = 2*dL/dg^mn - g_mn*L; off-diagonal symmetric variables use factor 1.
# Horndeski/EFT bridge only: X = -1/2 g^mn d_m Phi d_n Phi, so Y = -2X.

"""
PHASE 40: Why the first charged framed branch is h = 2

Status:
    Candidate selection theorem for the h index used in phase38.

Problem:
    phase38 uses h=2 in

        theta = h / 9 = 2 / 9.

    Why not h=1?

Key distinction:
    h=1 is a valid closure only for a projective/nematic director, where

        n ~ -n.

    But the charged RFG sector is framed and oriented: it uses the oriented
    material triad and the epsilon_ABC volume form/topological current.
    Therefore the lift must close as an oriented frame, not merely as a
    projective line.

Minimal closure:
    Let h count half-turns of the local director/framing.

        h=0: trivial oriented closure.
        h=1: projective closure only; n -> -n, orientation not closed.
        h=2: first non-trivial oriented closure; n -> n.

    Thus h=2 is the first non-trivial charged/framed branch.

Consequence:
    Combining phase39 and this file:

        C3(axis) x C3(strain phase) -> 9 slots,
        h_charged = 2,
        theta = 2/9.
"""

import math

from phase38_z9_theta_holonomy import (
    closure_lattice_order,
    positivity_edge_for_c3_operator,
)


def director_after_half_turns(h):
    """
    Unit director in a 2D slice after h half-turns.

    h=1 gives n -> -n, which is projectively closed but not oriented closed.
    h=2 gives n -> n.
    """
    angle = math.pi * h
    return (round(math.cos(angle), 12), round(math.sin(angle), 12))


def closure_type(h):
    """Classify the h branch by projective and oriented closure."""
    n0 = (1.0, 0.0)
    nh = director_after_half_turns(h)
    oriented_closed = nh == n0
    projective_closed = nh == n0 or nh == (-n0[0], -n0[1])
    nontrivial = h != 0
    charged_allowed = oriented_closed and nontrivial
    return {
        "h": h,
        "director": nh,
        "projective_closed": projective_closed,
        "oriented_closed": oriented_closed,
        "nontrivial": nontrivial,
        "charged_framed_allowed": charged_allowed,
    }


def c3_frequency_raw(theta):
    return [
        1.0 + math.sqrt(2.0) * math.cos(theta + 2.0 * math.pi * k / 3.0)
        for k in range(3)
    ]


def branch_selection_table(max_h=6):
    """
    Combine projective/oriented closure with C3-spectrum positivity.
    """
    n_slots = closure_lattice_order()
    edge = positivity_edge_for_c3_operator()
    rows = []
    for h in range(max_h + 1):
        theta = h / n_slots
        raw = c3_frequency_raw(theta)
        positive = all(value > 0.0 for value in raw)
        closure = closure_type(h)
        selected = closure["charged_framed_allowed"] and positive and h == 2
        rows.append(
            {
                **closure,
                "theta": theta,
                "below_positive_edge": theta < edge,
                "positive_spectrum": positive,
                "selected": selected,
            }
        )
    return rows


def h2_selection_summary():
    n_slots = closure_lattice_order()
    h = 2
    theta = h / n_slots
    return {
        "closure_slots": n_slots,
        "h_selected": h,
        "reason_h0_rejected": "trivial branch, not a generation-splitting defect",
        "reason_h1_rejected": "projective/nematic closure only: n -> -n",
        "reason_h2_selected": "first non-trivial oriented framed closure: n -> n",
        "theta": theta,
        "theta_formula": "theta = h / 9 = 2 / 9",
        "remaining_open_point": "derive the oriented-frame requirement from the charged RFG coupling, not as a selection rule.",
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PHASE 40: Projective vs oriented closure selects h = 2")
    print("=" * 72)

    print("\n1. Branch selection table")
    for row in branch_selection_table():
        marker = " <-- selected" if row["selected"] else ""
        print(
            f"  h={row['h']} theta={row['theta']:.9f} "
            f"n={row['director']} "
            f"projective={row['projective_closed']} "
            f"oriented={row['oriented_closed']} "
            f"positive={row['positive_spectrum']}{marker}"
        )

    print("\n2. h=2 summary")
    for key, value in h2_selection_summary().items():
        if isinstance(value, float):
            print(f"  {key:24s}: {value:.12f}")
        else:
            print(f"  {key:24s}: {value}")
