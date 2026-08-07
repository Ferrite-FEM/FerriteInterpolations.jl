"""Generate the literal Julia reference table for test/test_fortin_soulie.jl.

Tabulates symfem's Fortin-Soulie basis at a few points and prints the values
in Ferrite DOF order. The triangle coordinates agree pointwise between symfem
and Ferrite; only entity numbering differs.

Unlike generate_bernstein_table.py, the DOF matching is geometric: symfem's
`entity_dofs` edge numbering does not consistently follow `reference.edges`
across element families, but all Fortin-Soulie DOFs are point evaluations, so
each DOF is assigned to its Ferrite entity by locating its point (edge DOFs
ordered by the parameter along the Ferrite edge, from first to second vertex).
"""

from fractions import Fraction

import symfem

FERRITE_EDGES = [(1, 2), (2, 3), (3, 1)]
FERRITE_VERTICES = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)), (Fraction(0), Fraction(0))]

POINTS = [
    (Fraction(3, 10), Fraction(1, 5)),
    (Fraction(1, 10), Fraction(3, 5)),
    (Fraction(1, 4), Fraction(1, 4)),
    (Fraction(1, 20), Fraction(9, 10)),
]


def edge_parameter(a, b, p):
    """Return t with p == a + t * (b - a) if p is on segment [a, b], else None."""
    d = tuple(bi - ai for ai, bi in zip(a, b))
    ts = set()
    for pi, ai, di in zip(p, a, d):
        if di == 0:
            if pi != ai:
                return None
        else:
            ts.add(Fraction(pi - ai, 1) / di)
    if len(ts) != 1:
        return None
    (t,) = ts
    return t if 0 < t < 1 else None


def main():
    el = symfem.create_element("triangle", "Fortin-Soulie", 2)
    points = [tuple(Fraction(str(c)) for c in dof.point) for dof in el.dofs]
    order = []
    for fa, fb in FERRITE_EDGES:
        a, b = FERRITE_VERTICES[fa - 1], FERRITE_VERTICES[fb - 1]
        on_edge = [(t, i) for i, p in enumerate(points)
                   if (t := edge_parameter(a, b, p)) is not None]
        order.extend(i for _, i in sorted(on_edge))
    interior = [i for i in range(el.space_dim) if i not in order]
    order.extend(interior)
    assert sorted(order) == list(range(el.space_dim)), order

    basis = el.get_basis_functions()
    print(f"# symfem dofs in Ferrite order: {order}")
    print("# Ferrite reference coordinates of the DOFs:")
    for i in order:
        print(f"#   {points[i]}")
    for pt in POINTS:
        vals = [basis[i].subs(symfem.symbols.x, pt) for i in order]
        fcoords = ", ".join(str(float(c)) for c in pt)
        jvals = ", ".join(repr(float(v.as_sympy())) for v in vals)
        print(f"    (Vec(({fcoords})), [{jvals}]),")


if __name__ == "__main__":
    main()
