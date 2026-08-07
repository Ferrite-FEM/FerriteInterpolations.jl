"""Generate literal Julia reference tables for test/test_bernstein.jl.

Tabulates symfem's Bernstein basis at a few points and prints the values in
Ferrite's DOF order, with the point coordinates in Ferrite's reference cells.

Coordinate maps (Ferrite -> symfem):
- interval: Ferrite RefLine is [-1, 1], symfem [0, 1]: x_s = (x_f + 1) / 2
- triangle/tetrahedron: identical coordinates, but different vertex numbering
  (Ferrite RefTriangle v1=(1,0), v2=(0,1), v3=(0,0); symfem v0=(0,0),
  v1=(1,0), v2=(0,1); tetrahedron: Ferrite v_i = symfem v_{i-1}).

DOF order mapping is built structurally from symfem's entity_dofs:
Ferrite order is vertices -> edge interiors (edge DOFs from first towards
second vertex of the Ferrite reference edge) -> cell interior. symfem edge
DOFs are assumed ordered along the symfem edge from its first to its second
vertex; the list is reversed when the Ferrite edge maps onto the symfem edge
with opposite orientation.
"""

from fractions import Fraction

import symfem

# Ferrite reference edges (1-based vertex numbers), from Ferrite src/Grid/grid.jl.
FERRITE_EDGES = {
    "interval": [(1, 2)],
    "triangle": [(1, 2), (2, 3), (3, 1)],
    "tetrahedron": [(1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)],
}

# Ferrite vertex coordinates, from Ferrite src/docs.jl.
FERRITE_VERTICES = {
    "interval": [(-1,), (1,)],
    "triangle": [(1, 0), (0, 1), (0, 0)],
    "tetrahedron": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
}

# Test points in Ferrite reference coordinates (all strictly inside).
FERRITE_POINTS = {
    "interval": [(Fraction(-4, 5),), (Fraction(-1, 5),), (Fraction(1, 2),), (Fraction(9, 10),)],
    "triangle": [
        (Fraction(3, 10), Fraction(1, 5)),
        (Fraction(1, 10), Fraction(3, 5)),
        (Fraction(1, 4), Fraction(1, 4)),
        (Fraction(1, 20), Fraction(9, 10)),
    ],
    "tetrahedron": [
        (Fraction(1, 5), Fraction(3, 10), Fraction(2, 5)),
        (Fraction(1, 10), Fraction(1, 10), Fraction(1, 10)),
        (Fraction(1, 2), Fraction(1, 5), Fraction(1, 4)),
        (Fraction(1, 20), Fraction(3, 5), Fraction(3, 10)),
    ],
}


def ferrite_to_symfem_point(cell, pt):
    if cell == "interval":
        return tuple((x + 1) / 2 for x in pt)
    return pt


def vertex_map(cell, el):
    """Ferrite vertex number (1-based) -> symfem vertex index (0-based)."""
    svertices = [tuple(Fraction(c) for c in v) for v in el.reference.vertices]
    fmap = {}
    for fv, coord in enumerate(FERRITE_VERTICES[cell], start=1):
        scoord = ferrite_to_symfem_point(cell, tuple(Fraction(c) for c in coord))
        fmap[fv] = svertices.index(scoord)
    return fmap


def ferrite_dof_order(cell, el):
    """List of symfem dof indices in Ferrite DOF order."""
    vmap = vertex_map(cell, el)
    tdim = el.reference.tdim
    order = []
    # Vertices
    for fv in range(1, len(FERRITE_VERTICES[cell]) + 1):
        dofs = el.entity_dofs(0, vmap[fv])
        assert len(dofs) == 1, (cell, fv, dofs)
        order.extend(dofs)
    # Edge interiors (for the interval the single "edge" is the cell itself)
    sedges = [tuple(e) for e in el.reference.edges]
    for fa, fb in FERRITE_EDGES[cell]:
        sa, sb = vmap[fa], vmap[fb]
        if cell == "interval":
            dofs = el.entity_dofs(1, 0)
            reverse = (sa, sb) != (0, 1)
        else:
            matches = [i for i, e in enumerate(sedges) if set(e) == {sa, sb}]
            assert len(matches) == 1, (cell, (fa, fb), matches)
            (eidx,) = matches
            dofs = el.entity_dofs(1, eidx)
            reverse = sedges[eidx] != (sa, sb)
        order.extend(reversed(dofs) if reverse else dofs)
    # Cell interior (2D face / 3D volume)
    if tdim >= 2:
        for d, n in [(2, f) for f in range(len(el.reference.faces))] if tdim == 2 else []:
            order.extend(el.entity_dofs(d, n))
        if tdim == 2:
            pass
        else:
            for f in range(len(el.reference.faces)):
                order.extend(el.entity_dofs(2, f))
            order.extend(el.entity_dofs(3, 0))
    assert sorted(order) == list(range(el.space_dim)), (cell, order)
    return order


def main():
    for cell, degrees in [("interval", [1, 2, 3]), ("triangle", [1, 2, 3]), ("tetrahedron", [1, 2])]:
        for degree in degrees:
            el = symfem.create_element(cell, "Bernstein", degree)
            order = ferrite_dof_order(cell, el)
            basis = el.get_basis_functions()
            print(f"# {cell} degree {degree}: symfem dofs in Ferrite order: {order}")
            for fpt in FERRITE_POINTS[cell]:
                spt = ferrite_to_symfem_point(cell, fpt)
                vals = [basis[i].subs(symfem.symbols.x, spt) for i in order]
                fcoords = ", ".join(str(float(c)) for c in fpt)
                jvals = ", ".join(repr(float(v.as_sympy())) for v in vals)
                print(f"    (Vec(({fcoords},)), [{jvals}]),")
            print()


if __name__ == "__main__":
    main()
