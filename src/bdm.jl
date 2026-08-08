# Brezzi-Douglas-Marini element, degree 2
# (https://defelement.org/elements/brezzi-douglas-marini.html,
# DefElement: elements/brezzi-douglas-marini.def).
#
# Cells/degrees implemented: RefTriangle degree 2. Degree 1 exists in Ferrite
# as `Ferrite.BrezziDouglasMarini{RefTriangle, 1}` (hence the distinct type
# name here); degree 3+ and the tetrahedron would follow the same pattern
# (the tetrahedron beyond degree 1 additionally needs face-DOF orientation
# machinery that Ferrite does not have, cf. PLAN.md).
#
# Space: the full vector-valued P2 (12-dimensional). H(div) element with the
# contravariant Piola mapping. DOFs: three normal moments per edge (against
# the quadratic Lagrange basis of the edge, ordered along the FERRITE edge
# direction at the points s = 0, 1/2, 1 of the edge parametrization, with
# outward reference normal), plus three interior moments (symfem's, against
# the degree-1 Nedelec first-kind basis).
#
# Cross-cell orientation follows the in-tree pattern (RaviartThomas{RefTriangle,2},
# BrezziDouglasMarini{RefTriangle,1}): `adjust_dofs_during_distribution = true`
# reverses the intra-edge DOF order on edges whose global direction opposes
# the local one, and `get_direction` flips the sign of the edge shape
# functions (accounting for the opposite outward normal seen by the
# neighbor), so that shared DOFs represent the same global functional.
#
# Basis derivation: the action of the Ferrite-convention functionals above on
# symfem's "BDM" degree-2 basis is a signed permutation (verified exactly
# with sympy), so the basis below is the signed permutation of symfem's:
# Ferrite edge 1 = -(symfem edge 2 dofs, weight order s = 0, 1, 1/2),
# Ferrite edge 2 = +(symfem edge 1 reversed), Ferrite edge 3 = -(symfem edge
# 0), interior unchanged. Moment DOFs: not nodal, no `reference_coordinates`.

"""
    BDM{shape, order}()

Brezzi-Douglas-Marini H(div) element on the triangle, degree 2 (full vector
P2, contravariant Piola mapping): three normal-moment DOFs per edge and
three interior moments. Degree 1 lives in Ferrite itself as
`Ferrite.BrezziDouglasMarini`.
"""
struct BDM{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function BDM{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function BDM{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::BDM) = Ferrite.ContravariantPiolaMapping()
Ferrite.conformity(::BDM) = Ferrite.HdivConformity()

Ferrite.getnbasefunctions(::BDM{RefTriangle, 2}) = 12
Ferrite.edgedof_interior_indices(::BDM{RefTriangle, 2}) = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
Ferrite.facedof_interior_indices(::BDM{RefTriangle, 2}) = ((10, 11, 12),)
Ferrite.adjust_dofs_during_distribution(::BDM{RefTriangle, 2}) = true

function Ferrite.get_direction(::BDM{RefTriangle, 2}, shape_nr::Int, cell)
    shape_nr > 9 && return 1
    edge_nr = (shape_nr + 2) ÷ 3
    return Ferrite.get_edge_direction(cell, edge_nr)
end

function Ferrite.reference_shape_value(ip::BDM{RefTriangle, 2}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    # Edge 1: dofs at s = 0, 1/2, 1 along (1,0) -> (0,1)
    i == 1 && return Vec(18 * x^2 - 9 * x, -12 * x * y + 3 * y)
    i == 2 && return Vec((3 // 2) * x^2 + 9 * x * y - 3 * x, 9 * x * y + (3 // 2) * y^2 - 3 * y)
    i == 3 && return Vec(-12 * x * y + 3 * x, 18 * y^2 - 9 * y)
    # Edge 2: dofs at s = 0, 1/2, 1 along (0,1) -> (0,0)
    i == 4 && return Vec(-12 * x * y - 30 * y^2 + 3 * x + 24 * y - 3, 18 * y^2 - 9 * y)
    i == 5 && return Vec((-3 // 2) * x^2 + 15 * x * y + 15 * y^2 - 15 * y + 3 // 2, -9 * x * y + (-15 // 2) * y^2 + 6 * y)
    i == 6 && return Vec(-18 * x^2 - 48 * x * y - 30 * y^2 + 27 * x + 36 * y - 9, 12 * x * y + 12 * y^2 - 9 * y)
    # Edge 3: dofs at s = 0, 1/2, 1 along (0,0) -> (1,0)
    i == 7 && return Vec(12 * x^2 + 12 * x * y - 9 * x, -30 * x^2 - 48 * x * y - 18 * y^2 + 36 * x + 27 * y - 9)
    i == 8 && return Vec((-15 // 2) * x^2 - 9 * x * y + 6 * x, 15 * x^2 + 15 * x * y + (-3 // 2) * y^2 - 15 * x + 3 // 2)
    i == 9 && return Vec(18 * x^2 - 9 * x, -30 * x^2 - 12 * x * y + 24 * x + 3 * y - 3)
    # Interior moments
    i == 10 && return Vec(-36 * x^2 - 48 * x * y + 36 * x, 24 * x * y + 12 * y^2 - 12 * y)
    i == 11 && return Vec(12 * x^2 + 24 * x * y - 12 * x, -48 * x * y - 36 * y^2 + 36 * y)
    i == 12 && return Vec(-12 * x^2 - 48 * x * y + 12 * x, 48 * x * y + 12 * y^2 - 12 * y)
    return throw_out_of_range(ip, i)
end
