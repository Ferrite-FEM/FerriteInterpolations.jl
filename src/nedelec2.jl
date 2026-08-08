# Nedelec (second kind) element (https://defelement.org/elements/nedelec2.html,
# DefElement: elements/nedelec2.def).
#
# Cells/degrees implemented: RefTriangle degrees 1-2. The tetrahedron beyond
# degree 1 needs face-DOF orientation machinery Ferrite does not have
# (cf. PLAN.md); degree-1 tets and higher triangle degrees would follow the
# same pattern as below.
#
# Space: the full vector-valued P_k (unlike the first-kind family, which
# Ferrite ships as `Ferrite.Nedelec`). H(curl) element with the covariant
# Piola mapping. DOFs: k+1 tangential moments per edge (against the degree-k
# Lagrange basis of the edge, ordered along the FERRITE edge direction at the
# equispaced points of the edge parametrization, tangent = edge direction),
# plus (for degree 2) three interior moments (symfem's, against the degree-1
# Raviart-Thomas basis).
#
# Cross-cell orientation follows the in-tree Nedelec pattern:
# `adjust_dofs_during_distribution = true` reverses the intra-edge DOF order
# on edges whose global direction opposes the local one, and `get_direction`
# flips the sign (the neighbor traverses the edge, and hence sees the
# tangent, in the opposite direction).
#
# Basis derivation: as for src/bdm.jl, the action of the Ferrite-convention
# functionals on symfem's "N2curl" basis is exactly a signed permutation
# (verified with sympy), and the basis below is that signed permutation of
# symfem's. Moment DOFs: not nodal, no `reference_coordinates`.

"""
    NedelecSecondKind{shape, order}()

Nedelec element of the second kind on the triangle, degrees 1-2 (full vector
P_order, covariant Piola mapping): order+1 tangential-moment DOFs per edge
plus, for degree 2, three interior moments. H(curl)-conforming.
"""
struct NedelecSecondKind{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function NedelecSecondKind{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function NedelecSecondKind{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::NedelecSecondKind) = Ferrite.CovariantPiolaMapping()
Ferrite.conformity(::NedelecSecondKind) = Ferrite.HcurlConformity()
Ferrite.adjust_dofs_during_distribution(::NedelecSecondKind) = true

##########
# Degree 1
##########

Ferrite.getnbasefunctions(::NedelecSecondKind{RefTriangle, 1}) = 6
Ferrite.edgedof_interior_indices(::NedelecSecondKind{RefTriangle, 1}) = ((1, 2), (3, 4), (5, 6))
Ferrite.facedof_interior_indices(::NedelecSecondKind{RefTriangle, 1}) = ((),)

function Ferrite.get_direction(::NedelecSecondKind{RefTriangle, 1}, shape_nr::Int, cell)
    return Ferrite.get_edge_direction(cell, (shape_nr + 1) ÷ 2)
end

function Ferrite.reference_shape_value(ip::NedelecSecondKind{RefTriangle, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    # Edge 1: dofs at s = 0, 1 along (1,0) -> (0,1)
    i == 1 && return Vec(2 * y, 4 * x)
    i == 2 && return Vec(-4 * y, -2 * x)
    # Edge 2: dofs at s = 0, 1 along (0,1) -> (0,0)
    i == 3 && return Vec(-4 * y, -2 * x - 6 * y + 2)
    i == 4 && return Vec(2 * y, 4 * x + 6 * y - 4)
    # Edge 3: dofs at s = 0, 1 along (0,0) -> (1,0)
    i == 5 && return Vec(-6 * x - 4 * y + 4, -2 * x)
    i == 6 && return Vec(6 * x + 2 * y - 2, 4 * x)
    return throw_out_of_range(ip, i)
end

##########
# Degree 2
##########

Ferrite.getnbasefunctions(::NedelecSecondKind{RefTriangle, 2}) = 12
Ferrite.edgedof_interior_indices(::NedelecSecondKind{RefTriangle, 2}) = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
Ferrite.facedof_interior_indices(::NedelecSecondKind{RefTriangle, 2}) = ((10, 11, 12),)

function Ferrite.get_direction(::NedelecSecondKind{RefTriangle, 2}, shape_nr::Int, cell)
    shape_nr > 9 && return 1
    return Ferrite.get_edge_direction(cell, (shape_nr + 2) ÷ 3)
end

function Ferrite.reference_shape_value(ip::NedelecSecondKind{RefTriangle, 2}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    # Edge 1: dofs at s = 0, 1/2, 1 along (1,0) -> (0,1)
    i == 1 && return Vec(12 * x * y - 3 * y, 18 * x^2 - 9 * x)
    i == 2 && return Vec(-9 * x * y + (-3 // 2) * y^2 + 3 * y, (3 // 2) * x^2 + 9 * x * y - 3 * x)
    i == 3 && return Vec(-18 * y^2 + 9 * y, -12 * x * y + 3 * x)
    # Edge 2: dofs at s = 0, 1/2, 1 along (0,1) -> (0,0)
    i == 4 && return Vec(-18 * y^2 + 9 * y, -12 * x * y - 30 * y^2 + 3 * x + 24 * y - 3)
    i == 5 && return Vec(9 * x * y + (15 // 2) * y^2 - 6 * y, (-3 // 2) * x^2 + 15 * x * y + 15 * y^2 - 15 * y + 3 // 2)
    i == 6 && return Vec(-12 * x * y - 12 * y^2 + 9 * y, -18 * x^2 - 48 * x * y - 30 * y^2 + 27 * x + 36 * y - 9)
    # Edge 3: dofs at s = 0, 1/2, 1 along (0,0) -> (1,0)
    i == 7 && return Vec(30 * x^2 + 48 * x * y + 18 * y^2 - 36 * x - 27 * y + 9, 12 * x^2 + 12 * x * y - 9 * x)
    i == 8 && return Vec(-15 * x^2 - 15 * x * y + (3 // 2) * y^2 + 15 * x - 3 // 2, (-15 // 2) * x^2 - 9 * x * y + 6 * x)
    i == 9 && return Vec(30 * x^2 + 12 * x * y - 24 * x - 3 * y + 3, 18 * x^2 - 9 * x)
    # Interior moments
    i == 10 && return Vec(-24 * x * y - 12 * y^2 + 12 * y, -36 * x^2 - 48 * x * y + 36 * x)
    i == 11 && return Vec(48 * x * y + 36 * y^2 - 36 * y, 12 * x^2 + 24 * x * y - 12 * x)
    i == 12 && return Vec(-48 * x * y - 12 * y^2 + 12 * y, -12 * x^2 - 48 * x * y + 12 * x)
    return throw_out_of_range(ip, i)
end
