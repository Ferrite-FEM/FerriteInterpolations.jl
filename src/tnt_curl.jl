# TNTCurl element (https://defelement.org/elements/tnt-curl.html,
# DefElement: elements/tnt-curl.def).
#
# Cells/degrees implemented: RefQuadrilateral degree 1 (higher degrees and
# the hexahedron follow the same pattern; the hexahedron additionally needs
# face-DOF orientation machinery Ferrite does not have, cf. PLAN.md).
# CovariantPiola mapping. DOFs: two tangential moments per edge (linear Lagrange weights, ordered
# along the FERRITE edge direction, tangent = edge direction),
# plus 3 interior moments (symfem's).
#
# Orientation handling and basis derivation exactly as for src/bdm.jl /
# src/nedelec2.jl: the action of the Ferrite-convention functionals on
# symfem's "TNTcurl" basis is a signed permutation (verified with sympy);
# the dual basis is computed on symfem's [0, 1]^2 cell and pushed forward to
# Ferrite's [-1, 1]^2 with the Piola factor 2 (covariant) of x = (xi + 1)/2.
# Moment DOFs: not nodal, no reference_coordinates.

"""
    TNTCurl{shape, order}()

Tiniest-tensor H(curl) element on the quadrilateral, degree 1: two
tangential moments per edge plus three interior moments.
"""
struct TNTCurl{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function TNTCurl{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function TNTCurl{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::TNTCurl) = Ferrite.CovariantPiolaMapping()
Ferrite.conformity(::TNTCurl) = Ferrite.HcurlConformity()
Ferrite.adjust_dofs_during_distribution(::TNTCurl) = true

Ferrite.getnbasefunctions(::TNTCurl{RefQuadrilateral, 1}) = 11
Ferrite.edgedof_interior_indices(::TNTCurl{RefQuadrilateral, 1}) = ((1, 2), (3, 4), (5, 6), (7, 8))
Ferrite.facedof_interior_indices(::TNTCurl{RefQuadrilateral, 1}) = ((9, 10, 11),)

function Ferrite.get_direction(::TNTCurl{RefQuadrilateral, 1}, shape_nr::Int, cell)
    shape_nr > 8 && return 1
    return Ferrite.get_edge_direction(cell, (shape_nr + 1) ÷ 2)
end

function Ferrite.reference_shape_value(ip::TNTCurl{RefQuadrilateral, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return Vec((-9 // 4) * x * y^2 + 3 * x * y + (3 // 2) * y^2 + (-3 // 4) * x - y - 1 // 2, (9 // 4) * x^2 * y + (-9 // 4) * y)
    i == 2 && return Vec((9 // 4) * x * y^2 - 3 * x * y + (3 // 2) * y^2 + (3 // 4) * x - y - 1 // 2, (-9 // 4) * x^2 * y + (9 // 4) * y)
    i == 3 && return Vec((9 // 4) * x * y^2 + (-9 // 4) * x, (-9 // 4) * x^2 * y + (3 // 2) * x^2 - 3 * x * y + x + (-3 // 4) * y - 1 // 2)
    i == 4 && return Vec((-9 // 4) * x * y^2 + (9 // 4) * x, (9 // 4) * x^2 * y + (3 // 2) * x^2 + 3 * x * y + x + (3 // 4) * y - 1 // 2)
    i == 5 && return Vec((-9 // 4) * x * y^2 - 3 * x * y + (-3 // 2) * y^2 + (-3 // 4) * x - y + 1 // 2, (9 // 4) * x^2 * y + (-9 // 4) * y)
    i == 6 && return Vec((9 // 4) * x * y^2 + 3 * x * y + (-3 // 2) * y^2 + (3 // 4) * x - y + 1 // 2, (-9 // 4) * x^2 * y + (9 // 4) * y)
    i == 7 && return Vec((9 // 4) * x * y^2 + (-9 // 4) * x, (-9 // 4) * x^2 * y + (-3 // 2) * x^2 + 3 * x * y + x + (-3 // 4) * y + 1 // 2)
    i == 8 && return Vec((-9 // 4) * x * y^2 + (9 // 4) * x, (9 // 4) * x^2 * y + (-3 // 2) * x^2 - 3 * x * y + x + (3 // 4) * y + 1 // 2)
    i == 9 && return Vec((9 // 2) * x * y^2 - 3 * y^2 + (-9 // 2) * x + 3, (-9 // 2) * x^2 * y + (9 // 2) * y)
    i == 10 && return Vec((9 // 2) * x * y^2 + (-9 // 2) * x, (-9 // 2) * x^2 * y + 3 * x^2 + (9 // 2) * y - 3)
    i == 11 && return Vec(-9 * x * y^2 + 9 * x, 9 * x^2 * y - 9 * y)
    return throw_out_of_range(ip, i)
end
