# TrimmedSerendipityCurl element (https://defelement.org/elements/trimmed-serendipity-curl.html,
# DefElement: elements/trimmed-serendipity-curl.def).
#
# Cells/degrees implemented: RefQuadrilateral degree 1 (higher degrees and
# the hexahedron follow the same pattern; the hexahedron additionally needs
# face-DOF orientation machinery Ferrite does not have, cf. PLAN.md).
# CovariantPiola mapping. DOFs: two tangential moments per edge (linear Lagrange weights, ordered
# along the FERRITE edge direction, tangent = edge direction),
# plus 2 interior moments (symfem's).
#
# Orientation handling and basis derivation exactly as for src/bdm.jl /
# src/nedelec2.jl: the action of the Ferrite-convention functionals on
# symfem's "TScurl" basis is a signed permutation (verified with sympy);
# the dual basis is computed on symfem's [0, 1]^2 cell and pushed forward to
# Ferrite's [-1, 1]^2 with the Piola factor 2 (covariant) of x = (xi + 1)/2.
# Moment DOFs: not nodal, no reference_coordinates.

"""
    TrimmedSerendipityCurl{shape, order}()

Trimmed serendipity H(curl) element on the quadrilateral, degree 1 (S-minus
family): two tangential moments per edge plus two interior moments.
"""
struct TrimmedSerendipityCurl{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function TrimmedSerendipityCurl{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function TrimmedSerendipityCurl{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::TrimmedSerendipityCurl) = Ferrite.CovariantPiolaMapping()
Ferrite.conformity(::TrimmedSerendipityCurl) = Ferrite.HcurlConformity()
Ferrite.adjust_dofs_during_distribution(::TrimmedSerendipityCurl) = true

Ferrite.getnbasefunctions(::TrimmedSerendipityCurl{RefQuadrilateral, 1}) = 10
Ferrite.edgedof_interior_indices(::TrimmedSerendipityCurl{RefQuadrilateral, 1}) = ((1, 2), (3, 4), (5, 6), (7, 8))
Ferrite.facedof_interior_indices(::TrimmedSerendipityCurl{RefQuadrilateral, 1}) = ((9, 10),)

function Ferrite.get_direction(::TrimmedSerendipityCurl{RefQuadrilateral, 1}, shape_nr::Int, cell)
    shape_nr > 8 && return 1
    return Ferrite.get_edge_direction(cell, (shape_nr + 1) ÷ 2)
end

function Ferrite.reference_shape_value(ip::TrimmedSerendipityCurl{RefQuadrilateral, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return Vec(3 * x * y + (3 // 2) * y^2 - 3 * x - y - 1 // 2, zero(x))
    i == 2 && return Vec(-3 * x * y + (3 // 2) * y^2 + 3 * x - y - 1 // 2, zero(x))
    i == 3 && return Vec(zero(x), (3 // 2) * x^2 - 3 * x * y + x - 3 * y - 1 // 2)
    i == 4 && return Vec(zero(x), (3 // 2) * x^2 + 3 * x * y + x + 3 * y - 1 // 2)
    i == 5 && return Vec(-3 * x * y + (-3 // 2) * y^2 - 3 * x - y + 1 // 2, zero(x))
    i == 6 && return Vec(3 * x * y + (-3 // 2) * y^2 + 3 * x - y + 1 // 2, zero(x))
    i == 7 && return Vec(zero(x), (-3 // 2) * x^2 + 3 * x * y + x - 3 * y + 1 // 2)
    i == 8 && return Vec(zero(x), (-3 // 2) * x^2 - 3 * x * y + x + 3 * y + 1 // 2)
    i == 9 && return Vec(zero(x), 3 * x^2 - 3)
    i == 10 && return Vec(-3 * y^2 + 3, zero(x))
    return throw_out_of_range(ip, i)
end
