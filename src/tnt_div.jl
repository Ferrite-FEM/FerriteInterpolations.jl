# TNTDiv element (https://defelement.org/elements/tnt-div.html,
# DefElement: elements/tnt-div.def).
#
# Cells/degrees implemented: RefQuadrilateral degree 1 (higher degrees and
# the hexahedron follow the same pattern; the hexahedron additionally needs
# face-DOF orientation machinery Ferrite does not have, cf. PLAN.md).
# ContravariantPiola mapping. DOFs: two normal moments per edge (linear Lagrange weights, ordered
# along the FERRITE edge direction, outward reference normal),
# plus 3 interior moments (symfem's).
#
# Orientation handling and basis derivation exactly as for src/bdm.jl /
# src/nedelec2.jl: the action of the Ferrite-convention functionals on
# symfem's "TNTdiv" basis is a signed permutation (verified with sympy);
# the dual basis is computed on symfem's [0, 1]^2 cell and pushed forward to
# Ferrite's [-1, 1]^2 with the Piola factor 1/2 (contravariant) of x = (xi + 1)/2.
# Moment DOFs: not nodal, no reference_coordinates.

"""
    TNTDiv{shape, order}()

Tiniest-tensor H(div) element on the quadrilateral, degree 1: two normal
moments per edge plus three interior moments.
"""
struct TNTDiv{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function TNTDiv{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function TNTDiv{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::TNTDiv) = Ferrite.ContravariantPiolaMapping()
Ferrite.conformity(::TNTDiv) = Ferrite.HdivConformity()
Ferrite.adjust_dofs_during_distribution(::TNTDiv) = true

Ferrite.getnbasefunctions(::TNTDiv{RefQuadrilateral, 1}) = 11
Ferrite.edgedof_interior_indices(::TNTDiv{RefQuadrilateral, 1}) = ((1, 2), (3, 4), (5, 6), (7, 8))
Ferrite.facedof_interior_indices(::TNTDiv{RefQuadrilateral, 1}) = ((9, 10, 11),)

function Ferrite.get_direction(::TNTDiv{RefQuadrilateral, 1}, shape_nr::Int, cell)
    shape_nr > 8 && return 1
    return Ferrite.get_edge_direction(cell, (shape_nr + 1) ÷ 2)
end

function Ferrite.reference_shape_value(ip::TNTDiv{RefQuadrilateral, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return Vec((9 // 16) * x^2 * y + (-9 // 16) * y, (9 // 16) * x * y^2 + (-3 // 4) * x * y + (-3 // 8) * y^2 + (3 // 16) * x + (1 // 4) * y + 1 // 8)
    i == 2 && return Vec((-9 // 16) * x^2 * y + (9 // 16) * y, (-9 // 16) * x * y^2 + (3 // 4) * x * y + (-3 // 8) * y^2 + (-3 // 16) * x + (1 // 4) * y + 1 // 8)
    i == 3 && return Vec((-9 // 16) * x^2 * y + (3 // 8) * x^2 + (-3 // 4) * x * y + (1 // 4) * x + (-3 // 16) * y - 1 // 8, (-9 // 16) * x * y^2 + (9 // 16) * x)
    i == 4 && return Vec((9 // 16) * x^2 * y + (3 // 8) * x^2 + (3 // 4) * x * y + (1 // 4) * x + (3 // 16) * y - 1 // 8, (9 // 16) * x * y^2 + (-9 // 16) * x)
    i == 5 && return Vec((9 // 16) * x^2 * y + (-9 // 16) * y, (9 // 16) * x * y^2 + (3 // 4) * x * y + (3 // 8) * y^2 + (3 // 16) * x + (1 // 4) * y - 1 // 8)
    i == 6 && return Vec((-9 // 16) * x^2 * y + (9 // 16) * y, (-9 // 16) * x * y^2 + (-3 // 4) * x * y + (3 // 8) * y^2 + (-3 // 16) * x + (1 // 4) * y - 1 // 8)
    i == 7 && return Vec((-9 // 16) * x^2 * y + (-3 // 8) * x^2 + (3 // 4) * x * y + (1 // 4) * x + (-3 // 16) * y + 1 // 8, (-9 // 16) * x * y^2 + (9 // 16) * x)
    i == 8 && return Vec((9 // 16) * x^2 * y + (-3 // 8) * x^2 + (-3 // 4) * x * y + (1 // 4) * x + (3 // 16) * y + 1 // 8, (9 // 16) * x * y^2 + (-9 // 16) * x)
    i == 9 && return Vec((9 // 8) * x^2 * y + (-9 // 8) * y, (9 // 8) * x * y^2 + (-3 // 4) * y^2 + (-9 // 8) * x + 3 // 4)
    i == 10 && return Vec((9 // 8) * x^2 * y + (-3 // 4) * x^2 + (-9 // 8) * y + 3 // 4, (9 // 8) * x * y^2 + (-9 // 8) * x)
    i == 11 && return Vec((-9 // 4) * x^2 * y + (9 // 4) * y, (-9 // 4) * x * y^2 + (9 // 4) * x)
    return throw_out_of_range(ip, i)
end
