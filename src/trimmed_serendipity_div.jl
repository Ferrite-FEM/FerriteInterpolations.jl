# TrimmedSerendipityDiv element (https://defelement.org/elements/trimmed-serendipity-div.html,
# DefElement: elements/trimmed-serendipity-div.def).
#
# Cells/degrees implemented: RefQuadrilateral degree 1 (higher degrees and
# the hexahedron follow the same pattern; the hexahedron additionally needs
# face-DOF orientation machinery Ferrite does not have, cf. PLAN.md).
# ContravariantPiola mapping. DOFs: two normal moments per edge (linear Lagrange weights, ordered
# along the FERRITE edge direction, outward reference normal),
# plus 2 interior moments (symfem's).
#
# Orientation handling and basis derivation exactly as for src/bdm.jl /
# src/nedelec2.jl: the action of the Ferrite-convention functionals on
# symfem's "TSdiv" basis is a signed permutation (verified with sympy);
# the dual basis is computed on symfem's [0, 1]^2 cell and pushed forward to
# Ferrite's [-1, 1]^2 with the Piola factor 1/2 (contravariant) of x = (xi + 1)/2.
# Moment DOFs: not nodal, no reference_coordinates.

"""
    TrimmedSerendipityDiv{shape, order}()

Trimmed serendipity H(div) element on the quadrilateral, degree 1 (S-minus
family): two normal moments per edge plus two interior moments.
"""
struct TrimmedSerendipityDiv{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function TrimmedSerendipityDiv{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function TrimmedSerendipityDiv{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::TrimmedSerendipityDiv) = Ferrite.ContravariantPiolaMapping()
Ferrite.conformity(::TrimmedSerendipityDiv) = Ferrite.HdivConformity()
Ferrite.adjust_dofs_during_distribution(::TrimmedSerendipityDiv) = true

Ferrite.getnbasefunctions(::TrimmedSerendipityDiv{RefQuadrilateral, 1}) = 10
Ferrite.edgedof_interior_indices(::TrimmedSerendipityDiv{RefQuadrilateral, 1}) = ((1, 2), (3, 4), (5, 6), (7, 8))
Ferrite.facedof_interior_indices(::TrimmedSerendipityDiv{RefQuadrilateral, 1}) = ((9, 10),)

function Ferrite.get_direction(::TrimmedSerendipityDiv{RefQuadrilateral, 1}, shape_nr::Int, cell)
    shape_nr > 8 && return 1
    return Ferrite.get_edge_direction(cell, (shape_nr + 1) ÷ 2)
end

function Ferrite.reference_shape_value(ip::TrimmedSerendipityDiv{RefQuadrilateral, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return Vec(zero(x), (-3 // 4) * x * y + (-3 // 8) * y^2 + (3 // 4) * x + (1 // 4) * y + 1 // 8)
    i == 2 && return Vec(zero(x), (3 // 4) * x * y + (-3 // 8) * y^2 + (-3 // 4) * x + (1 // 4) * y + 1 // 8)
    i == 3 && return Vec((3 // 8) * x^2 + (-3 // 4) * x * y + (1 // 4) * x + (-3 // 4) * y - 1 // 8, zero(x))
    i == 4 && return Vec((3 // 8) * x^2 + (3 // 4) * x * y + (1 // 4) * x + (3 // 4) * y - 1 // 8, zero(x))
    i == 5 && return Vec(zero(x), (3 // 4) * x * y + (3 // 8) * y^2 + (3 // 4) * x + (1 // 4) * y - 1 // 8)
    i == 6 && return Vec(zero(x), (-3 // 4) * x * y + (3 // 8) * y^2 + (-3 // 4) * x + (1 // 4) * y - 1 // 8)
    i == 7 && return Vec((-3 // 8) * x^2 + (3 // 4) * x * y + (1 // 4) * x + (-3 // 4) * y + 1 // 8, zero(x))
    i == 8 && return Vec((-3 // 8) * x^2 + (-3 // 4) * x * y + (1 // 4) * x + (3 // 4) * y + 1 // 8, zero(x))
    i == 9 && return Vec((-3 // 4) * x^2 + 3 // 4, zero(x))
    i == 10 && return Vec(zero(x), (-3 // 4) * y^2 + 3 // 4)
    return throw_out_of_range(ip, i)
end
