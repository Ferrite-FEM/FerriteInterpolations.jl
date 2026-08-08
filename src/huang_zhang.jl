# Huang-Zhang element (https://defelement.org/elements/huang-zhang.html,
# DefElement: elements/huang-zhang.def).
#
# Cells/degrees implemented: RefQuadrilateral degree 1 (higher degrees follow
# the same pattern with more edge and interior moments; not transcribed).
# H(div) element with the contravariant Piola mapping (Huang & Zhang, "A
# family of stable mixed finite elements for the linear elasticity problem
# on quadrilateral grids").
#
# The 12 DOFs all sit on the edges: per edge, two normal moments (linear
# Lagrange weights) and one zeroth tangential moment, arranged palindromically
# as (normal at s=0, tangential, normal at s=1) exactly as for
# src/mardal_tai_winther.jl (see there for why this arrangement makes
# Ferrite's order reversal + sign flip work on reversed edges).
#
# Basis derivation as for src/bdm.jl, with one extra step: the functional
# matrix is computed on symfem's [0, 1]^2 reference square (with Ferrite's
# edge order/directions) and the resulting dual basis is pushed forward to
# Ferrite's [-1, 1]^2 with the contravariant Piola map of x = (ξ + 1)/2,
# i.e. v(ξ) = v01((ξ + 1)/2)/2. Moment DOFs: not nodal, no
# reference_coordinates.
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# Same obstruction as src/mardal_tai_winther.jl (see the detailed notes
# there): the element's tangential-moment DOFs are not preserved by the
# contravariant Piola map (they become metric-weighted), so neighboring
# cells with different Jacobians couple different physical functionals.
# Verified numerically with the commented-out implementation below (which
# passes the reference-cell symfem cross-check and pointwise normal
# continuity): the tangential-moment jump vanishes only for congruent
# translate cells. Needs Kirby-style per-cell DOF transformations upstream.

"""
    HuangZhang{shape, order}

Placeholder for the Huang-Zhang element. Not usable: the element needs
cell-dependent DOF transformations that Ferrite does not support, see the
`STATUS: BLOCKED` notes in `src/huang_zhang.jl` and
`src/mardal_tai_winther.jl`.
"""
struct HuangZhang{shape, order}
    function HuangZhang{shape, order}() where {shape, order}
        return error(
            "HuangZhang is not implemented: its tangential-moment DOFs need " *
                "cell-dependent (Kirby-style) DOF transformations that Ferrite does not " *
                "support; see src/huang_zhang.jl."
        )
    end
end

#=
"""
    HuangZhang{shape, order}()

Huang-Zhang H(div) element on the quadrilateral, degree 1: two normal
moments and one tangential moment per edge (12 DOFs, all shared with the
neighboring cells).
"""
struct HuangZhang{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function HuangZhang{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function HuangZhang{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::HuangZhang) = Ferrite.ContravariantPiolaMapping()
Ferrite.conformity(::HuangZhang) = Ferrite.HdivConformity()
Ferrite.adjust_dofs_during_distribution(::HuangZhang) = true

Ferrite.getnbasefunctions(::HuangZhang{RefQuadrilateral, 1}) = 12
Ferrite.edgedof_interior_indices(::HuangZhang{RefQuadrilateral, 1}) = ((1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12))
Ferrite.facedof_interior_indices(::HuangZhang{RefQuadrilateral, 1}) = ((),)

function Ferrite.get_direction(::HuangZhang{RefQuadrilateral, 1}, shape_nr::Int, cell)
    return Ferrite.get_edge_direction(cell, (shape_nr + 2) ÷ 3)
end

function Ferrite.reference_shape_value(ip::HuangZhang{RefQuadrilateral, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    # Edge 1: (normal at s=0, tangential, normal at s=1) along (-1,-1) -> (1,-1)
    i == 1 && return Vec(zero(x), (9 // 8) * x * y^2 + (-3 // 4) * x * y + (-3 // 8) * y^2 + (-3 // 8) * x + (1 // 4) * y + 1 // 8)
    i == 2 && return Vec((3 // 8) * x^2 * y + (-3 // 8) * x^2 + (-3 // 8) * y + 3 // 8, zero(x))
    i == 3 && return Vec(zero(x), (-9 // 8) * x * y^2 + (3 // 4) * x * y + (-3 // 8) * y^2 + (3 // 8) * x + (1 // 4) * y + 1 // 8)
    # Edge 2: along (1,-1) -> (1,1)
    i == 4 && return Vec((-9 // 8) * x^2 * y + (3 // 8) * x^2 + (-3 // 4) * x * y + (1 // 4) * x + (3 // 8) * y - 1 // 8, zero(x))
    i == 5 && return Vec(zero(x), (-3 // 8) * x * y^2 + (-3 // 8) * y^2 + (3 // 8) * x + 3 // 8)
    i == 6 && return Vec((9 // 8) * x^2 * y + (3 // 8) * x^2 + (3 // 4) * x * y + (1 // 4) * x + (-3 // 8) * y - 1 // 8, zero(x))
    # Edge 3: along (1,1) -> (-1,1)
    i == 7 && return Vec(zero(x), (9 // 8) * x * y^2 + (3 // 4) * x * y + (3 // 8) * y^2 + (-3 // 8) * x + (1 // 4) * y - 1 // 8)
    i == 8 && return Vec((3 // 8) * x^2 * y + (3 // 8) * x^2 + (-3 // 8) * y - 3 // 8, zero(x))
    i == 9 && return Vec(zero(x), (-9 // 8) * x * y^2 + (-3 // 4) * x * y + (3 // 8) * y^2 + (3 // 8) * x + (1 // 4) * y - 1 // 8)
    # Edge 4: along (-1,1) -> (-1,-1)
    i == 10 && return Vec((-9 // 8) * x^2 * y + (-3 // 8) * x^2 + (3 // 4) * x * y + (1 // 4) * x + (3 // 8) * y + 1 // 8, zero(x))
    i == 11 && return Vec(zero(x), (-3 // 8) * x * y^2 + (3 // 8) * y^2 + (3 // 8) * x - 3 // 8)
    i == 12 && return Vec((9 // 8) * x^2 * y + (-3 // 8) * x^2 + (-3 // 4) * x * y + (1 // 4) * x + (-3 // 8) * y + 1 // 8, zero(x))
    return throw_out_of_range(ip, i)
end
=#
