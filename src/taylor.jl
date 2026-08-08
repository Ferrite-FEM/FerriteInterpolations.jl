# Taylor element (https://defelement.org/elements/taylor.html,
# DefElement: elements/taylor.def).
#
# Cells/degrees implemented: RefLine and RefTriangle, degrees 1-2 (higher
# degrees follow the same pattern with higher-derivative DOF blocks; the
# tetrahedron likewise). All DOFs are cell DOFs: one integral moment over the
# cell plus point evaluations of all derivatives up to the polynomial degree
# at the cell midpoint. The element is a discontinuous ("DG Taylor") basis;
# nothing is shared between cells, so L2Conformity is used (DefElement's
# metadata says `sobolev: H1`, which cannot hold for an element without any
# boundary DOFs; this file follows the actual structure).
#
# Mapping: DefElement refers to Kirby's transformation theory ("A general
# approach to transforming finite elements") -- the physical DOFs (cell
# integral, physical derivatives at the midpoint) are related to the
# reference DOFs by a CELL-DEPENDENT block-diagonal matrix M(J):
#
#     integral:     d(v) = det(J) d̂(v̂)               (block  det J)
#     gradient:     ∇ₓv(mid) = J⁻ᵀ ∇̂v̂(midref)        (block  J⁻ᵀ)
#     hessian:      Hₓ = J⁻ᵀ Ĥ J⁻¹                    (block  S(J⁻¹))
#
# Because every DOF is cell-local, this transformation never has to be
# coordinated across cells, so it CAN be expressed inside Ferrite's mapping
# machinery: `TaylorMapping` below computes B = M(J)⁻ᵀ in `apply_mapping!`
# and maps the basis as N_phys = B N̂ (and gradients accordingly). Only valid
# for AFFINE geometries (constant J; exact for straight-edged simplex cells),
# which is also the only geometry the element is defined on here.
#
# The reference basis is the dual of the Ferrite-reference DOFs (derived with
# sympy); it relates to symfem's Taylor basis by a scaled permutation (the
# scaling is the reference-map derivative factor on the interval; symfem
# orders the derivative DOFs lexicographically in the multi-index, e.g.
# (0,1), (0,2), (1,0), (1,1), (2,0) for the degree-2 triangle).
#
# Not nodal (integral + derivative DOFs): no `reference_coordinates`.

"""
    TaylorMapping()

Cell-dependent (Kirby-style) mapping for the Taylor element: the reference
basis is recombined per cell with the block-diagonal matrix M(J)⁻ᵀ so that
the physical DOFs (cell integral and midpoint derivatives) remain dual to
the mapped basis. Valid for affine geometries only.
"""
struct TaylorMapping end

"""
    Taylor{shape, order}()

Taylor element on RefLine/RefTriangle, degrees 1-2: a discontinuous basis
with one cell-integral DOF and point evaluations of the derivatives up to
`order` at the cell midpoint. Uses [`TaylorMapping`](@ref).
"""
struct Taylor{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::Taylor) = Ferrite.L2Conformity()
Ferrite.adjust_dofs_during_distribution(::Taylor) = false
Ferrite.mapping_type(::Taylor) = TaylorMapping()

# All DOFs in the cell interior, like DiscontinuousLagrange.
Ferrite.volumedof_interior_indices(ip::Taylor) = ntuple(i -> i, Ferrite.getnbasefunctions(ip))

Ferrite.getnbasefunctions(::Taylor{RefLine, 1}) = 2
Ferrite.getnbasefunctions(::Taylor{RefLine, 2}) = 3
Ferrite.getnbasefunctions(::Taylor{RefTriangle, 1}) = 3
Ferrite.getnbasefunctions(::Taylor{RefTriangle, 2}) = 6

function Ferrite.reference_shape_value(ip::Taylor{RefLine, order}, ξ::Vec{1}, i::Int) where {order}
    x = ξ[1]
    i == 1 && return 1 // 2 + zero(x)
    i == 2 && return x
    order >= 2 && i == 3 && return x^2 / 2 - 1 // 6
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_shape_value(ip::Taylor{RefTriangle, order}, ξ::Vec{2}, i::Int) where {order}
    x, y = ξ[1], ξ[2]
    i == 1 && return 2 + zero(x)
    i == 2 && return x - 1 // 3
    i == 3 && return y - 1 // 3
    if order >= 2
        i == 4 && return x^2 / 2 - x / 3 + 1 // 36
        i == 5 && return x * y - x / 3 - y / 3 + 5 // 36
        i == 6 && return y^2 / 2 - y / 3 + 1 // 36
    end
    return throw_out_of_range(ip, i)
end

# ---- the mapping -------------------------------------------------------------

Ferrite.required_geo_diff_order(::TaylorMapping, fun_diff_order::Int) = max(1, fun_diff_order)

# B = M(J)⁻ᵀ as a dense (N x N) tuple-of-rows for the given interpolation.
@inline function _taylor_B(::Taylor{RefLine, order}, J::Tensors.Tensor{2, 1}) where {order}
    j = J[1, 1]
    order == 1 && return ((1 / j, 0.0), (0.0, j))
    return ((1 / j, 0.0, 0.0), (0.0, j, 0.0), (0.0, 0.0, j^2))
end

@inline function _taylor_B(::Taylor{RefTriangle, order}, J::Tensors.Tensor{2, 2}) where {order}
    detJ = Tensors.det(J)
    z = 0.0
    if order == 1
        return (
            (1 / detJ, z, z),
            (z, J[1, 1], J[1, 2]),
            (z, J[2, 1], J[2, 2]),
        )
    end
    # Hessian block: Hₓ = A Ĥ Aᵀ with A = J⁻ᵀ; in (xx, xy, yy) components
    # M2 = [A11²      2A11A12          A12²    ;
    #       A11A21    A11A22 + A12A21  A12A22  ;
    #       A21²      2A21A22          A22²    ], and the DOF block is M2⁻ᵀ.
    A = Tensors.inv(J)'
    M2 = Tensors.Tensor{2, 3}(
        (
            A[1, 1]^2, A[1, 1] * A[2, 1], A[2, 1]^2,                                       # first column
            2 * A[1, 1] * A[1, 2], A[1, 1] * A[2, 2] + A[1, 2] * A[2, 1], 2 * A[2, 1] * A[2, 2],
            A[1, 2]^2, A[1, 2] * A[2, 2], A[2, 2]^2,
        )
    )
    B2 = Tensors.inv(M2)'
    return (
        (1 / detJ, z, z, z, z, z),
        (z, J[1, 1], J[1, 2], z, z, z),
        (z, J[2, 1], J[2, 2], z, z, z),
        (z, z, z, B2[1, 1], B2[1, 2], B2[1, 3]),
        (z, z, z, B2[2, 1], B2[2, 2], B2[2, 3]),
        (z, z, z, B2[3, 1], B2[3, 2], B2[3, 3]),
    )
end

@inline function Ferrite.apply_mapping!(
        funvals::Ferrite.FunctionValues{0}, ::TaylorMapping, q_point::Int, mapping_values, cell
    )
    B = _taylor_B(funvals.ip, Ferrite.getjacobian(mapping_values))
    N = Ferrite.getnbasefunctions(funvals)
    @inbounds for i in 1:N
        funvals.Nx[i, q_point] = sum(B[i][j] * funvals.Nξ[j, q_point] for j in 1:N)
    end
    return nothing
end

@inline function Ferrite.apply_mapping!(
        funvals::Ferrite.FunctionValues{1}, ::TaylorMapping, q_point::Int, mapping_values, cell
    )
    J = Ferrite.getjacobian(mapping_values)
    Jinv = Tensors.inv(J)
    B = _taylor_B(funvals.ip, J)
    N = Ferrite.getnbasefunctions(funvals)
    @inbounds for i in 1:N
        funvals.Nx[i, q_point] = sum(B[i][j] * funvals.Nξ[j, q_point] for j in 1:N)
        funvals.dNdx[i, q_point] = sum(B[i][j] * (funvals.dNdξ[j, q_point] ⋅ Jinv) for j in 1:N)
    end
    return nothing
end

# The mapping does not need the cell, keep cell-less reinit! working.
Ferrite.reinit_needs_cell(::Ferrite.CellValues{<:Ferrite.FunctionValues{<:Any, <:Taylor}}) = false
Ferrite.reinit_needs_cell(::Ferrite.FacetValues{<:Ferrite.FunctionValues{<:Any, <:Taylor}}) = false
