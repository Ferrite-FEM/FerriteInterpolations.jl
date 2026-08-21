# Bogner-Fox-Schmit element (C1 tensor Hermite on quadrilaterals) (https://defelement.org/elements/bogner-fox-schmitt.html,
# DefElement: elements/bogner-fox-schmitt.def).
#
# Cells/degrees implemented: RefQuadrilateral, degree 3 (as defined): the
# tensor product of the cubic Hermite line with four DOFs per vertex,
# (u, du/dx, du/dy, d2u/dxdy).
#
# DOF transformation (unblocks #12/#2): the physical derivative DOFs relate to
# the reference ones by a cell-dependent matrix M(J) (Kirby, "A general
# approach to transforming finite elements"), which is coordinated across the
# mesh through the openly dispatched dof-transformation stage from
# Ferrite-FEM/Ferrite.jl#1391. On axis-aligned rectangles (diagonal J) M is
# diagonal, so the element reuses Ferrite's `DiagonalDofTransformation` and
# only declares the per-basis-function scaling in `get_dof_scaling`:
# (1, J11, J22, J11*J22) per vertex.
#
# The restriction to axis-aligned rectangles is fundamental, not an
# implementation limit: for a non-diagonal Jacobian the chain rule for the
# mixed DOF d2u/dxi1 dxi2 drags in second derivatives that are not DOFs, so no
# DOF transformation matrix exists (and C1 coupling between neighboring cells
# is lost). The type-level part of the guard lives in
# `check_geometry_compatibility` (bilinear non-embedded geometry only); the
# per-cell rectangularity check runs in `get_dof_scaling` at `reinit!` time.

"""
    BognerFoxSchmitt()

The Bogner-Fox-Schmit element on the quadrilateral: the C1-continuous tensor
product of the cubic Hermite line, with four dofs per vertex
`(u, ∂u/∂x, ∂u/∂y, ∂²u/∂x∂y)` (dof kinds `:value`, `:derivative_x`,
`:derivative_y`, `:derivative_xy`), usable for fourth-order problems such as
Kirchhoff plates or the Cahn-Hilliard equation in primal form.

The derivative dofs are *physical* derivatives (the corresponding basis
functions are scaled by the geometric Jacobian). Only axis-aligned rectangular
`Quadrilateral` cells are supported; C¹-continuity across cell boundaries
would otherwise be lost.

`Dirichlet` conditions on facet- or vertex-sets constrain the value dofs by
default; the derivative dofs are selected with the `kind` keyword, e.g. for
zero normal derivative ``∂u/∂n = 0`` on a left/right boundary of a rectangular
grid
```julia
add!(ch, Dirichlet(:u, ∂Ω, Returns(0.0); kind = :derivative_x))
add!(ch, Dirichlet(:u, ∂Ω, Returns(0.0); kind = :derivative_xy))
```
where clamping also the mixed dof makes the normal derivative vanish along the
whole facet, not only at the vertices.
"""
struct BognerFoxSchmitt <: ScalarInterpolation{RefQuadrilateral, 3} end

Ferrite.conformity(::BognerFoxSchmitt) = Ferrite.H1Conformity() # C1 ⊂ C0; finer classification is not needed
Ferrite.adjust_dofs_during_distribution(::BognerFoxSchmitt) = false
Ferrite.dof_transformation(::BognerFoxSchmitt) = Ferrite.DiagonalDofTransformation()

Ferrite.getnbasefunctions(::BognerFoxSchmitt) = 16
Ferrite.vertexdof_indices(::BognerFoxSchmitt) = ((1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16))
Ferrite.edgedof_interior_indices(::BognerFoxSchmitt) = ((), (), (), ())
Ferrite.facedof_interior_indices(::BognerFoxSchmitt) = ((),)
function Ferrite.dof_kinds(::BognerFoxSchmitt)
    return (
        :value, :derivative_x, :derivative_y, :derivative_xy,
        :value, :derivative_x, :derivative_y, :derivative_xy,
        :value, :derivative_x, :derivative_y, :derivative_xy,
        :value, :derivative_x, :derivative_y, :derivative_xy,
    )
end

function Ferrite.reference_coordinates(::BognerFoxSchmitt)
    return [
        x for x in (
                Vec{2, Float64}((-1.0, -1.0)),
                Vec{2, Float64}((1.0, -1.0)),
                Vec{2, Float64}((1.0, 1.0)),
                Vec{2, Float64}((-1.0, 1.0)),
            ) for _ in 1:4
    ]
end

# 1D cubic Hermite basis on [-1, 1]: value/derivative at x = -1 (i = 1, 2) and
# at x = +1 (i = 3, 4). Identical to Ferrite's `Hermite{RefLine, 3}`.
function _bfs_hermite1d(x, i::Int)
    i == 1 && return (1 - x)^2 * (2 + x) / 4
    i == 2 && return (1 - x)^2 * (1 + x) / 4
    i == 3 && return (1 + x)^2 * (2 - x) / 4
    return -(1 + x)^2 * (1 - x) / 4
end

function Ferrite.reference_shape_value(ip::BognerFoxSchmitt, ξ::Vec{2}, i::Int)
    1 <= i <= 16 || throw_out_of_range(ip, i)
    v, k = divrem(i - 1, 4) # vertex (0-based) and dof kind: 0 value, 1 ∂ξ₁, 2 ∂ξ₂, 3 ∂ξ₁∂ξ₂
    i1 = ((v == 1 || v == 2) ? 2 : 0) + ((k == 1 || k == 3) ? 2 : 1)
    i2 = (v >= 2 ? 2 : 0) + ((k == 2 || k == 3) ? 2 : 1)
    return _bfs_hermite1d(ξ[1], i1) * _bfs_hermite1d(ξ[2], i2)
end

function Ferrite.check_geometry_compatibility(ip::BognerFoxSchmitt, ip_geo::Ferrite.VectorizedInterpolation{sdim}) where {sdim}
    if !(sdim == 2 && ip_geo.ip isa Ferrite.Lagrange{RefQuadrilateral, 1})
        throw(
            ArgumentError(
                "$(ip) requires a 2D bilinear geometry (ip_geo = Lagrange{RefQuadrilateral, 1}()^2), got $(ip_geo). " *
                    "Embedded (sdim > 2) or higher-order quadrilateral elements are not supported."
            )
        )
    end
    return nothing
end

# Type-level checks (`check_geometry_compatibility`) cannot rule out all unsupported
# geometries: a bilinear quadrilateral must additionally be an axis-aligned rectangle
# for the dof scaling to be cell-constant (diagonal Jacobian). Checked here since
# `get_dof_scaling` is the per-cell entry point of the transformation. (For a bilinear
# cell, off-diagonal Jacobian entries that vanish at two or more quadrature points per
# direction vanish everywhere.)
@inline function _check_rectangular(ip::BognerFoxSchmitt, J::Tensors.Tensor{2, 2})
    # √eps relative tolerance: round-off in the node coordinates is absolute in the
    # coordinate magnitude and can exceed eps(J) for fine grids (small J entries), while
    # genuinely sheared/rotated cells have off-diagonal entries comparable to the diagonal.
    tol = sqrt(eps(typeof(J[1, 1]))) * max(abs(J[1, 1]), abs(J[2, 2]))
    if abs(J[1, 2]) > tol || abs(J[2, 1]) > tol
        _throw_nonrectangular_cell(ip, J)
    end
    return nothing
end
@noinline function _throw_nonrectangular_cell(ip::BognerFoxSchmitt, J::Tensors.Tensor{2})
    throw(ArgumentError("$(ip) requires axis-aligned rectangular cells (diagonal Jacobian), got J = $(J)"))
end

@inline function Ferrite.get_dof_scaling(ip::BognerFoxSchmitt, shape_nr::Int, J::Tensors.Tensor{2, 2})
    _check_rectangular(ip, J)
    k = (shape_nr - 1) % 4
    k == 0 && return one(J[1, 1])
    k == 1 && return J[1, 1]
    k == 2 && return J[2, 2]
    return J[1, 1] * J[2, 2]
end
