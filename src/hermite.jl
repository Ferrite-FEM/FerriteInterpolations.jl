# Hermite element (cubic, C1-at-vertices) (https://defelement.org/elements/hermite.html,
# DefElement: elements/hermite.def).
#
# Cells/degrees implemented: RefTriangle, degree 3 (as defined). The interval
# version is `Ferrite.Hermite{RefLine, 3}`; the tetrahedron follows the same
# pattern with 3-blocks and four face-centroid value DOFs.
#
# DOFs: (u, du/dx, du/dy) at each vertex plus the value at the centroid. The
# vertex values and gradients are shared between neighboring cells, making the
# element C1 at the vertices and C0 across edges (the cubic trace on an edge
# is determined by the endpoint values and tangential derivatives, both shared
# through the full vertex gradients).
#
# DOF transformation (unblocks #2): the physical gradient DOFs relate to the
# reference ones through the cell Jacobian, so the mapped basis functions must
# be recombined per cell with the block-diagonal matrix
# M = blockdiag(1, J, 1, J, 1, J, 1) (2x2 vertex-gradient blocks; Kirby, "A
# general approach to transforming finite elements"). Unlike the diagonal M of
# the tensor-product Hermite elements this mixes basis functions, which is
# expressed through the openly dispatched dof-transformation stage from
# Ferrite-FEM/Ferrite.jl#1391: `HermiteDofTransformation` below implements its
# own `apply_dof_transformation!`, which `reinit!` applies after the mapping,
# uniformly to values/gradients/hessians (valid because M is cell-constant on
# affine cells -- and every `Lagrange{RefTriangle, 1}` triangle is affine, so
# unlike the tensor-product elements there is no axis-aligned restriction).
#
# The reference basis is the P3 dual basis of the reference DOFs (derived
# exactly by inverting the generalized Vandermonde matrix in rational
# arithmetic); it relates to symfem's Hermite basis by the vertex permutation
# (identical reference coordinates, see the test file).
#
# NOTE: `Ferrite` also exports a (different) `Hermite` type for the interval
# and tensor-product elements; with `using Ferrite, FerriteInterpolations` the
# unqualified name is ambiguous and one of them must be imported explicitly.

"""
    HermiteDofTransformation()

Cell-dependent (Kirby-style) dof transformation for the [`Hermite`](@ref)
triangle: the mapped basis function pair of each vertex-gradient dof pair is
recombined with the cell Jacobian `J` so that the dofs are the *physical*
gradient components at the vertices, consistently between neighboring cells.
Valid for affine geometries only.
"""
struct HermiteDofTransformation end

"""
    Hermite{RefTriangle, 3}()

The cubic Hermite element on the triangle: three dofs per vertex,
`(u, ∂u/∂x, ∂u/∂y)` (dof kinds `:value`, `:derivative_x`, `:derivative_y`),
plus the value at the centroid. The vertex gradients are shared between
neighboring cells, making the interpolation C¹ at the vertices and
C⁰-continuous across edges.

The gradient dofs are *physical* derivatives: the basis functions are
recombined per cell with [`HermiteDofTransformation`](@ref). All straight-edged
(affine) triangles are supported; embedded or curved elements are not.

`Dirichlet` conditions on facet- or vertex-sets constrain the value dofs by
default; the gradient dofs are selected with the `kind` keyword
(`:derivative_x`, `:derivative_y`).

!!! note
    `Ferrite` exports a different `Hermite` (the interval element); qualify
    the name when both packages are loaded.
"""
struct Hermite{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::Hermite) = Ferrite.H1Conformity()
Ferrite.adjust_dofs_during_distribution(::Hermite) = false
Ferrite.dof_transformation(::Hermite) = HermiteDofTransformation()

Ferrite.getnbasefunctions(::Hermite{RefTriangle, 3}) = 10
Ferrite.vertexdof_indices(::Hermite{RefTriangle, 3}) = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
Ferrite.edgedof_interior_indices(::Hermite{RefTriangle, 3}) = ((), (), ())
Ferrite.facedof_interior_indices(::Hermite{RefTriangle, 3}) = ((10,),)
function Ferrite.dof_kinds(::Hermite{RefTriangle, 3})
    return (
        :value, :derivative_x, :derivative_y,
        :value, :derivative_x, :derivative_y,
        :value, :derivative_x, :derivative_y,
        :value,
    )
end

function Ferrite.reference_coordinates(::Hermite{RefTriangle, 3})
    return [
        Vec{2, Float64}((1.0, 0.0)),
        Vec{2, Float64}((1.0, 0.0)),
        Vec{2, Float64}((1.0, 0.0)),
        Vec{2, Float64}((0.0, 1.0)),
        Vec{2, Float64}((0.0, 1.0)),
        Vec{2, Float64}((0.0, 1.0)),
        Vec{2, Float64}((0.0, 0.0)),
        Vec{2, Float64}((0.0, 0.0)),
        Vec{2, Float64}((0.0, 0.0)),
        Vec{2, Float64}((1 / 3, 1 / 3)),
    ]
end

function Ferrite.reference_shape_value(ip::Hermite{RefTriangle, 3}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    # Vertex 1 = (1, 0): value, d/dx, d/dy
    i == 1 && return 3x^2 - 7x * y - 2x^3 + 7x^2 * y + 7x * y^2
    i == 2 && return -x^2 + 2x * y + x^3 - 2x^2 * y - 2x * y^2
    i == 3 && return -x * y + 2x^2 * y + x * y^2
    # Vertex 2 = (0, 1): value, d/dx, d/dy
    i == 4 && return 3y^2 - 7x * y + 7x^2 * y + 7x * y^2 - 2y^3
    i == 5 && return -x * y + x^2 * y + 2x * y^2
    i == 6 && return -y^2 + 2x * y - 2x^2 * y - 2x * y^2 + y^3
    # Vertex 3 = (0, 0): value, d/dx, d/dy
    i == 7 && return 1 - 3x^2 - 13x * y - 3y^2 + 2x^3 + 13x^2 * y + 13x * y^2 + 2y^3
    i == 8 && return x - 2x^2 - 3x * y + x^3 + 3x^2 * y + 2x * y^2
    i == 9 && return y - 3x * y - 2y^2 + 2x^2 * y + 3x * y^2 + y^3
    # Centroid value
    i == 10 && return 27x * y * (1 - x - y)
    return throw_out_of_range(ip, i)
end

function Ferrite.check_geometry_compatibility(ip::Hermite, ip_geo::Ferrite.VectorizedInterpolation{sdim}) where {sdim}
    if !(sdim == 2 && ip_geo.ip isa Ferrite.Lagrange{RefTriangle, 1})
        throw(
            ArgumentError(
                "$(ip) requires an affine 2D geometry (ip_geo = Lagrange{RefTriangle, 1}()^2), got $(ip_geo). " *
                    "Embedded (sdim > 2) or curved elements are not supported."
            )
        )
    end
    return nothing
end

# ---- the dof transformation ---------------------------------------------------

# The Jacobian is needed for the vertex-gradient blocks (but the cell is not).
Ferrite.required_geo_diff_order(::HermiteDofTransformation) = 1
Ferrite.dof_transformation_needs_cell(::HermiteDofTransformation) = false

# Recombine the basis function pair (i1, i2) = (3v + 2, 3v + 3) of each
# vertex-gradient dof pair with the 2x2 block J. Duality with the physical
# gradient functionals: applying grad_x = J^-T grad_ξ to the recombined pair
# gives J J⁻¹ = I at the vertex. M is cell-constant, so the same recombination
# applies to values, gradients and hessians.
@inline function _apply_vertex_gradient_blocks!(A, q_point::Int, J::Tensors.Tensor{2, 2})
    @inbounds for v in 0:2
        i1, i2 = 3v + 2, 3v + 3
        a1, a2 = A[i1, q_point], A[i2, q_point]
        A[i1, q_point] = J[1, 1] * a1 + J[1, 2] * a2
        A[i2, q_point] = J[2, 1] * a1 + J[2, 2] * a2
    end
    return nothing
end

@inline function Ferrite.apply_dof_transformation!(
        funvals::Ferrite.FunctionValues{0}, ::HermiteDofTransformation, q_point::Int, mapping_values, args...
    )
    J = Ferrite.getjacobian(mapping_values)
    _apply_vertex_gradient_blocks!(funvals.Nx, q_point, J)
    return nothing
end

@inline function Ferrite.apply_dof_transformation!(
        funvals::Ferrite.FunctionValues{1}, ::HermiteDofTransformation, q_point::Int, mapping_values, args...
    )
    J = Ferrite.getjacobian(mapping_values)
    _apply_vertex_gradient_blocks!(funvals.Nx, q_point, J)
    _apply_vertex_gradient_blocks!(funvals.dNdx, q_point, J)
    return nothing
end

@inline function Ferrite.apply_dof_transformation!(
        funvals::Ferrite.FunctionValues{2}, ::HermiteDofTransformation, q_point::Int, mapping_values, args...
    )
    J = Ferrite.getjacobian(mapping_values)
    _apply_vertex_gradient_blocks!(funvals.Nx, q_point, J)
    _apply_vertex_gradient_blocks!(funvals.dNdx, q_point, J)
    _apply_vertex_gradient_blocks!(funvals.d2Ndx2, q_point, J)
    return nothing
end
