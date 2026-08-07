# Bernstein element (https://defelement.org/elements/bernstein.html,
# DefElement: elements/bernstein.def).
#
# Cells/degrees implemented: RefLine 1-3, RefTriangle 1-3, RefTetrahedron 1-2.
# Conformity: H1, identity mapping.
#
# DOF layout: one DOF per vertex (point evaluation), then Bernstein-coefficient
# DOFs on edge interiors (ordered from the first towards the second vertex of
# the reference edge) and cell interiors. The interior DOFs are not point
# evaluations, so `Ferrite.reference_coordinates` is deliberately not defined:
# the element cannot be used with Ferrite's nodal machinery (`Dirichlet` via
# `BCValues`, `apply_analytical!`, `L2Projector`) for degree >= 2.
#
# Basis: Bernstein polynomials B_α = (|α|! / α!) λ^α over the barycentric
# coordinates λ of the cell, |α| = degree (Ainsworth, Andriamaro, Davydov,
# "Bernstein-Bezier finite elements of arbitrary order and optimal assembly
# procedures", SIAM J. Sci. Comput. 33 (2011), DOI 10.1137/11082539X).
# Cross-checked against symfem in test/test_bernstein.jl.
#
# NOT implemented (blocked upstream in Ferrite):
# - RefTetrahedron degree >= 3: two or more interior DOFs per (shared) face
#   need face-orientation permutations; Ferrite's `interior_facedofs_on_lattice`
#   machinery is specific to the nodal Lagrange lattice and does not describe
#   Bernstein coefficient DOFs.
# Higher degrees on RefLine/RefTriangle are straightforward extensions of the
# tables below (edge-interior DOFs reverse correctly under Ferrite's
# `adjust_dofs_during_distribution` edge-reversal semantics).

"""
    Bernstein{shape, order}()

Bernstein(-Bezier) element of degree `order` on `shape`. H1-conforming.
Vertex DOFs are point evaluations; the remaining DOFs are Bernstein
coefficients (not point evaluations).
"""
struct Bernstein{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::Bernstein) = Ferrite.H1Conformity()

Ferrite.adjust_dofs_during_distribution(::Bernstein) = true
Ferrite.adjust_dofs_during_distribution(::Bernstein{<:Any, 1}) = false
Ferrite.adjust_dofs_during_distribution(::Bernstein{<:Any, 2}) = false

# Barycentric coordinates of ξ, one entry per reference vertex in Ferrite's
# vertex order (see src/docs.jl in Ferrite for the coordinate conventions).
_barycentrics(::Type{RefLine}, ξ::Vec{1}) = ((1 - ξ[1]) / 2, (1 + ξ[1]) / 2)
_barycentrics(::Type{RefTriangle}, ξ::Vec{2}) = (ξ[1], ξ[2], 1 - ξ[1] - ξ[2])
_barycentrics(::Type{RefTetrahedron}, ξ::Vec{3}) = (1 - ξ[1] - ξ[2] - ξ[3], ξ[1], ξ[2], ξ[3])

# (multinomial coefficient, barycentric exponent tuple) for each shape function,
# in Ferrite DOF order (vertices -> edge interiors -> face/volume interiors,
# edge DOFs from the first towards the second vertex of the reference edge).
_bernstein_data(::Bernstein{RefLine, 1}) = ((1, (1, 0)), (1, (0, 1)))
_bernstein_data(::Bernstein{RefLine, 2}) = ((1, (2, 0)), (1, (0, 2)), (2, (1, 1)))
_bernstein_data(::Bernstein{RefLine, 3}) = (
    (1, (3, 0)), (1, (0, 3)),                # vertices
    (3, (2, 1)), (3, (1, 2)),                # edge (1, 2)
)

_bernstein_data(::Bernstein{RefTriangle, 1}) = ((1, (1, 0, 0)), (1, (0, 1, 0)), (1, (0, 0, 1)))
_bernstein_data(::Bernstein{RefTriangle, 2}) = (
    (1, (2, 0, 0)), (1, (0, 2, 0)), (1, (0, 0, 2)),  # vertices
    (2, (1, 1, 0)),                                  # edge (1, 2)
    (2, (0, 1, 1)),                                  # edge (2, 3)
    (2, (1, 0, 1)),                                  # edge (3, 1)
)
_bernstein_data(::Bernstein{RefTriangle, 3}) = (
    (1, (3, 0, 0)), (1, (0, 3, 0)), (1, (0, 0, 3)),  # vertices
    (3, (2, 1, 0)), (3, (1, 2, 0)),                  # edge (1, 2)
    (3, (0, 2, 1)), (3, (0, 1, 2)),                  # edge (2, 3)
    (3, (1, 0, 2)), (3, (2, 0, 1)),                  # edge (3, 1)
    (6, (1, 1, 1)),                                  # face interior
)

_bernstein_data(::Bernstein{RefTetrahedron, 1}) = (
    (1, (1, 0, 0, 0)), (1, (0, 1, 0, 0)), (1, (0, 0, 1, 0)), (1, (0, 0, 0, 1)),
)
_bernstein_data(::Bernstein{RefTetrahedron, 2}) = (
    (1, (2, 0, 0, 0)), (1, (0, 2, 0, 0)), (1, (0, 0, 2, 0)), (1, (0, 0, 0, 2)),  # vertices
    (2, (1, 1, 0, 0)),  # edge (1, 2)
    (2, (0, 1, 1, 0)),  # edge (2, 3)
    (2, (1, 0, 1, 0)),  # edge (3, 1)
    (2, (1, 0, 0, 1)),  # edge (1, 4)
    (2, (0, 1, 0, 1)),  # edge (2, 4)
    (2, (0, 0, 1, 1)),  # edge (3, 4)
)

Ferrite.getnbasefunctions(ip::Bernstein) = length(_bernstein_data(ip))

function Ferrite.reference_shape_value(ip::Bernstein{shape}, ξ::Vec, i::Int) where {shape}
    data = _bernstein_data(ip)
    1 <= i <= length(data) || throw_out_of_range(ip, i)
    c, α = data[i]
    λ = _barycentrics(shape, ξ)
    return c * prod(λ .^ α)
end

# RefLine
Ferrite.vertexdof_indices(::Bernstein{RefLine}) = ((1,), (2,))
Ferrite.edgedof_interior_indices(::Bernstein{RefLine, 1}) = (((),))
Ferrite.edgedof_interior_indices(::Bernstein{RefLine, 2}) = ((3,),)
Ferrite.edgedof_interior_indices(::Bernstein{RefLine, 3}) = ((3, 4),)
Ferrite.facedof_interior_indices(::Bernstein{RefLine}) = ()
Ferrite.volumedof_interior_indices(::Bernstein{RefLine}) = ()

# RefTriangle
Ferrite.vertexdof_indices(::Bernstein{RefTriangle}) = ((1,), (2,), (3,))
Ferrite.edgedof_interior_indices(::Bernstein{RefTriangle, 1}) = ((), (), ())
Ferrite.edgedof_interior_indices(::Bernstein{RefTriangle, 2}) = ((4,), (5,), (6,))
Ferrite.edgedof_interior_indices(::Bernstein{RefTriangle, 3}) = ((4, 5), (6, 7), (8, 9))
Ferrite.facedof_interior_indices(::Bernstein{RefTriangle, 1}) = ((),)
Ferrite.facedof_interior_indices(::Bernstein{RefTriangle, 2}) = ((),)
Ferrite.facedof_interior_indices(::Bernstein{RefTriangle, 3}) = ((10,),)
Ferrite.volumedof_interior_indices(::Bernstein{RefTriangle}) = ()

# RefTetrahedron
Ferrite.vertexdof_indices(::Bernstein{RefTetrahedron}) = ((1,), (2,), (3,), (4,))
Ferrite.edgedof_interior_indices(::Bernstein{RefTetrahedron, 1}) = ((), (), (), (), (), ())
Ferrite.edgedof_interior_indices(::Bernstein{RefTetrahedron, 2}) = ((5,), (6,), (7,), (8,), (9,), (10,))
Ferrite.facedof_interior_indices(::Bernstein{RefTetrahedron, 1}) = ((), (), (), ())
Ferrite.facedof_interior_indices(::Bernstein{RefTetrahedron, 2}) = ((), (), (), ())
Ferrite.volumedof_interior_indices(::Bernstein{RefTetrahedron}) = ()
