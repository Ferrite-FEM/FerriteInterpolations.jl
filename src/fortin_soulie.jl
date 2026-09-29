# Fortin-Soulie element (https://defelement.org/elements/fortin-soulie.html,
# DefElement: elements/fortin-soulie.def), implemented as in the original paper:
# M. Fortin and M. Soulie, "A non-conforming piecewise quadratic finite element
# on triangles", IJNME 19 (1983), DOI 10.1002/nme.1620190405.
#
# Cells/degrees implemented: RefTriangle degree 2 (the element only exists for
# degree 2). Conformity: L2 (nonconforming), identity mapping.
#
# The Fortin-Soulie space W_h consists of the piecewise quadratics that are
# continuous at the two Gauss-Legendre points of every interior edge (the patch
# test then holds: the jump across an edge is orthogonal to P1 there). The six
# Gauss points of a triangle lie on an ellipse, so they cannot be used as DOFs.
# Instead, following Proposition 1 of the paper, W_h = X_h + Phi_h, where X_h is
# the standard continuous P2 space and Phi_h contains one "neutral function" per
# triangle,
#
#     phi_0(x) = 2 - 3 (lambda_1^2 + lambda_2^2 + lambda_3^2),
#
# which vanishes at all six Gauss points and equals 1 at the centroid. So the
# interpolation is Ferrite's `Lagrange{RefTriangle, 2}` (DOFs 1-6, keeping their
# vertex/edge association and hence shared between cells as usual) plus phi_0 as
# a cell DOF (DOF 7). Adding phi_0 on a cell leaves the Gauss-point values
# unchanged, which is why the enriched space is only Gauss-Legendre continuous.
#
# The representation is not unique: dim(X_h ∩ Phi_h) = 1, since the sum of the
# neutral functions of all cells is itself continuous (its trace on every edge is
# the same quadratic, -1 at the vertices and 1/2 at the midpoint). Consequences
# (paper, p. 508):
#  * With Dirichlet conditions (on any part of the boundary) the representation
#    is unique, and the conditions are imposed on the X_h components, i.e. on the
#    vertex and midpoint DOFs -- which is what Ferrite's facet Dirichlet does
#    through the Lagrange entity DOFs.
#  * For pure Neumann problems, two X_h values must be fixed instead of one (one
#    at a vertex and one at a midpoint), since constants can be written in two
#    ways.
# Locally, the seven functions span P2 (dimension 6), so element matrices are
# singular on a single cell; this is not a Ciarlet element.
#
# Differences from DefElement (which cites the same paper): DefElement/symfem
# define a local element with six point-evaluation DOFs (two points on two
# edges, the midpoint of the third edge, and the centroid), which does not
# appear in the paper and cannot give the Fortin-Soulie space through DOF
# sharing. This file follows the paper, so it is not cross-checked against
# symfem's "Fortin-Soulie" element. (symfem up to version 2025.12 also used the
# points 1/3 and 2/3 along each edge instead of the Gauss points; this is fixed
# in https://github.com/mscroggs/symfem/pull/344.)
#
# `reference_coordinates` lists the Lagrange nodes and, for the neutral
# function, the centroid (as Ferrite's `BubbleEnrichedLagrange` does). Only the
# Lagrange nodes are used by facet Dirichlet conditions; nodal interpolation
# (e.g. `apply_analytical!`) is not an interpolant of this element, since the
# basis is not nodal.

"""
    FortinSoulie{RefTriangle, 2}()

Fortin-Soulie nonconforming quadratic element on the triangle, as constructed in
Fortin & Soulie (1983): continuous P2 (the six `Lagrange{RefTriangle, 2}` DOFs,
shared between cells) enriched with one neutral function per cell,
`2 - 3(λ₁² + λ₂² + λ₃²)`, which vanishes at the Gauss-Legendre points of the edges.
The resulting global space is the piecewise quadratics that are continuous at the
two Gauss-Legendre points of every interior edge.

Impose Dirichlet conditions as usual (they act on the vertex and midpoint DOFs).
For pure Neumann problems, fix one vertex and one midpoint value, since the
representation has one global degree of freedom too many.
"""
struct FortinSoulie{shape, order} <: ScalarInterpolation{shape, order} end

const _FortinSoulieP2 = Ferrite.Lagrange{RefTriangle, 2}()

Ferrite.conformity(::FortinSoulie) = Ferrite.L2Conformity()
Ferrite.adjust_dofs_during_distribution(::FortinSoulie{RefTriangle, 2}) =
    Ferrite.adjust_dofs_during_distribution(_FortinSoulieP2)

Ferrite.getnbasefunctions(::FortinSoulie{RefTriangle, 2}) = 7

function Ferrite.reference_shape_value(ip::FortinSoulie{RefTriangle, 2}, ξ::Vec{2}, i::Int)
    1 <= i <= 6 && return Ferrite.reference_shape_value(_FortinSoulieP2, ξ, i)
    if i == 7
        # Neutral function in barycentric coordinates (λ₁, λ₂, λ₃) = (ξ₁, ξ₂, 1 - ξ₁ - ξ₂)
        λ₁, λ₂ = ξ[1], ξ[2]
        λ₃ = 1 - λ₁ - λ₂
        return 2 - 3 * (λ₁^2 + λ₂^2 + λ₃^2)
    end
    return throw_out_of_range(ip, i)
end

Ferrite.vertexdof_indices(::FortinSoulie{RefTriangle, 2}) = Ferrite.vertexdof_indices(_FortinSoulieP2)
Ferrite.edgedof_interior_indices(::FortinSoulie{RefTriangle, 2}) = Ferrite.edgedof_interior_indices(_FortinSoulieP2)
Ferrite.facedof_interior_indices(::FortinSoulie{RefTriangle, 2}) = ((7,),)

function Ferrite.reference_coordinates(::FortinSoulie{RefTriangle, 2})
    return [Ferrite.reference_coordinates(_FortinSoulieP2)..., Vec{2, Float64}((1 / 3, 1 / 3))]
end
