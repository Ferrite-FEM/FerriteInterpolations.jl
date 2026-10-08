# Crouzeix-Falk element, implemented as in the original paper: M. Crouzeix and
# R. S. Falk, "Nonconforming finite elements for the Stokes problem", Math.
# Comp. 52 (1989), DOI 10.2307/2008475 (see also
# https://defelement.org/elements/crouzeix-falk.html).
#
# Cells/degrees implemented: RefTriangle degree 3 (the element only exists for
# degree 3). Conformity: L2 (nonconforming), identity mapping. Space: full P3.
#
# The paper (Theorem 2.1) defines the global space as the piecewise cubics that
# are continuous at the three Gauss-Legendre points of every edge (and, for
# homogeneous Dirichlet conditions, vanish at the Gauss points on the
# boundary). The ten DOFs are therefore point evaluations at the three
# Gauss-Legendre points of every edge, at s = 1/2 - sqrt(15)/10, 1/2 and
# 1/2 + sqrt(15)/10 along the edge (ordered from the first towards the second
# vertex of the reference edge), plus the centroid. The jump across an edge then
# is a multiple of the cubic Legendre polynomial and hence orthogonal to P2 on
# the edge, which gives the patch test and the O(h^3) energy-norm estimate of
# the paper.
#
# Differences from DefElement: DefElement/symfem place the edge DOFs at the
# equispaced points 1/4, 1/2 and 3/4. The jump across an edge is then only
# orthogonal to the even polynomials, the patch test fails for quadratic
# solutions and the method loses two orders of convergence. This file follows
# the paper, so it is not cross-checked against symfem's "Crouzeix-Falk".
#
# The point set is symmetric under edge reversal, so the edge DOFs are true
# edge DOFs (shared between neighboring cells, like Ferrite's CrouzeixRaviart),
# with `adjust_dofs_during_distribution` handling the order reversal on edges
# with opposite local orientation. Interior traces only agree at the three
# shared points per edge -- that is the element's nonconformity, hence
# L2Conformity. Dirichlet conditions through Ferrite's `Dirichlet` act on the
# edge DOFs, i.e. at the Gauss points of the boundary edges, as in the paper.
#
# All DOFs are point evaluations, so `reference_coordinates` is defined.
# DOF order (vertices carry no DOFs): edge 1: (1, 2, 3); edge 2: (4, 5, 6);
# edge 3: (7, 8, 9); centroid: (10,).
#
# Basis: the dual basis of the DOFs above, derived with sympy (the DOFs are
# unisolvent on P3: the Vandermonde determinant in the monomial basis is
# 3 sqrt(15) / 8e8).

"""
    CrouzeixFalk{shape, order}()

Crouzeix-Falk element on the triangle, degree 3, as in Crouzeix & Falk (1989):
full P3 with point-evaluation DOFs at the three Gauss-Legendre points of every
edge (shared between neighboring cells) and at the centroid. Nonconforming
(L2): the global space is continuous at the Gauss points of the edges only.
Dirichlet conditions are imposed at the Gauss points of the boundary edges.
"""
struct CrouzeixFalk{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::CrouzeixFalk) = Ferrite.L2Conformity()
Ferrite.adjust_dofs_during_distribution(::CrouzeixFalk) = true

Ferrite.getnbasefunctions(::CrouzeixFalk{RefTriangle, 3}) = 10

function Ferrite.reference_shape_value(ip::CrouzeixFalk{RefTriangle, 3}, ξ::Vec{2, T}, i::Int) where {T}
    x, y = ξ[1], ξ[2]
    r = sqrt(T(15))
    i == 1 && return 5 * (20x^3 + 2r * x^2 * y + 37x^2 * y - 30x^2 - 2r * x * y^2 + 37x * y^2 - 41x * y + 12x + 20y^3 - 30y^2 + 12y - 1) / 6
    i == 2 && return -2 * (20x^3 + 19x^2 * y - 30x^2 + 19x * y^2 - 29x * y + 12x + 20y^3 - 30y^2 + 12y - 1) / 3
    i == 3 && return 5 * (20x^3 - 2r * x^2 * y + 37x^2 * y - 30x^2 + 2r * x * y^2 + 37x * y^2 - 41x * y + 12x + 20y^3 - 30y^2 + 12y - 1) / 6
    i == 4 && return -5 * (20x^3 + 2r * x^2 * y + 23x^2 * y - 30x^2 + 23x * y^2 + 6r * x * y^2 - 27x * y - 4r * x * y + 12x + 4r * y^3 - 6r * y^2 - 4y^2 + 4y + 2r * y - 1) / 6
    i == 5 && return 2 * (20x^3 + 41x^2 * y - 30x^2 + 41x * y^2 - 51x * y + 12x - 10y^2 + 10y - 1) / 3
    i == 6 && return -5 * (20x^3 - 2r * x^2 * y + 23x^2 * y - 30x^2 - 6r * x * y^2 + 23x * y^2 - 27x * y + 4r * x * y + 12x - 4r * y^3 - 4y^2 + 6r * y^2 - 2r * y + 4y - 1) / 6
    i == 7 && return 5 * (4r * x^3 - 23x^2 * y + 6r * x^2 * y - 6r * x^2 + 4x^2 - 23x * y^2 + 2r * x * y^2 - 4r * x * y + 27x * y - 4x + 2r * x - 20y^3 + 30y^2 - 12y + 1) / 6
    i == 8 && return 2 * (41x^2 * y - 10x^2 + 41x * y^2 - 51x * y + 10x + 20y^3 - 30y^2 + 12y - 1) / 3
    i == 9 && return -5 * (4r * x^3 + 23x^2 * y + 6r * x^2 * y - 6r * x^2 - 4x^2 + 2r * x * y^2 + 23x * y^2 - 27x * y - 4r * x * y + 4x + 2r * x + 20y^3 - 30y^2 + 12y - 1) / 6
    i == 10 && return -27x * y * (x + y - 1)
    return throw_out_of_range(ip, i)
end

# Gauss-Legendre points of the reference edge parametrization s in [0, 1].
const _CF_GAUSS = (1 / 2 - sqrt(15) / 10, 1 / 2, 1 / 2 + sqrt(15) / 10)

function Ferrite.reference_coordinates(::CrouzeixFalk{RefTriangle, 3})
    return [
        (Vec((1 - s, s)) for s in _CF_GAUSS)..., # edge 1
        (Vec((0.0, 1 - s)) for s in _CF_GAUSS)..., # edge 2
        (Vec((s, 0.0)) for s in _CF_GAUSS)..., # edge 3
        Vec((1 / 3, 1 / 3)), # centroid
    ]
end

Ferrite.vertexdof_indices(::CrouzeixFalk{RefTriangle, 3}) = ((), (), ())
Ferrite.edgedof_interior_indices(::CrouzeixFalk{RefTriangle, 3}) = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
Ferrite.facedof_interior_indices(::CrouzeixFalk{RefTriangle, 3}) = ((10,),)
Ferrite.volumedof_interior_indices(::CrouzeixFalk{RefTriangle, 3}) = ()
