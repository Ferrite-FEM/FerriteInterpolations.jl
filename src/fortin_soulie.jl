# Fortin-Soulie element (https://defelement.org/elements/fortin-soulie.html,
# DefElement: elements/fortin-soulie.def).
#
# Cells/degrees implemented: RefTriangle degree 2 (the element only exists for
# degree 2). Conformity: L2 (nonconforming), identity mapping.
#
# The six point-evaluation DOFs at the two Gauss-like third-points of every
# edge are linearly dependent on P2 (the nonconforming P2 bubble vanishes at
# all six), so the element is necessarily asymmetric (Fortin & Soulie, "A
# non-conforming piecewise quadratic finite element on triangles", IJNME 19
# (1983), DOI 10.1002/nme.1620190405): two point evaluations on edges 1 and 2
# (at the third-points, ordered from the first towards the second vertex of
# the reference edge), ONE on edge 3 (its midpoint), and one at the centroid.
#
# Following DiscontinuousLagrange, all DOFs are cell DOFs
# (volumedof_interior_indices): in Ferrite, entity-attached DOFs are
# identified across neighboring cells, and the Fortin-Soulie global
# nonconforming space (agreement at two Gauss points per interior edge)
# cannot be realized that way -- each triangle has exactly one single-DOF
# edge, so on a general mesh some shared edge would pair two DOFs of one
# cell with one of the other. Realizing the coupled Fortin-Soulie space
# would need mesh-level machinery (a mesh-dependent choice of local DOF
# layout) that Ferrite does not have; with cell DOFs the element is usable
# as a broken/DG-style basis, with any coupling imposed weakly. The
# geometric edge association of the DOFs is exposed through
# `dirichlet_edgedof_indices` so that facet Dirichlet BCs work.
#
# All DOFs are point evaluations, so `reference_coordinates` is defined.
# DOF order: edge 1 points: (1, 2); edge 2 points: (3, 4); edge 3 point: (5,);
# centroid: (6,).
#
# Basis: transcribed from symfem ("Fortin-Soulie", degree 2); the reference
# triangles agree pointwise (same coordinates, different vertex numbering),
# see test/test_fortin_soulie.jl for the cross-check.

"""
    FortinSoulie{shape, order}()

Fortin-Soulie element on the triangle, degree 2. Nonconforming (L2), with two
point-evaluation DOFs on the first two edges, one on the third, and one at the
centroid. All DOFs are cell DOFs (not shared between cells); see the source
file for why the coupled Fortin-Soulie space is not representable in Ferrite.
"""
struct FortinSoulie{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::FortinSoulie) = Ferrite.L2Conformity()
Ferrite.adjust_dofs_during_distribution(::FortinSoulie) = false

Ferrite.getnbasefunctions(::FortinSoulie{RefTriangle, 2}) = 6

function Ferrite.reference_shape_value(ip::FortinSoulie{RefTriangle, 2}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return 9x^2 / 2 - 3x / 2 - 27y^2 / 16 + 27y / 16 - 3 // 8
    i == 2 && return 9x * y - 3x + 27y^2 / 4 - 27y / 4 + 3 // 2
    i == 3 && return -9x * y + 3x - 9y^2 / 4 + 21y / 4 - 3 // 2
    i == 4 && return 9x^2 / 2 + 9x * y - 15x / 2 + 45y^2 / 16 - 93y / 16 + 21 // 8
    i == 5 && return 9y^2 / 2 - 9y / 2 + 1
    i == 6 && return -9x^2 - 9x * y + 9x - 81y^2 / 8 + 81y / 8 - 9 // 4
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::FortinSoulie{RefTriangle, 2})
    return [
        Vec((2 / 3, 1 / 3)), Vec((1 / 3, 2 / 3)), # on edge 1
        Vec((0.0, 2 / 3)), Vec((0.0, 1 / 3)),     # on edge 2
        Vec((1 / 2, 0.0)),                        # on edge 3
        Vec((1 / 3, 1 / 3)),                      # centroid
    ]
end

# All DOFs in the cell interior (not shared between cells), like
# DiscontinuousLagrange.
Ferrite.volumedof_interior_indices(ip::FortinSoulie) = ntuple(i -> i, Ferrite.getnbasefunctions(ip))

# Geometric edge association of the point-evaluation DOFs, for facet
# Dirichlet BCs.
Ferrite.dirichlet_edgedof_indices(::FortinSoulie{RefTriangle, 2}) = ((1, 2), (3, 4), (5,))
