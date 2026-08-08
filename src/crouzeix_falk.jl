# Crouzeix-Falk element (https://defelement.org/elements/crouzeix-falk.html,
# DefElement: elements/crouzeix-falk.def).
#
# Cells/degrees implemented: RefTriangle degree 3 (the element only exists for
# degree 3). Conformity: L2 (nonconforming), identity mapping. Space: full P3.
#
# The ten DOFs are point evaluations at three equispaced interior points per
# edge (at 1/4, 1/2 and 3/4 along the edge, ordered from the first towards the
# second vertex of the reference edge) plus the centroid (Crouzeix & Falk,
# "Nonconforming finite elements for the Stokes problem", Math. Comp. 52
# (1989), DOI 10.2307/2008475).
#
# Unlike Fortin-Soulie, the DOF layout is edge-symmetric (every edge carries
# three DOFs, and the point set {1/4, 1/2, 3/4} maps onto itself under edge
# reversal), so the coupled Crouzeix-Falk space IS realizable through
# Ferrite's entity-based DOF identification: the edge DOFs are true edge DOFs
# (shared between neighboring cells, like Ferrite's CrouzeixRaviart), with
# `adjust_dofs_during_distribution` handling the order reversal on edges with
# opposite local orientation. Interior traces still only agree at the three
# shared points per edge -- that is the element's nonconformity, hence
# L2Conformity.
#
# All DOFs are point evaluations, so `reference_coordinates` is defined.
# DOF order (vertices carry no DOFs): edge 1: (1, 2, 3); edge 2: (4, 5, 6);
# edge 3: (7, 8, 9); centroid: (10,).
#
# Basis: transcribed from symfem ("Crouzeix-Falk", degree 3); the reference
# triangles agree pointwise (same coordinates, different vertex/edge
# numbering), see test/test_crouzeix_falk.jl for the cross-check.

"""
    CrouzeixFalk{shape, order}()

Crouzeix-Falk element on the triangle, degree 3. Nonconforming (L2), with
three point-evaluation DOFs per edge (shared between neighboring cells) and
one at the centroid.
"""
struct CrouzeixFalk{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::CrouzeixFalk) = Ferrite.L2Conformity()
Ferrite.adjust_dofs_during_distribution(::CrouzeixFalk) = true

Ferrite.getnbasefunctions(::CrouzeixFalk{RefTriangle, 3}) = 10

function Ferrite.reference_shape_value(ip::CrouzeixFalk{RefTriangle, 3}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return 64x^3 / 3 + 42x^2 * y - 32x^2 + 94x * y^2 / 3 - 134x * y / 3 + 44x / 3 + 64y^3 / 3 - 32y^2 + 44y / 3 - 2
    i == 2 && return -32x^3 - 43x^2 * y + 48x^2 - 43x * y^2 + 59x * y - 22x - 32y^3 + 48y^2 - 22y + 3
    i == 3 && return 64x^3 / 3 + 94x^2 * y / 3 - 32x^2 + 42x * y^2 - 134x * y / 3 + 44x / 3 + 64y^3 / 3 - 32y^2 + 44y / 3 - 2
    i == 4 && return -64x^3 / 3 - 98x^2 * y / 3 + 32x^2 - 130x * y^2 / 3 + 46x * y - 44x / 3 - 32y^3 / 3 + 24y^2 - 40y / 3 + 2
    i == 5 && return 32x^3 + 53x^2 * y - 48x^2 + 53x * y^2 - 69x * y + 22x - 16y^2 + 16y - 3
    i == 6 && return -64x^3 / 3 - 22x^2 * y + 32x^2 - 34x * y^2 / 3 + 74x * y / 3 - 44x / 3 + 32y^3 / 3 - 8y^2 - 8y / 3 + 2
    i == 7 && return 32x^3 / 3 - 34x^2 * y / 3 - 8x^2 - 22x * y^2 + 74x * y / 3 - 8x / 3 - 64y^3 / 3 + 32y^2 - 44y / 3 + 2
    i == 8 && return 53x^2 * y - 16x^2 + 53x * y^2 - 69x * y + 16x + 32y^3 - 48y^2 + 22y - 3
    i == 9 && return -32x^3 / 3 - 130x^2 * y / 3 + 24x^2 - 98x * y^2 / 3 + 46x * y - 40x / 3 - 64y^3 / 3 + 32y^2 - 44y / 3 + 2
    i == 10 && return -27x^2 * y - 27x * y^2 + 27x * y
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::CrouzeixFalk{RefTriangle, 3})
    return [
        Vec((3 / 4, 1 / 4)), Vec((1 / 2, 1 / 2)), Vec((1 / 4, 3 / 4)), # edge 1
        Vec((0.0, 3 / 4)), Vec((0.0, 1 / 2)), Vec((0.0, 1 / 4)),       # edge 2
        Vec((1 / 4, 0.0)), Vec((1 / 2, 0.0)), Vec((3 / 4, 0.0)),       # edge 3
        Vec((1 / 3, 1 / 3)),                                           # centroid
    ]
end

Ferrite.vertexdof_indices(::CrouzeixFalk{RefTriangle, 3}) = ((), (), ())
Ferrite.edgedof_interior_indices(::CrouzeixFalk{RefTriangle, 3}) = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
Ferrite.facedof_interior_indices(::CrouzeixFalk{RefTriangle, 3}) = ((10,),)
Ferrite.volumedof_interior_indices(::CrouzeixFalk{RefTriangle, 3}) = ()
