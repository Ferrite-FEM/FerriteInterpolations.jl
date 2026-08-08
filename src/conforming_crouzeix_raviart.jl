# Conforming Crouzeix-Raviart element
# (https://defelement.org/elements/conforming-crouzeix-raviart.html,
# DefElement: elements/conforming-crouzeix-raviart.def).
#
# Cells/degrees implemented: RefTriangle degrees 2-4 (matching DefElement's
# examples; the element exists for any degree k >= 1 and higher degrees follow
# the same pattern, they are just not transcribed here). Degree 1 coincides
# with P1 Lagrange and is therefore not implemented -- use
# `Ferrite.Lagrange{RefTriangle, 1}`. Identity mapping.
#
# Space: P_k enriched with the k-1 degree-(k+1) polynomials
# {x^i y^(k-i) (x+y), i = 1, ..., k-1} (Crouzeix & Raviart, "Conforming and
# nonconforming finite element methods for solving the stationary Stokes
# equations I", RAIRO 7 (1973), DOI 10.1051/m2an/197307R300331). DOFs are
# point evaluations at the k+1 equispaced lattice points of every edge
# (vertices + k-1 interior points) plus interior points.
#
# Conformity: H1. DefElement tags the element `sobolev: L2`, but the shared
# vertex/edge DOFs do determine the traces (hence the name "conforming"): the
# enrichment functions vanish on the two reference legs and restrict to
# degree-k polynomials on the hypotenuse, so the trace of the full space on
# every edge has degree <= k and is fixed by the k+1 shared point values.
# The two-cell test verifies continuity along the whole shared edge.
#
# All DOFs are point evaluations, so `reference_coordinates` is defined.
# DOF order: vertices (1, 2, 3), then k-1 DOFs per edge (ordered from the
# first towards the second vertex of the reference edge), then the interior.
#
# Basis: transcribed from symfem ("conforming Crouzeix-Raviart", degrees
# 2-4); the reference triangles agree pointwise (same coordinates, different
# vertex/edge numbering), see test/test_conforming_crouzeix_raviart.jl for
# the cross-check.

"""
    ConformingCrouzeixRaviart{shape, order}()

Conforming Crouzeix-Raviart element on the triangle, degrees 2-4: P_order
enriched with degree-(order+1) polynomials, with point-evaluation DOFs at the
equispaced edge lattice points and in the interior. H1-conforming (degree 1
coincides with P1 Lagrange and is not implemented).
"""
struct ConformingCrouzeixRaviart{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::ConformingCrouzeixRaviart) = Ferrite.H1Conformity()

Ferrite.adjust_dofs_during_distribution(::ConformingCrouzeixRaviart) = true
Ferrite.adjust_dofs_during_distribution(::ConformingCrouzeixRaviart{<:Any, 2}) = false

Ferrite.vertexdof_indices(::ConformingCrouzeixRaviart{RefTriangle}) = ((1,), (2,), (3,))
Ferrite.volumedof_interior_indices(::ConformingCrouzeixRaviart{RefTriangle}) = ()

##########
# Degree 2
##########

Ferrite.getnbasefunctions(::ConformingCrouzeixRaviart{RefTriangle, 2}) = 7
Ferrite.edgedof_interior_indices(::ConformingCrouzeixRaviart{RefTriangle, 2}) = ((4,), (5,), (6,))
Ferrite.facedof_interior_indices(::ConformingCrouzeixRaviart{RefTriangle, 2}) = ((7,),)

function Ferrite.reference_shape_value(ip::ConformingCrouzeixRaviart{RefTriangle, 2}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return -3 * x^2 * y - 3 * x * y^2 + 2 * x^2 + 3 * x * y - x
    i == 2 && return -3 * x^2 * y - 3 * x * y^2 + 3 * x * y + 2 * y^2 - y
    i == 3 && return -3 * x^2 * y - 3 * x * y^2 + 2 * x^2 + 7 * x * y + 2 * y^2 - 3 * x - 3 * y + 1
    i == 4 && return 12 * x^2 * y + 12 * x * y^2 - 8 * x * y
    i == 5 && return 12 * x^2 * y + 12 * x * y^2 - 16 * x * y - 4 * y^2 + 4 * y
    i == 6 && return 12 * x^2 * y + 12 * x * y^2 - 4 * x^2 - 16 * x * y + 4 * x
    i == 7 && return -27 * x^2 * y - 27 * x * y^2 + 27 * x * y
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::ConformingCrouzeixRaviart{RefTriangle, 2})
    return [
        Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)), # vertices
        Vec((1 / 2, 1 / 2)), Vec((0.0, 1 / 2)), Vec((1 / 2, 0.0)), # edge midpoints
        Vec((1 / 3, 1 / 3)), # interior
    ]
end

##########
# Degree 3
##########

Ferrite.getnbasefunctions(::ConformingCrouzeixRaviart{RefTriangle, 3}) = 12
Ferrite.edgedof_interior_indices(::ConformingCrouzeixRaviart{RefTriangle, 3}) = ((4, 5), (6, 7), (8, 9))
Ferrite.facedof_interior_indices(::ConformingCrouzeixRaviart{RefTriangle, 3}) = ((10, 11, 12),)

function Ferrite.reference_shape_value(ip::ConformingCrouzeixRaviart{RefTriangle, 3}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return (-243 // 20) * x^3 * y + (-243 // 20) * x^2 * y^2 + (9 // 2) * x^3 + (333 // 20) * x^2 * y + (9 // 2) * x * y^2 + (-9 // 2) * x^2 + (-9 // 2) * x * y + x
    i == 2 && return (-243 // 20) * x^2 * y^2 + (-243 // 20) * x * y^3 + (9 // 2) * x^2 * y + (333 // 20) * x * y^2 + (9 // 2) * y^3 + (-9 // 2) * x * y + (-9 // 2) * y^2 + y
    i == 3 && return (243 // 20) * x^3 * y + (243 // 10) * x^2 * y^2 + (243 // 20) * x * y^3 + (-9 // 2) * x^3 + (-333 // 10) * x^2 * y + (-333 // 10) * x * y^2 + (-9 // 2) * y^3 + 9 * x^2 + (513 // 20) * x * y + 9 * y^2 + (-11 // 2) * x + (-11 // 2) * y + 1
    i == 4 && return (243 // 5) * x^3 * y + (729 // 20) * x^2 * y^2 + (-243 // 20) * x * y^3 + (-459 // 10) * x^2 * y + (27 // 20) * x * y^2 + (63 // 10) * x * y
    i == 5 && return (-243 // 20) * x^3 * y + (729 // 20) * x^2 * y^2 + (243 // 5) * x * y^3 + (27 // 20) * x^2 * y + (-459 // 10) * x * y^2 + (63 // 10) * x * y
    i == 6 && return (243 // 20) * x^3 * y + (729 // 10) * x^2 * y^2 + (243 // 4) * x * y^3 + (-351 // 10) * x^2 * y + (-486 // 5) * x * y^2 + (-27 // 2) * y^3 + (549 // 20) * x * y + 18 * y^2 + (-9 // 2) * y
    i == 7 && return (-243 // 5) * x^3 * y + (-2187 // 20) * x^2 * y^2 + (-243 // 4) * x * y^3 + (999 // 10) * x^2 * y + (2511 // 20) * x * y^2 + (27 // 2) * y^3 + (-603 // 10) * x * y + (-45 // 2) * y^2 + 9 * y
    i == 8 && return (-243 // 4) * x^3 * y + (-2187 // 20) * x^2 * y^2 + (-243 // 5) * x * y^3 + (27 // 2) * x^3 + (2511 // 20) * x^2 * y + (999 // 10) * x * y^2 + (-45 // 2) * x^2 + (-603 // 10) * x * y + 9 * x
    i == 9 && return (243 // 4) * x^3 * y + (729 // 10) * x^2 * y^2 + (243 // 20) * x * y^3 + (-27 // 2) * x^3 + (-486 // 5) * x^2 * y + (-351 // 10) * x * y^2 + 18 * x^2 + (549 // 20) * x * y + (-9 // 2) * x
    i == 10 && return (2187 // 20) * x^3 * y + (2187 // 10) * x^2 * y^2 + (2187 // 20) * x * y^3 + (-972 // 5) * x^2 * y + (-972 // 5) * x * y^2 + (1701 // 20) * x * y
    i == 11 && return (-2187 // 20) * x^2 * y^2 + (-2187 // 20) * x * y^3 + (243 // 10) * x^2 * y + (2673 // 20) * x * y^2 + (-243 // 10) * x * y
    i == 12 && return (-2187 // 20) * x^3 * y + (-2187 // 20) * x^2 * y^2 + (2673 // 20) * x^2 * y + (243 // 10) * x * y^2 + (-243 // 10) * x * y
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::ConformingCrouzeixRaviart{RefTriangle, 3})
    return [
        Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)),           # vertices
        Vec((2 / 3, 1 / 3)), Vec((1 / 3, 2 / 3)),                    # edge 1
        Vec((0.0, 2 / 3)), Vec((0.0, 1 / 3)),                        # edge 2
        Vec((1 / 3, 0.0)), Vec((2 / 3, 0.0)),                        # edge 3
        Vec((2 / 9, 2 / 9)), Vec((2 / 9, 5 / 9)), Vec((5 / 9, 2 / 9)), # interior
    ]
end

##########
# Degree 4
##########

Ferrite.getnbasefunctions(::ConformingCrouzeixRaviart{RefTriangle, 4}) = 18
Ferrite.edgedof_interior_indices(::ConformingCrouzeixRaviart{RefTriangle, 4}) = ((4, 5, 6), (7, 8, 9), (10, 11, 12))
Ferrite.facedof_interior_indices(::ConformingCrouzeixRaviart{RefTriangle, 4}) = ((13, 14, 15, 16, 17, 18),)

function Ferrite.reference_shape_value(ip::ConformingCrouzeixRaviart{RefTriangle, 4}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return (-208 // 5) * x^4 * y + (-1264 // 25) * x^3 * y^2 + (-448 // 25) * x^2 * y^3 + (-224 // 25) * x * y^4 + (32 // 3) * x^4 + (5732 // 75) * x^3 * y + (3956 // 75) * x^2 * y^2 + (448 // 25) * x * y^3 - 16 * x^3 + (-3152 // 75) * x^2 * y + (-404 // 25) * x * y^2 + (22 // 3) * x^2 + (36 // 5) * x * y - x
    i == 2 && return (-224 // 25) * x^4 * y + (-448 // 25) * x^3 * y^2 + (-1264 // 25) * x^2 * y^3 + (-208 // 5) * x * y^4 + (448 // 25) * x^3 * y + (3956 // 75) * x^2 * y^2 + (5732 // 75) * x * y^3 + (32 // 3) * y^4 + (-404 // 25) * x^2 * y + (-3152 // 75) * x * y^2 - 16 * y^3 + (36 // 5) * x * y + (22 // 3) * y^2 - y
    i == 3 && return (-208 // 5) * x^4 * y + (-2896 // 25) * x^3 * y^2 + (-2896 // 25) * x^2 * y^3 + (-208 // 5) * x * y^4 + (32 // 3) * x^4 + (3316 // 25) * x^3 * y + (17624 // 75) * x^2 * y^2 + (3316 // 25) * x * y^3 + (32 // 3) * y^4 + (-80 // 3) * x^3 + (-10676 // 75) * x^2 * y + (-10676 // 75) * x * y^2 + (-80 // 3) * y^3 + (70 // 3) * x^2 + (1516 // 25) * x * y + (70 // 3) * y^2 + (-25 // 3) * x + (-25 // 3) * y + 1
    i == 4 && return (896 // 5) * x^4 * y + 128 * x^3 * y^2 + (-128 // 5) * x^2 * y^3 + (128 // 5) * x * y^4 + (-736 // 3) * x^3 * y + (-896 // 15) * x^2 * y^2 + (-416 // 15) * x * y^3 + (464 // 5) * x^2 * y + (272 // 15) * x * y^2 + (-32 // 3) * x * y
    i == 5 && return (-336 // 5) * x^4 * y + (912 // 5) * x^3 * y^2 + (912 // 5) * x^2 * y^3 + (-336 // 5) * x * y^4 + (268 // 5) * x^3 * y + (-1064 // 5) * x^2 * y^2 + (268 // 5) * x * y^3 + (-2 // 5) * x^2 * y + (-2 // 5) * x * y^2 + 2 * x * y
    i == 6 && return (128 // 5) * x^4 * y + (-128 // 5) * x^3 * y^2 + 128 * x^2 * y^3 + (896 // 5) * x * y^4 + (-416 // 15) * x^3 * y + (-896 // 15) * x^2 * y^2 + (-736 // 3) * x * y^3 + (272 // 15) * x^2 * y + (464 // 5) * x * y^2 + (-32 // 3) * x * y
    i == 7 && return (128 // 5) * x^4 * y + 128 * x^3 * y^2 + (1792 // 5) * x^2 * y^3 + 256 * x * y^4 + (-224 // 3) * x^3 * y + (-5408 // 15) * x^2 * y^2 + (-7616 // 15) * x * y^3 + (-128 // 3) * y^4 + (1328 // 15) * x^2 * y + (4208 // 15) * x * y^2 + (224 // 3) * y^3 + (-224 // 5) * x * y + (-112 // 3) * y^2 + (16 // 3) * y
    i == 8 && return (-336 // 5) * x^4 * y + (-2256 // 5) * x^3 * y^2 - 768 * x^2 * y^3 - 384 * x * y^4 + (1076 // 5) * x^3 * y + 980 * x^2 * y^2 + 896 * x * y^3 + 64 * y^4 + (-1214 // 5) * x^2 * y + (-3034 // 5) * x * y^2 - 128 * y^3 + (534 // 5) * x * y + 76 * y^2 - 12 * y
    i == 9 && return (896 // 5) * x^4 * y + (2944 // 5) * x^3 * y^2 + (3328 // 5) * x^2 * y^3 + 256 * x * y^4 + (-7072 // 15) * x^3 * y + (-16352 // 15) * x^2 * y^2 + (-10304 // 15) * x * y^3 + (-128 // 3) * y^4 + 432 * x^2 * y + (1744 // 3) * x * y^2 + 96 * y^3 + (-2336 // 15) * x * y + (-208 // 3) * y^2 + 16 * y
    i == 10 && return 256 * x^4 * y + (3328 // 5) * x^3 * y^2 + (2944 // 5) * x^2 * y^3 + (896 // 5) * x * y^4 + (-128 // 3) * x^4 + (-10304 // 15) * x^3 * y + (-16352 // 15) * x^2 * y^2 + (-7072 // 15) * x * y^3 + 96 * x^3 + (1744 // 3) * x^2 * y + 432 * x * y^2 + (-208 // 3) * x^2 + (-2336 // 15) * x * y + 16 * x
    i == 11 && return -384 * x^4 * y - 768 * x^3 * y^2 + (-2256 // 5) * x^2 * y^3 + (-336 // 5) * x * y^4 + 64 * x^4 + 896 * x^3 * y + 980 * x^2 * y^2 + (1076 // 5) * x * y^3 - 128 * x^3 + (-3034 // 5) * x^2 * y + (-1214 // 5) * x * y^2 + 76 * x^2 + (534 // 5) * x * y - 12 * x
    i == 12 && return 256 * x^4 * y + (1792 // 5) * x^3 * y^2 + 128 * x^2 * y^3 + (128 // 5) * x * y^4 + (-128 // 3) * x^4 + (-7616 // 15) * x^3 * y + (-5408 // 15) * x^2 * y^2 + (-224 // 3) * x * y^3 + (224 // 3) * x^3 + (4208 // 15) * x^2 * y + (1328 // 15) * x * y^2 + (-112 // 3) * x^2 + (-224 // 5) * x * y + (16 // 3) * x
    i == 13 && return -432 * x^4 * y - 1296 * x^3 * y^2 - 1296 * x^2 * y^3 - 432 * x * y^4 + 1044 * x^3 * y + 2088 * x^2 * y^2 + 1044 * x * y^3 - 822 * x^2 * y - 822 * x * y^2 + 210 * x * y
    i == 14 && return (13824 // 25) * x^3 * y^2 + (27648 // 25) * x^2 * y^3 + (13824 // 25) * x * y^4 + (-2304 // 25) * x^3 * y + (-29952 // 25) * x^2 * y^2 + (-27648 // 25) * x * y^3 + (4224 // 25) * x^2 * y + (15744 // 25) * x * y^2 + (-384 // 5) * x * y
    i == 15 && return -432 * x^2 * y^3 - 432 * x * y^4 + 252 * x^2 * y^2 + 684 * x * y^3 - 30 * x^2 * y - 282 * x * y^2 + 30 * x * y
    i == 16 && return (13824 // 25) * x^4 * y + (27648 // 25) * x^3 * y^2 + (13824 // 25) * x^2 * y^3 + (-27648 // 25) * x^3 * y + (-29952 // 25) * x^2 * y^2 + (-2304 // 25) * x * y^3 + (15744 // 25) * x^2 * y + (4224 // 25) * x * y^2 + (-384 // 5) * x * y
    i == 17 && return (-13824 // 25) * x^3 * y^2 + (-13824 // 25) * x^2 * y^3 + (2304 // 25) * x^3 * y + (18432 // 25) * x^2 * y^2 + (2304 // 25) * x * y^3 + (-2688 // 25) * x^2 * y + (-2688 // 25) * x * y^2 + (384 // 25) * x * y
    i == 18 && return -432 * x^4 * y - 432 * x^3 * y^2 + 684 * x^3 * y + 252 * x^2 * y^2 - 282 * x^2 * y - 30 * x * y^2 + 30 * x * y
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::ConformingCrouzeixRaviart{RefTriangle, 4})
    return [
        Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)),              # vertices
        Vec((3 / 4, 1 / 4)), Vec((1 / 2, 1 / 2)), Vec((1 / 4, 3 / 4)),  # edge 1
        Vec((0.0, 3 / 4)), Vec((0.0, 1 / 2)), Vec((0.0, 1 / 4)),        # edge 2
        Vec((1 / 4, 0.0)), Vec((1 / 2, 0.0)), Vec((3 / 4, 0.0)),        # edge 3
        Vec((1 / 6, 1 / 6)), Vec((1 / 6, 5 / 12)), Vec((1 / 6, 2 / 3)), # interior
        Vec((5 / 12, 1 / 6)), Vec((5 / 12, 5 / 12)), Vec((2 / 3, 1 / 6)),
    ]
end
