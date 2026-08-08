# Radau element (https://defelement.org/elements/radau.html,
# DefElement: elements/radau.def).
#
# Cells/degrees implemented: RefLine degree 2 only. Identity mapping,
# H1-conforming (vertex DOFs shared as usual in 1D). Nodal: Lagrange basis at
# the points {-1, 1, (1 - sqrt(6))/5}, i.e. the endpoints plus the interior
# Radau quadrature point (basis transcribed from symfem
# "Lagrange variant=radau" composed with the [-1, 1] -> [0, 1] map).
#
# NOT implemented, and what would have to change:
#  * Degree 1 coincides with `Ferrite.Lagrange{RefLine, 1}` (use that).
#  * Degrees >= 3: symfem's Radau quadrature (symfem/quadrature.py, radau())
#    raises NotImplementedError for more than 3 points, so there is no
#    reference to verify against; extend symfem first.
#  * RefQuadrilateral/RefHexahedron (any degree): the tensor-product Radau
#    lattice is ASYMMETRIC on each edge (e.g. the interior point at
#    3/5 - sqrt(6)/10 on [0, 1]), so the edge point set does not map onto
#    itself under edge reversal. Ferrite identifies entity DOFs across cells
#    (with at most an order reversal), which would pair DOFs sitting at
#    DIFFERENT physical points on a shared edge with opposite orientation --
#    the resulting "continuity" constraint is wrong, exactly the
#    Fortin-Soulie obstruction. Supporting this conformingly would need
#    orientation-dependent DOF transformations (Kirby style) upstream in
#    Ferrite; a cell-DOF (broken/DG) variant would be possible but would not
#    be the H1 element DefElement describes.

"""
    Radau{shape, order}()

Radau element on the interval, degree 2: Lagrange basis at the endpoints and
the interior Radau point (1 - sqrt(6))/5. H1-conforming. See the source file
for why other cells/degrees are not implemented.
"""
struct Radau{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::Radau) = Ferrite.H1Conformity()
Ferrite.adjust_dofs_during_distribution(::Radau) = false

Ferrite.getnbasefunctions(::Radau{RefLine, 2}) = 3

Ferrite.vertexdof_indices(::Radau{RefLine, 2}) = ((1,), (2,))
Ferrite.edgedof_interior_indices(::Radau{RefLine, 2}) = ((3,),)
Ferrite.facedof_interior_indices(::Radau{RefLine, 2}) = ()
Ferrite.volumedof_interior_indices(::Radau{RefLine, 2}) = ()

function Ferrite.reference_shape_value(ip::Radau{RefLine, 2}, ξ::Vec{1, T}, i::Int) where {T}
    x = ξ[1]
    s = sqrt(6 * one(T))
    i == 1 && return (s / 12 + 1 // 2) * x^2 - x / 2 - s / 12
    i == 2 && return (1 - s / 4) * x^2 + x / 2 - 1 // 2 + s / 4
    i == 3 && return (s / 6 - 3 // 2) * x^2 - s / 6 + 3 // 2
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::Radau{RefLine, 2})
    return [Vec((-1.0,)), Vec((1.0,)), Vec(((1 - sqrt(6)) / 5,))]
end
