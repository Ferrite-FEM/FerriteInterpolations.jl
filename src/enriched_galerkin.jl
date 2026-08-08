# Enriched Galerkin element
# (https://defelement.org/elements/enriched-galerkin.html,
# DefElement: elements/enriched-galerkin.def).
#
# Cells/degrees implemented: every cell/degree combination supported by
# `Ferrite.Lagrange` (all seven reference cells), since the implementation
# fully delegates to it: continuous Lagrange of the given degree enriched
# with one discontinuous piecewise-constant function per cell (symfem's basis
# is exactly the Lagrange basis plus the constant 1). Identity mapping.
#
# The Lagrange DOFs keep their usual entity association and are shared
# between neighboring cells exactly as for continuous Lagrange (including
# distribution-time adjustment); the enrichment constant is appended as a
# cell DOF (last index, never shared). The global space therefore contains
# the discontinuous per-cell constants: L2 conformity. Dirichlet facet BCs
# fall back to the Lagrange entity DOFs through Ferrite's defaults.
#
# `reference_coordinates` is not defined: the constant's DOF is an integral
# moment, not a point evaluation, so the element is not nodal.
#
# Note for the pyramid: Ferrite's `Lagrange{RefPyramid}` is a polynomial
# variant while DefElement/symfem use the rational pyramid Lagrange basis,
# so on pyramids this element enriches Ferrite's variant (no symfem
# cross-check there, see the test file).

"""
    EnrichedGalerkin{shape, order}()

Enriched Galerkin element: continuous Lagrange of degree `order` enriched
with one discontinuous piecewise-constant function per cell (appended as the
last, cell-local DOF). Available wherever `Ferrite.Lagrange{shape, order}`
is; discontinuous (L2) as a global space.
"""
struct EnrichedGalerkin{shape, order} <: ScalarInterpolation{shape, order} end

_lagrange(::EnrichedGalerkin{shape, order}) where {shape, order} = Ferrite.Lagrange{shape, order}()

Ferrite.conformity(::EnrichedGalerkin) = Ferrite.L2Conformity()
Ferrite.adjust_dofs_during_distribution(ip::EnrichedGalerkin) = Ferrite.adjust_dofs_during_distribution(_lagrange(ip))

Ferrite.getnbasefunctions(ip::EnrichedGalerkin) = Ferrite.getnbasefunctions(_lagrange(ip)) + 1

function Ferrite.reference_shape_value(ip::EnrichedGalerkin, ξ::Vec{dim, T}, i::Int) where {dim, T}
    lag = _lagrange(ip)
    N = Ferrite.getnbasefunctions(lag)
    i <= N && return Ferrite.reference_shape_value(lag, ξ, i)
    i == N + 1 && return one(T)
    return throw_out_of_range(ip, i)
end

Ferrite.vertexdof_indices(ip::EnrichedGalerkin) = Ferrite.vertexdof_indices(_lagrange(ip))
Ferrite.edgedof_interior_indices(ip::EnrichedGalerkin) = Ferrite.edgedof_interior_indices(_lagrange(ip))
Ferrite.facedof_interior_indices(ip::EnrichedGalerkin) = Ferrite.facedof_interior_indices(_lagrange(ip))
function Ferrite.volumedof_interior_indices(ip::EnrichedGalerkin)
    return (Ferrite.volumedof_interior_indices(_lagrange(ip))..., Ferrite.getnbasefunctions(ip))
end
