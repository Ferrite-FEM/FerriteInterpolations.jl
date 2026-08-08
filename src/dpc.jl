# dPc element ("discontinuous polynomial cubical",
# https://defelement.org/elements/dpc.html, DefElement: elements/dpc.def).
#
# Cells/degrees implemented: RefQuadrilateral degrees 1-3, RefHexahedron
# degrees 1-2 (higher degrees follow the same pattern). The element also
# exists on the interval, where it coincides with
# `Ferrite.DiscontinuousLagrange{RefLine, k}` and is therefore not
# implemented here.
#
# Space: the TOTAL-degree polynomials P_k on the tensor-product cell (not
# Q_k). DOFs are point evaluations at the simplex-equispaced lattice in the
# cell's first-vertex corner (the far corner carries no points), following
# symfem/basix's `simplex_equispaced` variant; all DOFs are cell DOFs
# (DiscontinuousLagrange pattern), L2 conformity.
#
# Mapping deviation: DefElement specifies `mapping: L2 Piola` (values scale
# by 1/det J). Ferrite has no such mapping, so the identity mapping is used,
# as for Ferrite's own DiscontinuousLagrange. On affinely-mapped cells
# (parallelograms/parallelepipeds) the spanned space is identical (the
# 1/det J factor is a per-cell constant); on non-affinely-mapped cells the
# two variants span different physical spaces.
#
# Unlike DiscontinuousLagrange, no `dirichlet_*dof_indices` are defined: the
# asymmetric lattice leaves most facets with an incomplete trace point set
# (e.g. the far edges of the quadrilateral carry a single point), so facet
# Dirichlet conditions cannot be represented anyway.
#
# All DOFs are point evaluations, so `reference_coordinates` is defined.
#
# Basis: transcribed from symfem ("dPc") composed with the [-1, 1]^d ->
# [0, 1]^d coordinate map x = (ξ + 1)/2, see test/test_dpc.jl for the
# cross-check.

"""
    DPC{shape, order}()

dPc element ("discontinuous polynomial cubical"): the total-degree polynomial
space P_order on tensor-product cells, with point-evaluation cell DOFs at a
corner-based simplex lattice. Discontinuous (L2). RefQuadrilateral degrees
1-3, RefHexahedron degrees 1-2. The vector-valued "vector dPc" element is
`DPC{shape, order}()^vdim` (`Ferrite.VectorizedInterpolation`).
"""
struct DPC{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::DPC) = Ferrite.L2Conformity()
Ferrite.adjust_dofs_during_distribution(::DPC) = false

# All DOFs in the cell interior (not shared between cells), like
# DiscontinuousLagrange.
Ferrite.volumedof_interior_indices(ip::DPC) = ntuple(i -> i, Ferrite.getnbasefunctions(ip))

####################
# RefQuadrilateral #
####################

Ferrite.getnbasefunctions(::DPC{RefQuadrilateral, 1}) = 3

function Ferrite.reference_shape_value(ip::DPC{RefQuadrilateral, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return (-1 // 2) * x + (-1 // 2) * y
    i == 2 && return (1 // 2) * x + 1 // 2
    i == 3 && return (1 // 2) * y + 1 // 2
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::DPC{RefQuadrilateral, 1})
    return [Vec((-1.0, -1.0)), Vec((1.0, -1.0)), Vec((-1.0, 1.0))]
end

Ferrite.getnbasefunctions(::DPC{RefQuadrilateral, 2}) = 6

function Ferrite.reference_shape_value(ip::DPC{RefQuadrilateral, 2}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return (1 // 2) * x^2 + x * y + (1 // 2) * y^2 + (1 // 2) * x + (1 // 2) * y
    i == 2 && return -x^2 - x * y - x - y
    i == 3 && return (1 // 2) * x^2 + (1 // 2) * x
    i == 4 && return -x * y - y^2 - x - y
    i == 5 && return x * y + x + y + 1
    i == 6 && return (1 // 2) * y^2 + (1 // 2) * y
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::DPC{RefQuadrilateral, 2})
    return [
        Vec((-1.0, -1.0)), Vec((0.0, -1.0)), Vec((1.0, -1.0)),
        Vec((-1.0, 0.0)), Vec((0.0, 0.0)), Vec((-1.0, 1.0)),
    ]
end

Ferrite.getnbasefunctions(::DPC{RefQuadrilateral, 3}) = 10

function Ferrite.reference_shape_value(ip::DPC{RefQuadrilateral, 3}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return (-9 // 16) * x^3 + (-27 // 16) * x^2 * y + (-27 // 16) * x * y^2 + (-9 // 16) * y^3 + (-9 // 8) * x^2 + (-9 // 4) * x * y + (-9 // 8) * y^2 + (-1 // 2) * x + (-1 // 2) * y
    i == 2 && return (27 // 16) * x^3 + (27 // 8) * x^2 * y + (27 // 16) * x * y^2 + (45 // 16) * x^2 + (9 // 2) * x * y + (27 // 16) * y^2 + (9 // 8) * x + (9 // 8) * y
    i == 3 && return (-27 // 16) * x^3 + (-27 // 16) * x^2 * y + (-9 // 4) * x^2 + (-9 // 4) * x * y + (-9 // 16) * x + (-9 // 16) * y
    i == 4 && return (9 // 16) * x^3 + (9 // 16) * x^2 + (-1 // 16) * x - 1 // 16
    i == 5 && return (27 // 16) * x^2 * y + (27 // 8) * x * y^2 + (27 // 16) * y^3 + (27 // 16) * x^2 + (9 // 2) * x * y + (45 // 16) * y^2 + (9 // 8) * x + (9 // 8) * y
    i == 6 && return (-27 // 8) * x^2 * y + (-27 // 8) * x * y^2 + (-27 // 8) * x^2 + (-27 // 4) * x * y + (-27 // 8) * y^2 + (-27 // 8) * x + (-27 // 8) * y
    i == 7 && return (27 // 16) * x^2 * y + (27 // 16) * x^2 + (9 // 4) * x * y + (9 // 4) * x + (9 // 16) * y + 9 // 16
    i == 8 && return (-27 // 16) * x * y^2 + (-27 // 16) * y^3 + (-9 // 4) * x * y + (-9 // 4) * y^2 + (-9 // 16) * x + (-9 // 16) * y
    i == 9 && return (27 // 16) * x * y^2 + (9 // 4) * x * y + (27 // 16) * y^2 + (9 // 16) * x + (9 // 4) * y + 9 // 16
    i == 10 && return (9 // 16) * y^3 + (9 // 16) * y^2 + (-1 // 16) * y - 1 // 16
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::DPC{RefQuadrilateral, 3})
    return [
        Vec((-1.0, -1.0)), Vec((-1 / 3, -1.0)), Vec((1 / 3, -1.0)), Vec((1.0, -1.0)),
        Vec((-1.0, -1 / 3)), Vec((-1 / 3, -1 / 3)), Vec((1 / 3, -1 / 3)),
        Vec((-1.0, 1 / 3)), Vec((-1 / 3, 1 / 3)),
        Vec((-1.0, 1.0)),
    ]
end

#################
# RefHexahedron #
#################

Ferrite.getnbasefunctions(::DPC{RefHexahedron, 1}) = 4

function Ferrite.reference_shape_value(ip::DPC{RefHexahedron, 1}, ξ::Vec{3}, i::Int)
    x, y, z = ξ[1], ξ[2], ξ[3]
    i == 1 && return (-1 // 2) * x + (-1 // 2) * y + (-1 // 2) * z - 1 // 2
    i == 2 && return (1 // 2) * x + 1 // 2
    i == 3 && return (1 // 2) * y + 1 // 2
    i == 4 && return (1 // 2) * z + 1 // 2
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::DPC{RefHexahedron, 1})
    return [
        Vec((-1.0, -1.0, -1.0)), Vec((1.0, -1.0, -1.0)),
        Vec((-1.0, 1.0, -1.0)), Vec((-1.0, -1.0, 1.0)),
    ]
end

Ferrite.getnbasefunctions(::DPC{RefHexahedron, 2}) = 10

function Ferrite.reference_shape_value(ip::DPC{RefHexahedron, 2}, ξ::Vec{3}, i::Int)
    x, y, z = ξ[1], ξ[2], ξ[3]
    i == 1 && return (1 // 2) * x^2 + x * y + x * z + (1 // 2) * y^2 + y * z + (1 // 2) * z^2 + (3 // 2) * x + (3 // 2) * y + (3 // 2) * z + 1
    i == 2 && return -x^2 - x * y - x * z - 2 * x - y - z - 1
    i == 3 && return (1 // 2) * x^2 + (1 // 2) * x
    i == 4 && return -x * y - y^2 - y * z - x - 2 * y - z - 1
    i == 5 && return x * y + x + y + 1
    i == 6 && return (1 // 2) * y^2 + (1 // 2) * y
    i == 7 && return -x * z - y * z - z^2 - x - y - 2 * z - 1
    i == 8 && return x * z + x + z + 1
    i == 9 && return y * z + y + z + 1
    i == 10 && return (1 // 2) * z^2 + (1 // 2) * z
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::DPC{RefHexahedron, 2})
    return [
        Vec((-1.0, -1.0, -1.0)), Vec((0.0, -1.0, -1.0)), Vec((1.0, -1.0, -1.0)),
        Vec((-1.0, 0.0, -1.0)), Vec((0.0, 0.0, -1.0)), Vec((-1.0, 1.0, -1.0)),
        Vec((-1.0, -1.0, 0.0)), Vec((0.0, -1.0, 0.0)), Vec((-1.0, 0.0, 0.0)),
        Vec((-1.0, -1.0, 1.0)),
    ]
end
