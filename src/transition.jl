# Transition element (https://defelement.org/elements/transition.html,
# DefElement: elements/transition.def).
#
# Cells/degrees implemented: RefTriangle with interior order 2 and per-edge
# orders 1 or 2, given as a tuple type parameter in FERRITE edge order, e.g.
# `Transition{RefTriangle, 2, (1, 2, 1)}()` has a quadratic (midpoint) DOF
# only on edge 2. All seven nontrivial combinations are provided; the
# all-ones combination coincides with `Ferrite.Lagrange{RefTriangle, 1}` and
# is not implemented (while (2, 2, 2) coincides with quadratic Lagrange but
# is kept for completeness of the family). Higher orders and the tetrahedron
# follow the same pattern (symfem: `transition` with `edge_orders` /
# `face_orders`) but are not transcribed.
#
# The element enables CONFORMING meshes with mixed polynomial degree: on an
# order-1 edge the trace is linear (determined by the vertex DOFs), so the
# cell can neighbor a P1 cell; on an order-2 edge the trace is quadratic and
# couples with P2/transition neighbors through the shared midpoint DOF.
# Nodal, identity mapping, H1 conformity, at most one DOF per edge so no
# distribution-time adjustment is needed.
#
# Basis: transcribed from symfem ("transition", order 2, with the edge
# orders permuted to symfem's edge numbering: symfem edge_orders =
# [ferrite_e3, ferrite_e2, ferrite_e1]); see test/test_transition.jl.

"""
    Transition{shape, order, edge_orders}()

Transition element on the triangle with interior order 2 and per-edge orders
1 or 2 (`edge_orders` is a tuple in Ferrite edge order): quadratic Lagrange
with the quadratic edge DOF omitted on order-1 edges, whose traces are then
linear. Enables conforming degree transitions in mixed P1/P2 meshes.
"""
struct Transition{shape, order, EO} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::Transition) = Ferrite.H1Conformity()
Ferrite.adjust_dofs_during_distribution(::Transition) = false

Ferrite.vertexdof_indices(::Transition{RefTriangle}) = ((1,), (2,), (3,))
Ferrite.facedof_interior_indices(::Transition{RefTriangle, 2}) = ((),)
Ferrite.volumedof_interior_indices(::Transition{RefTriangle}) = ()

# Edge orders (1, 1, 2); Ferrite DOF -> 0-based symfem DOF perm (with symfem
# edge_orders [2, 1, 1]): [1, 2, 0, 3]
Ferrite.getnbasefunctions(::Transition{RefTriangle, 2, (1, 1, 2)}) = 4
Ferrite.edgedof_interior_indices(::Transition{RefTriangle, 2, (1, 1, 2)}) = ((), (), (4,))
Ferrite.reference_coordinates(::Transition{RefTriangle, 2, (1, 1, 2)}) = [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)), Vec((1 / 2, 0.0))]
function Ferrite.reference_shape_value(ip::Transition{RefTriangle, 2, (1, 1, 2)}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return 2 * x^2 + 2 * x * y - x
    i == 2 && return y
    i == 3 && return 2 * x^2 + 2 * x * y - 3 * x - y + 1
    i == 4 && return -4 * x^2 - 4 * x * y + 4 * x
    return throw_out_of_range(ip, i)
end

# Edge orders (1, 2, 1); Ferrite DOF -> 0-based symfem DOF perm (with symfem
# edge_orders [1, 2, 1]): [1, 2, 0, 3]
Ferrite.getnbasefunctions(::Transition{RefTriangle, 2, (1, 2, 1)}) = 4
Ferrite.edgedof_interior_indices(::Transition{RefTriangle, 2, (1, 2, 1)}) = ((), (4,), ())
Ferrite.reference_coordinates(::Transition{RefTriangle, 2, (1, 2, 1)}) = [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)), Vec((0.0, 1 / 2))]
function Ferrite.reference_shape_value(ip::Transition{RefTriangle, 2, (1, 2, 1)}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return x
    i == 2 && return 2 * x * y + 2 * y^2 - y
    i == 3 && return 2 * x * y + 2 * y^2 - x - 3 * y + 1
    i == 4 && return -4 * x * y - 4 * y^2 + 4 * y
    return throw_out_of_range(ip, i)
end

# Edge orders (1, 2, 2); Ferrite DOF -> 0-based symfem DOF perm (with symfem
# edge_orders [2, 2, 1]): [1, 2, 0, 4, 3]
Ferrite.getnbasefunctions(::Transition{RefTriangle, 2, (1, 2, 2)}) = 5
Ferrite.edgedof_interior_indices(::Transition{RefTriangle, 2, (1, 2, 2)}) = ((), (4,), (5,))
Ferrite.reference_coordinates(::Transition{RefTriangle, 2, (1, 2, 2)}) = [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)), Vec((0.0, 1 / 2)), Vec((1 / 2, 0.0))]
function Ferrite.reference_shape_value(ip::Transition{RefTriangle, 2, (1, 2, 2)}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return 2 * x^2 + 2 * x * y - x
    i == 2 && return 2 * x * y + 2 * y^2 - y
    i == 3 && return 2 * x^2 + 4 * x * y + 2 * y^2 - 3 * x - 3 * y + 1
    i == 4 && return -4 * x * y - 4 * y^2 + 4 * y
    i == 5 && return -4 * x^2 - 4 * x * y + 4 * x
    return throw_out_of_range(ip, i)
end

# Edge orders (2, 1, 1); Ferrite DOF -> 0-based symfem DOF perm (with symfem
# edge_orders [1, 1, 2]): [1, 2, 0, 3]
Ferrite.getnbasefunctions(::Transition{RefTriangle, 2, (2, 1, 1)}) = 4
Ferrite.edgedof_interior_indices(::Transition{RefTriangle, 2, (2, 1, 1)}) = ((4,), (), ())
Ferrite.reference_coordinates(::Transition{RefTriangle, 2, (2, 1, 1)}) = [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)), Vec((1 / 2, 1 / 2))]
function Ferrite.reference_shape_value(ip::Transition{RefTriangle, 2, (2, 1, 1)}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return -2 * x * y + x
    i == 2 && return -2 * x * y + y
    i == 3 && return -x - y + 1
    i == 4 && return 4 * x * y
    return throw_out_of_range(ip, i)
end

# Edge orders (2, 1, 2); Ferrite DOF -> 0-based symfem DOF perm (with symfem
# edge_orders [2, 1, 2]): [1, 2, 0, 4, 3]
Ferrite.getnbasefunctions(::Transition{RefTriangle, 2, (2, 1, 2)}) = 5
Ferrite.edgedof_interior_indices(::Transition{RefTriangle, 2, (2, 1, 2)}) = ((4,), (), (5,))
Ferrite.reference_coordinates(::Transition{RefTriangle, 2, (2, 1, 2)}) = [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)), Vec((1 / 2, 1 / 2)), Vec((1 / 2, 0.0))]
function Ferrite.reference_shape_value(ip::Transition{RefTriangle, 2, (2, 1, 2)}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return 2 * x^2 - x
    i == 2 && return -2 * x * y + y
    i == 3 && return 2 * x^2 + 2 * x * y - 3 * x - y + 1
    i == 4 && return 4 * x * y
    i == 5 && return -4 * x^2 - 4 * x * y + 4 * x
    return throw_out_of_range(ip, i)
end

# Edge orders (2, 2, 1); Ferrite DOF -> 0-based symfem DOF perm (with symfem
# edge_orders [1, 2, 2]): [1, 2, 0, 4, 3]
Ferrite.getnbasefunctions(::Transition{RefTriangle, 2, (2, 2, 1)}) = 5
Ferrite.edgedof_interior_indices(::Transition{RefTriangle, 2, (2, 2, 1)}) = ((4,), (5,), ())
Ferrite.reference_coordinates(::Transition{RefTriangle, 2, (2, 2, 1)}) = [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)), Vec((1 / 2, 1 / 2)), Vec((0.0, 1 / 2))]
function Ferrite.reference_shape_value(ip::Transition{RefTriangle, 2, (2, 2, 1)}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return -2 * x * y + x
    i == 2 && return 2 * y^2 - y
    i == 3 && return 2 * x * y + 2 * y^2 - x - 3 * y + 1
    i == 4 && return 4 * x * y
    i == 5 && return -4 * x * y - 4 * y^2 + 4 * y
    return throw_out_of_range(ip, i)
end

# Edge orders (2, 2, 2); Ferrite DOF -> 0-based symfem DOF perm (with symfem
# edge_orders [2, 2, 2]): [1, 2, 0, 5, 4, 3]
Ferrite.getnbasefunctions(::Transition{RefTriangle, 2, (2, 2, 2)}) = 6
Ferrite.edgedof_interior_indices(::Transition{RefTriangle, 2, (2, 2, 2)}) = ((4,), (5,), (6,))
Ferrite.reference_coordinates(::Transition{RefTriangle, 2, (2, 2, 2)}) = [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0)), Vec((1 / 2, 1 / 2)), Vec((0.0, 1 / 2)), Vec((1 / 2, 0.0))]
function Ferrite.reference_shape_value(ip::Transition{RefTriangle, 2, (2, 2, 2)}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return 2 * x^2 - x
    i == 2 && return 2 * y^2 - y
    i == 3 && return 2 * x^2 + 4 * x * y + 2 * y^2 - 3 * x - 3 * y + 1
    i == 4 && return 4 * x * y
    i == 5 && return -4 * x * y - 4 * y^2 + 4 * y
    i == 6 && return -4 * x^2 - 4 * x * y + 4 * x
    return throw_out_of_range(ip, i)
end
