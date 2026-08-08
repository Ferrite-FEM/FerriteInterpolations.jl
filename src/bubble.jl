# Bubble element (https://defelement.org/elements/bubble.html,
# DefElement: elements/bubble.def).
#
# Cells/degrees implemented: RefLine degrees 2-3, RefTriangle degrees 3-4,
# RefTetrahedron degree 4 (DefElement's example range; the element exists for
# every degree above the cell's minimum and higher degrees follow the same
# pattern). Identity mapping.
#
# The basis functions are the cell bubble (the product of the barycentric
# coordinates) times a polynomial basis of the complementary degree; all of
# them vanish on the entire cell boundary. The DOFs are point evaluations at
# interior lattice points, so the element is nodal. Since every basis
# function vanishes on the boundary, the global space (extended by zero
# continuity-wise) is continuous: H1Conformity, with all DOFs interior to the
# cell -- nothing is ever shared between cells and no distribution-time
# adjustment is needed.
#
# DOF placement per Ferrite's entity scheme: the cell's own top-dimensional
# entity, i.e. edge-interior for RefLine, face-interior for RefTriangle and
# volume-interior for RefTetrahedron.
#
# Basis: transcribed from symfem ("bubble"); simplex cells agree pointwise
# with symfem's reference cells, the RefLine basis is symfem's composed with
# the [-1, 1] -> [0, 1] coordinate map x = (ξ + 1)/2. See test/test_bubble.jl
# for the cross-check.

"""
    Bubble{shape, order}()

Bubble element: interior-only basis functions vanishing on the whole cell
boundary (the cell bubble times a complementary polynomial basis). RefLine
degrees 2-3, RefTriangle degrees 3-4, RefTetrahedron degree 4.
"""
struct Bubble{shape, order} <: ScalarInterpolation{shape, order} end

Ferrite.conformity(::Bubble) = Ferrite.H1Conformity()
Ferrite.adjust_dofs_during_distribution(::Bubble) = false

###########
# RefLine #
###########

Ferrite.vertexdof_indices(::Bubble{RefLine}) = ((), ())
Ferrite.facedof_interior_indices(::Bubble{RefLine}) = ()
Ferrite.volumedof_interior_indices(::Bubble{RefLine}) = ()

Ferrite.getnbasefunctions(::Bubble{RefLine, 2}) = 1
Ferrite.edgedof_interior_indices(::Bubble{RefLine, 2}) = ((1,),)

function Ferrite.reference_shape_value(ip::Bubble{RefLine, 2}, ξ::Vec{1}, i::Int)
    x = ξ[1]
    i == 1 && return 1 - x^2
    return throw_out_of_range(ip, i)
end

Ferrite.reference_coordinates(::Bubble{RefLine, 2}) = [Vec((0.0,))]

Ferrite.getnbasefunctions(::Bubble{RefLine, 3}) = 2
Ferrite.edgedof_interior_indices(::Bubble{RefLine, 3}) = ((1, 2),)

function Ferrite.reference_shape_value(ip::Bubble{RefLine, 3}, ξ::Vec{1}, i::Int)
    x = ξ[1]
    i == 1 && return (27 // 16) * x^3 + (-9 // 16) * x^2 + (-27 // 16) * x + 9 // 16
    i == 2 && return (-27 // 16) * x^3 + (-9 // 16) * x^2 + (27 // 16) * x + 9 // 16
    return throw_out_of_range(ip, i)
end

Ferrite.reference_coordinates(::Bubble{RefLine, 3}) = [Vec((-1 / 3,)), Vec((1 / 3,))]

###############
# RefTriangle #
###############

Ferrite.vertexdof_indices(::Bubble{RefTriangle}) = ((), (), ())
Ferrite.edgedof_interior_indices(::Bubble{RefTriangle}) = ((), (), ())
Ferrite.volumedof_interior_indices(::Bubble{RefTriangle}) = ()

Ferrite.getnbasefunctions(::Bubble{RefTriangle, 3}) = 1
Ferrite.facedof_interior_indices(::Bubble{RefTriangle, 3}) = ((1,),)

function Ferrite.reference_shape_value(ip::Bubble{RefTriangle, 3}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return -27 * x^2 * y - 27 * x * y^2 + 27 * x * y
    return throw_out_of_range(ip, i)
end

Ferrite.reference_coordinates(::Bubble{RefTriangle, 3}) = [Vec((1 / 3, 1 / 3))]

Ferrite.getnbasefunctions(::Bubble{RefTriangle, 4}) = 3
Ferrite.facedof_interior_indices(::Bubble{RefTriangle, 4}) = ((1, 2, 3),)

function Ferrite.reference_shape_value(ip::Bubble{RefTriangle, 4}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    i == 1 && return 128 * x^3 * y + 256 * x^2 * y^2 + 128 * x * y^3 - 224 * x^2 * y - 224 * x * y^2 + 96 * x * y
    i == 2 && return -128 * x^3 * y - 128 * x^2 * y^2 + 160 * x^2 * y + 32 * x * y^2 - 32 * x * y
    i == 3 && return -128 * x^2 * y^2 - 128 * x * y^3 + 32 * x^2 * y + 160 * x * y^2 - 32 * x * y
    return throw_out_of_range(ip, i)
end

function Ferrite.reference_coordinates(::Bubble{RefTriangle, 4})
    return [Vec((1 / 4, 1 / 4)), Vec((1 / 2, 1 / 4)), Vec((1 / 4, 1 / 2))]
end

##################
# RefTetrahedron #
##################

Ferrite.vertexdof_indices(::Bubble{RefTetrahedron}) = ((), (), (), ())
Ferrite.edgedof_interior_indices(::Bubble{RefTetrahedron}) = ((), (), (), (), (), ())
Ferrite.facedof_interior_indices(::Bubble{RefTetrahedron}) = ((), (), (), ())

Ferrite.getnbasefunctions(::Bubble{RefTetrahedron, 4}) = 1
Ferrite.volumedof_interior_indices(::Bubble{RefTetrahedron, 4}) = (1,)

function Ferrite.reference_shape_value(ip::Bubble{RefTetrahedron, 4}, ξ::Vec{3}, i::Int)
    x, y, z = ξ[1], ξ[2], ξ[3]
    i == 1 && return -256 * x^2 * y * z - 256 * x * y^2 * z - 256 * x * y * z^2 + 256 * x * y * z
    return throw_out_of_range(ip, i)
end

Ferrite.reference_coordinates(::Bubble{RefTetrahedron, 4}) = [Vec((1 / 4, 1 / 4, 1 / 4))]
