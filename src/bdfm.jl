# Brezzi-Douglas-Fortin-Marini element
# (https://defelement.org/elements/brezzi-douglas-fortin-marini.html,
# DefElement: elements/brezzi-douglas-fortin-marini.def).
#
# Cells/degrees implemented: RefTriangle degrees 1-2 (higher degrees follow
# the same pattern; the quadrilateral/tetrahedron variants would additionally
# need the same face-orientation caveats as the other H(div) elements,
# cf. PLAN.md).
#
# Space: {p in P_(k+1)^2 : p.n in P_k on every edge} (dimension 9 for k = 1,
# 17 for k = 2). H(div) element with the contravariant Piola mapping. DOFs:
# k+1 normal moments per edge (degree-k Lagrange weights, ordered along the
# FERRITE edge direction, outward reference normal) plus interior moments
# (symfem's; 3 for k = 1, 8 for k = 2).
#
# Orientation handling and basis derivation exactly as for src/bdm.jl: the
# action of the Ferrite-convention functionals on symfem's "BDFM" basis is a
# signed permutation (verified with sympy) and the basis below is that signed
# permutation of symfem's. Moment DOFs: not nodal, no reference_coordinates.

"""
    BDFM{shape, order}()

Brezzi-Douglas-Fortin-Marini H(div) element on the triangle, degrees 1-2:
vector polynomials of degree order+1 whose normal trace has degree order,
with order+1 normal-moment DOFs per edge plus interior moments.
"""
struct BDFM{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function BDFM{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function BDFM{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::BDFM) = Ferrite.ContravariantPiolaMapping()
Ferrite.conformity(::BDFM) = Ferrite.HdivConformity()
Ferrite.adjust_dofs_during_distribution(::BDFM) = true

##########
# Degree 1
##########

Ferrite.getnbasefunctions(::BDFM{RefTriangle, 1}) = 9
Ferrite.edgedof_interior_indices(::BDFM{RefTriangle, 1}) = ((1, 2), (3, 4), (5, 6))
Ferrite.facedof_interior_indices(::BDFM{RefTriangle, 1}) = ((7, 8, 9),)

function Ferrite.get_direction(::BDFM{RefTriangle, 1}, shape_nr::Int, cell)
    shape_nr > 6 && return 1
    return Ferrite.get_edge_direction(cell, (shape_nr + 1) ÷ 2)
end

function Ferrite.reference_shape_value(ip::BDFM{RefTriangle, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    # Edge 1: dofs at s = 0, 1 along (1,0) -> (0,1)
    i == 1 && return Vec(13 * x^2 + 10 * x * y - 9 * x, -2 * x * y - 5 * y^2 + 3 * y)
    i == 2 && return Vec(-5 * x^2 - 2 * x * y + 3 * x, 10 * x * y + 13 * y^2 - 9 * y)
    # Edge 2: dofs at s = 0, 1 along (0,1) -> (0,0)
    i == 3 && return Vec(5 * x^2 + 18 * x * y - 7 * x - 6 * y + 2, -10 * x * y + 3 * y^2 + y)
    i == 4 && return Vec(-13 * x^2 - 18 * x * y + 17 * x + 6 * y - 4, 2 * x * y - 3 * y^2 + y)
    # Edge 3: dofs at s = 0, 1 along (0,0) -> (1,0)
    i == 5 && return Vec(-3 * x^2 + 2 * x * y + x, -18 * x * y - 13 * y^2 + 6 * x + 17 * y - 4)
    i == 6 && return Vec(3 * x^2 - 10 * x * y + x, 18 * x * y + 5 * y^2 - 6 * x - 7 * y + 2)
    # Interior moments
    i == 7 && return Vec(-36 * x^2 - 48 * x * y + 36 * x, 24 * x * y + 12 * y^2 - 12 * y)
    i == 8 && return Vec(12 * x^2 + 24 * x * y - 12 * x, -48 * x * y - 36 * y^2 + 36 * y)
    i == 9 && return Vec(-12 * x^2 - 48 * x * y + 12 * x, 48 * x * y + 12 * y^2 - 12 * y)
    return throw_out_of_range(ip, i)
end

##########
# Degree 2
##########

Ferrite.getnbasefunctions(::BDFM{RefTriangle, 2}) = 17
Ferrite.edgedof_interior_indices(::BDFM{RefTriangle, 2}) = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
Ferrite.facedof_interior_indices(::BDFM{RefTriangle, 2}) = ((10, 11, 12, 13, 14, 15, 16, 17),)

function Ferrite.get_direction(::BDFM{RefTriangle, 2}, shape_nr::Int, cell)
    shape_nr > 9 && return 1
    return Ferrite.get_edge_direction(cell, (shape_nr + 2) ÷ 3)
end

function Ferrite.reference_shape_value(ip::BDFM{RefTriangle, 2}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    # Edge 1: dofs at s = 0, 1/2, 1 along (1,0) -> (0,1)
    i == 1 && return Vec(73 * x^3 + 42 * x^2 * y - 21 * x * y^2 - 80 * x^2 + 16 * x, -39 * x^2 * y - 42 * x * y^2 + 7 * y^3 + 40 * x * y - 4 * y)
    i == 2 && return Vec((5 // 2) * x^3 + 45 * x^2 * y + (75 // 2) * x * y^2 - 10 * x^2 - 40 * x * y + 6 * x, (75 // 2) * x^2 * y + 45 * x * y^2 + (5 // 2) * y^3 - 40 * x * y - 10 * y^2 + 6 * y)
    i == 3 && return Vec(7 * x^3 - 42 * x^2 * y - 39 * x * y^2 + 40 * x * y - 4 * x, -21 * x^2 * y + 42 * x * y^2 + 73 * y^3 - 80 * y^2 + 16 * y)
    # Edge 2: dofs at s = 0, 1/2, 1 along (0,1) -> (0,0)
    i == 4 && return Vec(7 * x^3 + 84 * x^2 * y + 150 * x * y^2 - 21 * x^2 - 128 * x * y - 30 * y^2 + 17 * x + 24 * y - 3, -21 * x^2 * y - 84 * x * y^2 + 10 * y^3 + 42 * x * y + 4 * y^2 - 5 * y)
    i == 5 && return Vec((5 // 2) * x^3 - 75 * x^2 * y - 75 * x * y^2 + (5 // 2) * x^2 + 90 * x * y + 15 * y^2 + (-13 // 2) * x - 15 * y + 3 // 2, (75 // 2) * x^2 * y + 30 * x * y^2 - 5 * y^3 - 35 * x * y + (7 // 2) * y)
    i == 6 && return Vec(73 * x^3 + 216 * x^2 * y + 150 * x * y^2 - 139 * x^2 - 232 * x * y - 30 * y^2 + 75 * x + 36 * y - 9, -39 * x^2 * y - 36 * x * y^2 + 10 * y^3 + 38 * x * y - 4 * y^2 - 3 * y)
    # Edge 3: dofs at s = 0, 1/2, 1 along (0,0) -> (1,0)
    i == 7 && return Vec(10 * x^3 - 36 * x^2 * y - 39 * x * y^2 - 4 * x^2 + 38 * x * y - 3 * x, 150 * x^2 * y + 216 * x * y^2 + 73 * y^3 - 30 * x^2 - 232 * x * y - 139 * y^2 + 36 * x + 75 * y - 9)
    i == 8 && return Vec(-5 * x^3 + 30 * x^2 * y + (75 // 2) * x * y^2 - 35 * x * y + (7 // 2) * x, -75 * x^2 * y - 75 * x * y^2 + (5 // 2) * y^3 + 15 * x^2 + 90 * x * y + (5 // 2) * y^2 - 15 * x + (-13 // 2) * y + 3 // 2)
    i == 9 && return Vec(10 * x^3 - 84 * x^2 * y - 21 * x * y^2 + 4 * x^2 + 42 * x * y - 5 * x, 150 * x^2 * y + 84 * x * y^2 + 7 * y^3 - 30 * x^2 - 128 * x * y - 21 * y^2 + 24 * x + 17 * y - 3)
    # Interior moments
    i == 10 && return Vec(80 * x^3 + 300 * x^2 * y + 240 * x * y^2 - 160 * x^2 - 300 * x * y + 80 * x, -60 * x^2 * y - 120 * x * y^2 - 40 * y^3 + 80 * x * y + 60 * y^2 - 20 * y)
    i == 11 && return Vec(-80 * x^3 + 120 * x * y^2 + 80 * x^2 - 100 * x * y, 60 * x^2 * y - 20 * y^3 - 40 * x * y + 20 * y^2)
    i == 12 && return Vec(-40 * x^3 - 120 * x^2 * y - 60 * x * y^2 + 60 * x^2 + 80 * x * y - 20 * x, 240 * x^2 * y + 300 * x * y^2 + 80 * y^3 - 300 * x * y - 160 * y^2 + 80 * y)
    i == 13 && return Vec(-20 * x^3 + 60 * x * y^2 + 20 * x^2 - 40 * x * y, 120 * x^2 * y - 80 * y^3 - 100 * x * y + 80 * y^2)
    i == 14 && return Vec(-40 * x^3 - 240 * x^2 * y - 120 * x * y^2 + 60 * x^2 + 140 * x * y - 20 * x, 240 * x^2 * y + 180 * x * y^2 + 20 * y^3 - 180 * x * y - 40 * y^2 + 20 * y)
    i == 15 && return Vec(-20 * x^3 - 180 * x^2 * y - 240 * x * y^2 + 40 * x^2 + 180 * x * y - 20 * x, 120 * x^2 * y + 240 * x * y^2 + 40 * y^3 - 140 * x * y - 60 * y^2 + 20 * y)
    i == 16 && return Vec(-30 * x^2 * y - 15 * x * y^2 - 5 * x^2 + 20 * x * y + 5 * x, -30 * x * y^2 - 15 * y^3 + 10 * x * y + 20 * y^2 - 5 * y)
    i == 17 && return Vec(-15 * x^3 - 30 * x^2 * y + 20 * x^2 + 10 * x * y - 5 * x, -15 * x^2 * y - 30 * x * y^2 + 20 * x * y - 5 * y^2 + 5 * y)
    return throw_out_of_range(ip, i)
end
