# Mardal-Tai-Winther element
# (https://defelement.org/elements/mardal-tai-winther.html,
# DefElement: elements/mardal-tai-winther.def).
#
# Cells/degrees implemented: RefTriangle degree 1 (the element only exists
# for degree 1 on the triangle in DefElement). H(div) element with the
# contravariant Piola mapping, designed as a Stokes velocity element that is
# H(div)-conforming and H1-nonconforming (Mardal, Tai & Winther, "A robust
# finite element method for Darcy-Stokes flow", SINUM 40 (2002)).
#
# The 9 DOFs all sit on the edges: per edge, two normal moments (linear
# Lagrange weights) and one zeroth tangential moment. In Ferrite ordering the
# intra-edge order is (normal at s=0, tangential, normal at s=1) -- a
# PALINDROMIC arrangement, so `adjust_dofs_during_distribution`'s order
# reversal maps the tuple onto itself correctly on reversed edges (the two
# normal moments swap, the tangential moment is its own mirror image), and
# `get_direction` flips the sign of all three (both the outward normal and
# the traversal direction flip for the neighbor).
#
# Basis derivation as for src/bdm.jl: the action of these functionals on
# symfem's "MTW" basis is a signed permutation (verified with sympy). The
# basis functions are cubics. Moment DOFs: not nodal, no
# reference_coordinates.
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# The complete implementation below (kept commented out) is correct on the
# reference cell (it passes the symfem cross-check) and its NORMAL continuity
# works across cells, but the defining TANGENTIAL-moment continuity of the
# MTW method does not survive Ferrite's mapping machinery: the contravariant
# Piola map preserves normal moments exactly, while tangential moments become
# METRIC-WEIGHTED reference moments, detJ * v' * (J J')^-1 * t. Two
# neighboring cells therefore couple different physical functionals whenever
# their Jacobians differ. Verified numerically: the tangential-moment jump of
# a shared-DOF field vanishes to machine precision when the two cells are
# congruent translates (identical J) and is O(1) otherwise. This is exactly
# the class of elements that needs general per-cell DOF transformations
# (R. Kirby, "A general approach to transforming finite elements", SMAI-JCM
# 2018): the physical DOFs are related to the reference DOFs by a
# cell-geometry-dependent matrix, not by a signed permutation. Ferrite's
# orientation machinery (adjust_dofs_during_distribution + get_direction's
# +-1) cannot express that; upstream support for Kirby-style transformations
# would unblock MTW (and Huang-Zhang, Arnold-Winther, Hermite, ...).

"""
    MardalTaiWinther{shape, order}

Placeholder for the Mardal-Tai-Winther element. Not usable: the element
needs cell-dependent DOF transformations that Ferrite does not support, see
the `STATUS: BLOCKED` notes in `src/mardal_tai_winther.jl`.
"""
struct MardalTaiWinther{shape, order}
    function MardalTaiWinther{shape, order}() where {shape, order}
        return error(
            "MardalTaiWinther is not implemented: its tangential-moment DOFs need " *
                "cell-dependent (Kirby-style) DOF transformations that Ferrite does not " *
                "support; see src/mardal_tai_winther.jl."
        )
    end
end

#=
"""
    MardalTaiWinther{shape, order}()

Mardal-Tai-Winther element on the triangle (degree 1): an H(div)-conforming
Stokes velocity element with two normal moments and one tangential moment
per edge (9 DOFs, all shared with the neighboring cells).
"""
struct MardalTaiWinther{shape, order, vdim} <: VectorInterpolation{vdim, shape, order}
    function MardalTaiWinther{shape, order}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}()
    end
    function MardalTaiWinther{shape, order, rdim}() where {rdim, shape <: Ferrite.AbstractRefShape{rdim}, order}
        return new{shape, order, rdim}() # Support construction from `typeof(ip)()`
    end
end

Ferrite.mapping_type(::MardalTaiWinther) = Ferrite.ContravariantPiolaMapping()
Ferrite.conformity(::MardalTaiWinther) = Ferrite.HdivConformity()
Ferrite.adjust_dofs_during_distribution(::MardalTaiWinther) = true

Ferrite.getnbasefunctions(::MardalTaiWinther{RefTriangle, 1}) = 9
Ferrite.edgedof_interior_indices(::MardalTaiWinther{RefTriangle, 1}) = ((1, 2, 3), (4, 5, 6), (7, 8, 9))
Ferrite.facedof_interior_indices(::MardalTaiWinther{RefTriangle, 1}) = ((),)

function Ferrite.get_direction(::MardalTaiWinther{RefTriangle, 1}, shape_nr::Int, cell)
    return Ferrite.get_edge_direction(cell, (shape_nr + 2) ÷ 3)
end

function Ferrite.reference_shape_value(ip::MardalTaiWinther{RefTriangle, 1}, ξ::Vec{2}, i::Int)
    x, y = ξ[1], ξ[2]
    # Edge 1: (normal at s=0, tangential, normal at s=1) along (1,0) -> (0,1)
    i == 1 && return Vec(-30 * x^3 - 144 * x^2 * y - 126 * x * y^2 + 57 * x^2 + 138 * x * y - 23 * x, 90 * x^2 * y + 144 * x * y^2 + 42 * y^3 - 114 * x * y - 69 * y^2 + 25 * y)
    i == 2 && return Vec(-6 * x^3 - 24 * x^2 * y - 18 * x * y^2 + 9 * x^2 + 18 * x * y - 3 * x, 18 * x^2 * y + 24 * x * y^2 + 6 * y^3 - 18 * x * y - 9 * y^2 + 3 * y)
    i == 3 && return Vec(42 * x^3 + 144 * x^2 * y + 90 * x * y^2 - 69 * x^2 - 114 * x * y + 25 * x, -126 * x^2 * y - 144 * x * y^2 - 30 * y^3 + 138 * x * y + 57 * y^2 - 23 * y)
    # Edge 2: (normal at s=0, tangential, normal at s=1) along (0,1) -> (0,0)
    i == 4 && return Vec(48 * x^3 + 120 * x^2 * y + 36 * x * y^2 - 66 * x^2 - 60 * x * y + 16 * x - 6 * y + 2, -144 * x^2 * y - 120 * x * y^2 - 12 * y^3 + 132 * x * y + 30 * y^2 - 14 * y)
    i == 5 && return Vec(12 * x^3 + 24 * x^2 * y - 18 * x^2 - 12 * x * y + 6 * x, -36 * x^2 * y - 24 * x * y^2 + 36 * x * y + 6 * y^2 - 6 * y)
    i == 6 && return Vec(-24 * x^3 - 24 * x^2 * y + 36 * x * y^2 + 24 * x^2 - 24 * x * y + 4 * x + 6 * y - 4, 72 * x^2 * y + 24 * x * y^2 - 12 * y^3 - 48 * x * y + 12 * y^2 - 2 * y)
    # Edge 3: (normal at s=0, tangential, normal at s=1) along (0,0) -> (1,0)
    i == 7 && return Vec(-12 * x^3 + 24 * x^2 * y + 72 * x * y^2 + 12 * x^2 - 48 * x * y - 2 * x, 36 * x^2 * y - 24 * x * y^2 - 24 * y^3 - 24 * x * y + 24 * y^2 + 6 * x + 4 * y - 4)
    i == 8 && return Vec(24 * x^2 * y + 36 * x * y^2 - 6 * x^2 - 36 * x * y + 6 * x, -24 * x * y^2 - 12 * y^3 + 12 * x * y + 18 * y^2 - 6 * y)
    i == 9 && return Vec(-12 * x^3 - 120 * x^2 * y - 144 * x * y^2 + 30 * x^2 + 132 * x * y - 14 * x, 36 * x^2 * y + 120 * x * y^2 + 48 * y^3 - 60 * x * y - 66 * y^2 - 6 * x + 16 * y + 2)
    return throw_out_of_range(ip, i)
end
=#
