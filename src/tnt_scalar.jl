# Scalar tiniest-tensor (TNT) element on quadrilaterals/hexahedra (https://defelement.org/elements/tnt.html,
# DefElement: elements/tnt.def).
#
# Scalar tiniest-tensor (TNT) element on quadrilaterals/hexahedra. (The H(div)/H(curl) TNT variants ARE implemented, src/tnt_div.jl and src/tnt_curl.jl.)
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# The scalar tiniest-tensor element (tnt.def; degree >= 2, symfem raises for degree 1) has edge DOFs that are integral moments against DERIVATIVES of an edge Lagrange basis. Under edge reversal these weights transform with SIGNS (odd/even
# parity), not by a permutation. Ferrite's scalar-interpolation path has only
# the order-reversal of `adjust_dofs_during_distribution`; per-DOF sign flips
# (`get_direction`) exist only for Piola-mapped vector elements. Supporting
# this needs sign-aware distribution for scalar elements (or full Kirby-style
# transformations) upstream.

"""
    TNTScalar

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/tnt_scalar.jl`.
"""
struct TNTScalar
    function TNTScalar()
        return error("TNTScalar is not implemented; see the STATUS: BLOCKED notes in src/tnt_scalar.jl.")
    end
end
