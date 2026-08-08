# P1-iso-P2 element (https://defelement.org/elements/p1-iso-p2.html,
# DefElement: elements/p1-iso-p2.def).
#
# P1-iso-P2 element.
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This is a MACRO element: its basis functions are piecewise polynomials on a
# barycentric split of the cell (piecewise linears on the uniform refinement of the cell). Ferrite interpolations evaluate a
# single smooth expression per basis function via `reference_shape_value`,
# and quadrature rules know nothing about interior sub-cell boundaries, so
# piecewise bases cannot be represented (values would be wrong at/near the
# internal interfaces and quadrature would not resolve the kinks). Upstream
# support for composite/macro reference cells would be needed.

"""
    P1IsoP2

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/p1_iso_p2.jl`.
"""
struct P1IsoP2
    function P1IsoP2()
        return error("P1IsoP2 is not implemented; see the STATUS: BLOCKED notes in src/p1_iso_p2.jl.")
    end
end
