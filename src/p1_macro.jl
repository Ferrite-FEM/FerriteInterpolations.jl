# P1 macro element (https://defelement.org/elements/p1-macro.html,
# DefElement: elements/p1-macro.def).
#
# P1 macro element.
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This is a MACRO element: its basis functions are piecewise polynomials on a
# barycentric split of the cell (piecewise linears on the barycentric split plus a bubble). Ferrite interpolations evaluate a
# single smooth expression per basis function via `reference_shape_value`,
# and quadrature rules know nothing about interior sub-cell boundaries, so
# piecewise bases cannot be represented (values would be wrong at/near the
# internal interfaces and quadrature would not resolve the kinks). Upstream
# support for composite/macro reference cells would be needed.

"""
    P1Macro

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/p1_macro.jl`.
"""
struct P1Macro
    function P1Macro()
        return error("P1Macro is not implemented; see the STATUS: BLOCKED notes in src/p1_macro.jl.")
    end
end
