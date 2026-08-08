# Hsieh-Clough-Tocher element (C1 macro triangle) (https://defelement.org/elements/hsieh-clough-tocher.html,
# DefElement: elements/hsieh-clough-tocher.def).
#
# Hsieh-Clough-Tocher element (C1 macro triangle). Also covers the reduced variant (reduced-hsieh-clough-tocher.def).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This is a MACRO element: its basis functions are piecewise polynomials on a
# barycentric split of the cell (cubics on the three subtriangles of the barycentric split, C1 across the internal interfaces). Ferrite interpolations evaluate a
# single smooth expression per basis function via `reference_shape_value`,
# and quadrature rules know nothing about interior sub-cell boundaries, so
# piecewise bases cannot be represented (values would be wrong at/near the
# internal interfaces and quadrature would not resolve the kinks). Upstream
# support for composite/macro reference cells would be needed.

"""
    HsiehCloughTocher

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/hsieh_clough_tocher.jl`.
"""
struct HsiehCloughTocher
    function HsiehCloughTocher()
        return error("HsiehCloughTocher is not implemented; see the STATUS: BLOCKED notes in src/hsieh_clough_tocher.jl.")
    end
end
