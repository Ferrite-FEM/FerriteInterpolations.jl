# Johnson-Mercier element (macro stress element) (https://defelement.org/elements/johnson-mercier.html,
# DefElement: elements/johnson-mercier.def).
#
# Johnson-Mercier element (macro stress element).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This is a MACRO element: its basis functions are piecewise polynomials on a
# barycentric split of the cell (symmetric-matrix-valued piecewise polynomials on the Alfeld split (also matrix-valued, cf. the Regge/Arnold-Winther notes)). Ferrite interpolations evaluate a
# single smooth expression per basis function via `reference_shape_value`,
# and quadrature rules know nothing about interior sub-cell boundaries, so
# piecewise bases cannot be represented (values would be wrong at/near the
# internal interfaces and quadrature would not resolve the kinks). Upstream
# support for composite/macro reference cells would be needed.

"""
    JohnsonMercier

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/johnson_mercier.jl`.
"""
struct JohnsonMercier
    function JohnsonMercier()
        return error("JohnsonMercier is not implemented; see the STATUS: BLOCKED notes in src/johnson_mercier.jl.")
    end
end
