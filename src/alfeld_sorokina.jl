# Alfeld-Sorokina element (C0 vector macro element with continuous divergence) (https://defelement.org/elements/alfeld-sorokina.html,
# DefElement: elements/alfeld-sorokina.def).
#
# Alfeld-Sorokina element (C0 vector macro element with continuous divergence).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This is a MACRO element: its basis functions are piecewise polynomials on a
# barycentric split of the cell (the Alfeld/barycentric split, combined with derivative DOFs). Ferrite interpolations evaluate a
# single smooth expression per basis function via `reference_shape_value`,
# and quadrature rules know nothing about interior sub-cell boundaries, so
# piecewise bases cannot be represented (values would be wrong at/near the
# internal interfaces and quadrature would not resolve the kinks). Upstream
# support for composite/macro reference cells would be needed.

"""
    AlfeldSorokina

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/alfeld_sorokina.jl`.
"""
struct AlfeldSorokina
    function AlfeldSorokina()
        return error("AlfeldSorokina is not implemented; see the STATUS: BLOCKED notes in src/alfeld_sorokina.jl.")
    end
end
