# Guzman-Neilan elements (first and second kind) (https://defelement.org/elements/guzman-neilan.html,
# DefElement: elements/guzman-neilan.def).
#
# Guzman-Neilan elements (first and second kind).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This is a MACRO element: its basis functions are piecewise polynomials on a
# barycentric split of the cell (divergence-free piecewise-polynomial bubbles on the Alfeld split (both first and second kind, guzman-neilan2.def)). Ferrite interpolations evaluate a
# single smooth expression per basis function via `reference_shape_value`,
# and quadrature rules know nothing about interior sub-cell boundaries, so
# piecewise bases cannot be represented (values would be wrong at/near the
# internal interfaces and quadrature would not resolve the kinks). Upstream
# support for composite/macro reference cells would be needed.

"""
    GuzmanNeilan

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/guzman_neilan.jl`.
"""
struct GuzmanNeilan
    function GuzmanNeilan()
        return error("GuzmanNeilan is not implemented; see the STATUS: BLOCKED notes in src/guzman_neilan.jl.")
    end
end
