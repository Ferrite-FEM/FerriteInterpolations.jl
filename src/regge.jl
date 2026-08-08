# Regge element (metric/strain element) (https://defelement.org/elements/regge.html,
# DefElement: elements/regge.def).
#
# Regge element (metric/strain element).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This element is MATRIX-VALUED (symmetric-matrix fields in H(curl curl), double-covariant Piola mapping). Ferrite has scalar and vector
# interpolations only (ScalarInterpolation / VectorInterpolation) and no
# matrix-valued interpolation type, no symmetric-matrix FEValues storage, and
# no double-Piola mappings. All of that would have to be added upstream
# before this element (and the other stress/strain elements) can exist.

"""
    Regge

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/regge.jl`.
"""
struct Regge
    function Regge()
        return error("Regge is not implemented; see the STATUS: BLOCKED notes in src/regge.jl.")
    end
end
