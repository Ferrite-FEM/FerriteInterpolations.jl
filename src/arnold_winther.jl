# Arnold-Winther elements (conforming and nonconforming stress elements) (https://defelement.org/elements/arnold-winther.html,
# DefElement: elements/arnold-winther.def).
#
# Arnold-Winther elements (conforming and nonconforming stress elements).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This element is MATRIX-VALUED (symmetric-matrix stress fields in H(div div); also covers the nonconforming variant (nonconforming-arnold-winther.def)). Ferrite has scalar and vector
# interpolations only (ScalarInterpolation / VectorInterpolation) and no
# matrix-valued interpolation type, no symmetric-matrix FEValues storage, and
# no double-Piola mappings. All of that would have to be added upstream
# before this element (and the other stress/strain elements) can exist.

"""
    ArnoldWinther

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/arnold_winther.jl`.
"""
struct ArnoldWinther
    function ArnoldWinther()
        return error("ArnoldWinther is not implemented; see the STATUS: BLOCKED notes in src/arnold_winther.jl.")
    end
end
