# Gopalakrishnan-Lederer-Schoberl element (https://defelement.org/elements/gopalakrishnan-lederer-schoberl.html,
# DefElement: elements/gopalakrishnan-lederer-schoberl.def).
#
# Gopalakrishnan-Lederer-Schoberl element.
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This element is MATRIX-VALUED (traceless-matrix fields with normal-tangential continuity). Ferrite has scalar and vector
# interpolations only (ScalarInterpolation / VectorInterpolation) and no
# matrix-valued interpolation type, no symmetric-matrix FEValues storage, and
# no double-Piola mappings. All of that would have to be added upstream
# before this element (and the other stress/strain elements) can exist.

"""
    GopalakrishnanLedererSchoberl

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/gopalakrishnan_lederer_schoberl.jl`.
"""
struct GopalakrishnanLedererSchoberl
    function GopalakrishnanLedererSchoberl()
        return error("GopalakrishnanLedererSchoberl is not implemented; see the STATUS: BLOCKED notes in src/gopalakrishnan_lederer_schoberl.jl.")
    end
end
