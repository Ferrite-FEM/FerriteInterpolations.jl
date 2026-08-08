# Hellan-Herrmann-Johnson element (plate bending) (https://defelement.org/elements/hellan-herrmann-johnson.html,
# DefElement: elements/hellan-herrmann-johnson.def).
#
# Hellan-Herrmann-Johnson element (plate bending).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# This element is MATRIX-VALUED (symmetric-matrix moment fields with normal-normal continuity). Ferrite has scalar and vector
# interpolations only (ScalarInterpolation / VectorInterpolation) and no
# matrix-valued interpolation type, no symmetric-matrix FEValues storage, and
# no double-Piola mappings. All of that would have to be added upstream
# before this element (and the other stress/strain elements) can exist.

"""
    HellanHerrmannJohnson

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/hellan_herrmann_johnson.jl`.
"""
struct HellanHerrmannJohnson
    function HellanHerrmannJohnson()
        return error("HellanHerrmannJohnson is not implemented; see the STATUS: BLOCKED notes in src/hellan_herrmann_johnson.jl.")
    end
end
