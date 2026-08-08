# Arnold-Boffi-Falk element (quadrilateral H(div)) (https://defelement.org/elements/arnold-boffi-falk.html,
# DefElement: elements/arnold-boffi-falk.def).
#
# Arnold-Boffi-Falk element (quadrilateral H(div)).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# The element's DOFs include derivative-type functionals (divergence moments among the cell DOFs (the def's mapping field refers to Kirby's transformation theory rather than a plain Piola map), with shared edge DOFs), which are
# not preserved by any of Ferrite's mappings: the physical functionals relate
# to the reference ones by a CELL-GEOMETRY-DEPENDENT matrix (chain rule
# through J), and since these DOFs are shared between neighboring cells the
# transformation must be coordinated across the mesh. That requires general
# per-cell DOF transformations (R. Kirby, "A general approach to transforming
# finite elements", SMAI-JCM 2018) -- Ferrite's orientation machinery
# (adjust_dofs_during_distribution + get_direction's +-1) cannot express it.
# (Compare src/taylor.jl, where the same mathematics IS implementable only
# because all Taylor DOFs are cell-local.)

"""
    ArnoldBoffiFalk

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/arnold_boffi_falk.jl`.
"""
struct ArnoldBoffiFalk
    function ArnoldBoffiFalk()
        return error("ArnoldBoffiFalk is not implemented; see the STATUS: BLOCKED notes in src/arnold_boffi_falk.jl.")
    end
end
