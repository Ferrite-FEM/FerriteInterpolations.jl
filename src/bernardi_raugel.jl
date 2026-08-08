# Bernardi-Raugel element (P1^d + facet normal bubbles for Stokes) (https://defelement.org/elements/bernardi-raugel.html,
# DefElement: elements/bernardi-raugel.def).
#
# Bernardi-Raugel element (P1^d + facet normal bubbles for Stokes).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# The element's DOFs include derivative-type functionals (vector vertex values mixed with facet-normal bubble moments (the normal direction is cell-geometry dependent under the identity map)), which are
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
    BernardiRaugel

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/bernardi_raugel.jl`.
"""
struct BernardiRaugel
    function BernardiRaugel()
        return error("BernardiRaugel is not implemented; see the STATUS: BLOCKED notes in src/bernardi_raugel.jl.")
    end
end
