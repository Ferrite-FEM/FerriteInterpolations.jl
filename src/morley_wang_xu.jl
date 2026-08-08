# Morley-Wang-Xu element (https://defelement.org/elements/morley-wang-xu.html,
# DefElement: elements/morley-wang-xu.def).
#
# Morley-Wang-Xu element.
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# The element's DOFs include derivative-type functionals (facet integral moments of normal derivatives), which are
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
    MorleyWangXu

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/morley_wang_xu.jl`.
"""
struct MorleyWangXu
    function MorleyWangXu()
        return error("MorleyWangXu is not implemented; see the STATUS: BLOCKED notes in src/morley_wang_xu.jl.")
    end
end
