# Brezzi-Douglas-Duran-Fortin element (hexahedral H(div)) (https://defelement.org/elements/brezzi-douglas-duran-fortin.html,
# DefElement: elements/brezzi-douglas-duran-fortin.def).
#
# Brezzi-Douglas-Duran-Fortin element (hexahedral H(div)).
#
# STATUS: BLOCKED -- what must change upstream in Ferrite
#
# The element lives on the hexahedron with order+1 >= 3 normal moments per (quadrilateral) face. Ferrite's face
# orientation handling covers only the Lagrange point-lattice permutations;
# moment DOFs on faces (more than one DOF per face) need general
# orientation-dependent face-DOF transformations, which do not exist yet.
# (This is the same reason the in-tree 3D H(div)/H(curl) families stop at
# the lowest degree.)

"""
    BrezziDouglasDuranFortin

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/brezzi_douglas_duran_fortin.jl`.
"""
struct BrezziDouglasDuranFortin
    function BrezziDouglasDuranFortin()
        return error("BrezziDouglasDuranFortin is not implemented; see the STATUS: BLOCKED notes in src/brezzi_douglas_duran_fortin.jl.")
    end
end
