using Ferrite
using Test
include("test_utils.jl")

# symfem cross-check of Ferrite's own CrouzeixRaviart and RannacherTurek
# (DefElement: crouzeix-raviart, which covers both simplices and hypercubes).

@testset "Ferrite CrouzeixRaviart/RannacherTurek (symfem)" begin
    # The bases agree exactly, up to the facet numbering.
    test_symfem_reference(CrouzeixRaviart{RefTriangle, 1}(), "Crouzeix-Raviart", 1, [2, 1, 0])
    test_symfem_reference(CrouzeixRaviart{RefTetrahedron, 1}(), "Crouzeix-Raviart", 1, [0, 1, 3, 2])
    test_symfem_reference(RannacherTurek{RefQuadrilateral, 1}(), "Rannacher-Turek", 1, [0, 2, 3, 1])
    test_symfem_reference(RannacherTurek{RefHexahedron, 1}(), "Rannacher-Turek", 1, [0, 1, 3, 4, 2, 5])
end
