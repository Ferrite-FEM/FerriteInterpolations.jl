using Ferrite
using Test
include("test_utils.jl")

# symfem cross-check of Ferrite's own RaviartThomas. Ferrite degree k is
# symfem/DefElement degree k - 1 ("Raviart-Thomas" on simplices, "Qdiv" on
# hypercubes). The Ferrite bases are signed permutations of symfem's; on the
# hypercubes additionally scaled by 2^(1 - dim), the contravariant Piola factor
# J / det(J) of the map [0, 1]^dim -> [-1, 1]^dim between the reference cells.

@testset "Ferrite RaviartThomas (symfem)" begin
    test_symfem_reference_vector(
        RaviartThomas{RefTriangle, 1}(), "Raviart-Thomas", 0,
        [(-1, 2), (1, 1), (-1, 0)],
    )
    test_symfem_reference_vector(
        RaviartThomas{RefTriangle, 2}(), "Raviart-Thomas", 1,
        [(-1, 4), (-1, 5), (1, 3), (1, 2), (-1, 0), (-1, 1), (1, 6), (1, 7)],
    )
    test_symfem_reference_vector(
        RaviartThomas{RefTetrahedron, 1}(), "Raviart-Thomas", 0,
        [(-1, 0), (1, 1), (1, 3), (-1, 2)],
    )
    test_symfem_reference_vector(
        RaviartThomas{RefQuadrilateral, 1}(), "Qdiv", 0,
        [(-1 // 2, 0), (-1 // 2, 2), (1 // 2, 3), (1 // 2, 1)],
    )
    test_symfem_reference_vector(
        RaviartThomas{RefHexahedron, 1}(), "Qdiv", 0,
        [(-1 // 4, 0), (1 // 4, 1), (1 // 4, 3), (-1 // 4, 4), (-1 // 4, 2), (1 // 4, 5)],
    )
end
