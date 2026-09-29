using Ferrite
using Test
include("test_utils.jl")

# symfem cross-check of Ferrite's own Nedelec (first kind). Ferrite degree k
# is symfem/DefElement degree k - 1 ("Nedelec" on simplices, "Qcurl" on
# hypercubes). The Ferrite bases are signed permutations of symfem's; on the
# hypercubes additionally scaled by 1/2, the covariant Piola factor J^-T of
# the map [0, 1]^dim -> [-1, 1]^dim between the reference cells.

@testset "Ferrite Nedelec (symfem)" begin
    test_symfem_reference_vector(
        Nedelec{RefTriangle, 1}(), "Nedelec", 0,
        [(1, 2), (-1, 1), (1, 0)],
    )
    test_symfem_reference_vector(
        Nedelec{RefTriangle, 2}(), "Nedelec", 1,
        [(1, 4), (1, 5), (-1, 3), (-1, 2), (1, 0), (1, 1), (1, 6), (1, 7)],
    )
    test_symfem_reference_vector(
        Nedelec{RefTetrahedron, 1}(), "Nedelec", 0,
        [(1, 0), (1, 3), (-1, 1), (1, 2), (1, 4), (1, 5)],
    )
    test_symfem_reference_vector(
        Nedelec{RefQuadrilateral, 1}(), "Qcurl", 0,
        [(1 // 2, 0), (1 // 2, 2), (-1 // 2, 3), (-1 // 2, 1)],
    )
    test_symfem_reference_vector(
        Nedelec{RefHexahedron, 1}(), "Qcurl", 0,
        [
            (1 // 2, 0), (1 // 2, 3), (-1 // 2, 5), (-1 // 2, 1), (1 // 2, 8), (1 // 2, 10),
            (-1 // 2, 11), (-1 // 2, 9), (1 // 2, 2), (1 // 2, 4), (1 // 2, 7), (1 // 2, 6),
        ],
    )
end
