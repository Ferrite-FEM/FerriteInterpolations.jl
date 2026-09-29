using Ferrite
using Test
include("test_utils.jl")

# symfem cross-check of Ferrite's own BrezziDouglasMarini (degree 1; degree 2
# is `BDM` in this package, see test_bdm.jl). Same conventions as `BDM`:
# edge DOFs along the Ferrite edge direction with outward normals, giving a
# signed permutation of symfem's basis.

@testset "Ferrite BrezziDouglasMarini (symfem)" begin
    test_symfem_reference_vector(
        BrezziDouglasMarini{RefTriangle, 1}(), "BDM", 1,
        [(-1, 4), (-1, 5), (1, 3), (1, 2), (-1, 0), (-1, 1)],
    )
end
