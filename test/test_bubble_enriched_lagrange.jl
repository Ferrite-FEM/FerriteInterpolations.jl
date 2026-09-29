using Ferrite
using Test
include("test_utils.jl")

# symfem cross-check of Ferrite's own BubbleEnrichedLagrange.

@testset "Ferrite BubbleEnrichedLagrange (symfem)" begin
    # Vertex functions up to the vertex numbering, bubble last in both.
    test_symfem_reference(BubbleEnrichedLagrange{RefTriangle, 1}(), "bubble enriched Lagrange", 1, [1, 2, 0, 3])

    # Vectorized ("bubble enriched vector Lagrange"), components innermost.
    test_symfem_reference_vector(
        BubbleEnrichedLagrange{RefTriangle, 1}()^2, "bubble enriched vector Lagrange", 1,
        [(1, 2), (1, 3), (1, 4), (1, 5), (1, 0), (1, 1), (1, 6), (1, 7)],
    )
end
