using Ferrite
using Test
include("test_utils.jl")

# symfem cross-check of Ferrite's own Serendipity. Ferrite's edge DOFs are
# point evaluations at the edge midpoints, DefElement's are integral moments,
# so the bases differ (only the edge functions coincide, up to a factor 2/3).
# The span check plus nodality (Ferrite's own test suite checks
# N_i(ξ_j) = δ_ij at reference_coordinates) pins down the basis completely.

@testset "Ferrite Serendipity (symfem)" begin
    for ip in (Serendipity{RefQuadrilateral, 2}(), Serendipity{RefHexahedron, 2}())
        test_symfem_span(ip, "serendipity", 2)
        test_kronecker_delta(ip)
    end
end
