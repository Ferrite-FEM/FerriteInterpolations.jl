using FerriteInterpolations
using Ferrite
using Test

include("test_utils.jl")

@testset "FerriteInterpolations" begin
    # One test file per element, mirroring src/
    include("test_bernstein.jl")
end
