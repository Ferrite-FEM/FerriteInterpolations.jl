# Shared test helpers. Element-specific data (symfem reference tables etc.)
# lives in the per-element test files.

using Ferrite
using Ferrite: getrefshape, getnbasefunctions, reference_shape_value
using LinearAlgebra: norm
using Test

# Ferrite's generic interpolation property checker (test_interpolation_properties).
include(joinpath(pkgdir(Ferrite), "test", "interpolation_test_utils.jl"))

# Sample a (not necessarily uniformly distributed) random point strictly inside
# the reference cell.
function sample_reference_point(::Type{Ferrite.RefHypercube{dim}}) where {dim}
    return Vec{dim}(ntuple(_ -> 2 * rand() - 1, dim))
end
function sample_reference_point(::Type{Ferrite.RefSimplex{dim}}) where {dim}
    # Dirichlet(1, ..., 1) via normalized exponentials; the reference simplex is
    # {ξ ≥ 0, sum(ξ) ≤ 1} for both RefTriangle and RefTetrahedron.
    w = ntuple(_ -> -log(rand()), dim + 1)
    s = sum(w)
    return Vec{dim}(ntuple(i -> w[i] / s, dim))
end
function sample_reference_point(::Type{RefPrism})
    tri = sample_reference_point(RefTriangle)
    return Vec{3}((tri[1], tri[2], rand()))
end
function sample_reference_point(::Type{RefPyramid})
    z = rand()
    return Vec{3}((rand() * (1 - z), rand() * (1 - z), z))
end

# Σᵢ Nᵢ(ξ) == 1 at random reference points (any partition-of-unity element).
function test_partition_of_unity(ip; npoints = 20)
    @testset "partition of unity: $ip" begin
        for _ in 1:npoints
            ξ = sample_reference_point(getrefshape(ip))
            s = sum(i -> reference_shape_value(ip, ξ, i), 1:getnbasefunctions(ip))
            @test s ≈ 1 atol = 1.0e-12
        end
    end
end

# Nᵢ(ξⱼ) == δᵢⱼ at reference_coordinates (nodal elements only).
function test_kronecker_delta(ip)
    @testset "Kronecker delta: $ip" begin
        coords = Ferrite.reference_coordinates(ip)
        N = getnbasefunctions(ip)
        @test length(coords) == N
        for (j, ξ) in pairs(coords), i in 1:N
            @test reference_shape_value(ip, ξ, i) ≈ (i == j ? 1.0 : 0.0) atol = 1.0e-12
        end
    end
end

# Exponent tuples for the monomial basis of P_degree (:P, total degree) or
# Q_degree (:Q, per-variable degree) in `dim` variables.
function monomial_exponents(space::Symbol, dim::Int, degree::Int)
    ranges = ntuple(_ -> 0:degree, dim)
    exps = vec(collect(Iterators.product(ranges...)))
    space === :Q && return exps
    space === :P && return filter(e -> sum(e) <= degree, exps)
    return error("unknown polynomial space $space")
end

# Check that every monomial of P_degree/Q_degree lies in the span of the shape
# functions: least-squares fit at random points, residual ≈ 0. Catches basis
# transcription typos without needing the dual basis.
function test_polynomial_reproduction(ip, degree::Int; space::Symbol)
    @testset "$space$degree ⊆ span: $ip" begin
        shape = getrefshape(ip)
        N = getnbasefunctions(ip)
        npts = 3 * N + 10
        points = [sample_reference_point(shape) for _ in 1:npts]
        A = [reference_shape_value(ip, ξ, i) for ξ in points, i in 1:N]
        for e in monomial_exponents(space, Ferrite.getrefdim(ip), degree)
            b = [prod(ξ .^ e) for ξ in points]
            c = A \ b
            @test norm(A * c - b) < 1.0e-10 * max(1, norm(b))
        end
    end
end

# Evaluate value and AD gradient with different float types (smoke test for
# type-generic implementations).
function test_type_genericity(ip)
    @testset "type genericity: $ip" begin
        shape = getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        for T in (Float32, Float64)
            ξ = Vec{dim, T}(sample_reference_point(shape))
            for i in 1:getnbasefunctions(ip)
                v = reference_shape_value(ip, ξ, i)
                @test v isa Ferrite.shape_value_type(ip, T)
                g = Ferrite.reference_shape_gradient(ip, ξ, i)
                @test all(isfinite, g)
            end
        end
    end
end
