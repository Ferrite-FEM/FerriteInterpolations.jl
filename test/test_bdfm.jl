using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "BDFM" begin
    ips = (BDFM{RefTriangle, 1}(), BDFM{RefTriangle, 2}())

    # (i) Interpolation-level tests
    for ip in ips
        test_interpolation_properties(ip)
        test_type_genericity(ip)
    end

    # Live symfem cross-check ("BDFM"): signed permutations derived with
    # sympy from the Ferrite DOF conventions (see src/bdfm.jl).
    test_symfem_reference_vector(
        ips[1], "BDFM", 1,
        [(-1, 4), (-1, 5), (1, 3), (1, 2), (-1, 0), (-1, 1), (1, 6), (1, 7), (1, 8)],
    )
    test_symfem_reference_vector(
        ips[2], "BDFM", 2,
        [
            (-1, 6), (-1, 8), (-1, 7), (1, 4), (1, 5), (1, 3), (-1, 0), (-1, 2), (-1, 1),
            (1, 9), (1, 10), (1, 11), (1, 12), (1, 13), (1, 14), (1, 15), (1, 16),
        ],
    )

    # The defining property of BDFM: the space is P_(k+1)^2 restricted to
    # normal traces of degree k; the vector P_k space is contained.
    @testset "P$k vector span: $ip" for (k, ip) in zip(1:2, ips)
        N = getnbasefunctions(ip)
        pts = [sample_reference_point(RefTriangle) for _ in 1:4N]
        A = zeros(2 * length(pts), N)
        for (p, ξ) in pairs(pts), i in 1:N
            v = reference_shape_value(ip, ξ, i)
            A[2p - 1, i] = v[1]
            A[2p, i] = v[2]
        end
        for e in monomial_exponents(:P, 2, k), d in 1:2
            b = zeros(2 * length(pts))
            for (p, ξ) in pairs(pts)
                b[2p - 2 + d] = prod(ξ .^ e)
            end
            c = A \ b
            @test norm(A * c - b) < 1.0e-10 * max(1, norm(b))
        end
        # Normal trace on each edge has degree <= k: check that N_i ⋅ n is
        # orthogonal to nothing extra -- practically, that on edge 3 (y = 0,
        # n = (0,-1)) the trace s -> -N_i[2](s, 0) has vanishing degree-(k+1)
        # Legendre moment.
        for i in 1:N
            # Sample the trace and fit a polynomial of degree k+1: the
            # leading coefficient must vanish.
            ss = range(0, 1; length = k + 3)
            vals = [-reference_shape_value(ip, Vec((s, 0.0)), i)[2] for s in ss]
            V = [s^p for s in ss, p in 0:(k + 2 - 1)]
            coeffs = V \ vals
            @test abs(coeffs[end]) < 1.0e-10
        end
    end

    # (ii) Integration tests
    nodes = [
        Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
        Node(Vec((0.0, 1.0))), Node(Vec((1.2, 1.1))),
    ]
    grid = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes)
    for (k, ip) in zip(1:2, ips)
        test_hdiv_two_cell(ip, grid, 2, 3, k + 1)
        test_divergence_theorem(ip, getcells(grid, 2), getcoordinates(grid, 2))
    end
end
