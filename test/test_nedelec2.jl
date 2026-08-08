using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "NedelecSecondKind" begin
    ips = (NedelecSecondKind{RefTriangle, 1}(), NedelecSecondKind{RefTriangle, 2}())

    # (i) Interpolation-level tests
    for ip in ips
        test_interpolation_properties(ip)
        test_type_genericity(ip)
    end

    # Live symfem cross-check ("N2curl"): the Ferrite bases are signed
    # permutations of symfem's (derived exactly with sympy, see
    # src/nedelec2.jl).
    test_symfem_reference_vector(
        ips[1], "N2curl", 1,
        [(1, 4), (1, 5), (-1, 3), (-1, 2), (1, 0), (1, 1)],
    )
    test_symfem_reference_vector(
        ips[2], "N2curl", 2,
        [(1, 6), (1, 8), (1, 7), (-1, 4), (-1, 5), (-1, 3), (1, 0), (1, 2), (1, 1), (1, 9), (1, 10), (1, 11)],
    )

    # Full vector P_k span
    @testset "P$k vector span" for (k, ip) in zip(1:2, ips)
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
    end

    # (ii) Integration tests on a two-cell mesh with an opposite-orientation
    # shared edge.
    nodes = [
        Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
        Node(Vec((0.0, 1.0))), Node(Vec((1.2, 1.1))),
    ]
    grid = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes)

    @testset "DofHandler + tangential continuity: $ip" for (k, ip) in zip(1:2, ips)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        N = getnbasefunctions(ip)
        @test ndofs(dh) == 2N - (k + 1)
        @test length(intersect(celldofs(dh, 1), celldofs(dh, 2))) == k + 1
        u = rand(ndofs(dh))
        fqr = FacetQuadratureRule{RefTriangle}(4)
        fv1 = FacetValues(fqr, ip, Lagrange{RefTriangle, 1}())
        fv2 = FacetValues(fqr, ip, Lagrange{RefTriangle, 1}())
        coords1 = getcoordinates(grid, 1)
        coords2 = getcoordinates(grid, 2)
        # Shared edge: facet 2 of cell 1, facet 3 of cell 2
        reinit!(fv1, getcells(grid, 1), coords1, 2)
        reinit!(fv2, getcells(grid, 2), coords2, 3)
        u1 = u[celldofs(dh, 1)]
        u2 = u[celldofs(dh, 2)]
        for qp1 in 1:getnquadpoints(fv1)
            x1 = spatial_coordinate(fv1, qp1, coords1)
            qp2 = findfirst(qp -> norm(spatial_coordinate(fv2, qp, coords2) - x1) < 1.0e-12, 1:getnquadpoints(fv2))
            @test qp2 !== nothing
            # H(curl): the TANGENTIAL component is continuous (the normal
            # component may jump). Same physical tangent for both sides.
            n1 = getnormal(fv1, qp1)
            t = Vec(-n1[2], n1[1])
            v1 = function_value(fv1, qp1, u1)
            v2 = function_value(fv2, qp2, u2)
            @test v1 ⋅ t ≈ v2 ⋅ t rtol = 1.0e-10 atol = 1.0e-12
        end
    end

    # Stokes/curl theorem on a single (affine) cell: int curl(v) dA equals
    # the counterclockwise circulation, for every basis function (exercises
    # the covariant Piola mapping of values and gradients consistently).
    @testset "curl theorem: $ip" for ip in ips
        coords = getcoordinates(grid, 2)
        cell = getcells(grid, 2)
        N = getnbasefunctions(ip)
        cv = CellValues(QuadratureRule{RefTriangle}(4), ip, Lagrange{RefTriangle, 1}())
        reinit!(cv, cell, coords)
        fv = FacetValues(FacetQuadratureRule{RefTriangle}(4), ip, Lagrange{RefTriangle, 1}())
        for i in 1:N
            curlint = sum(
                (g = shape_gradient(cv, qp, i); (g[2, 1] - g[1, 2]) * getdetJdV(cv, qp))
                    for qp in 1:getnquadpoints(cv)
            )
            circ = 0.0
            for facet in 1:3
                reinit!(fv, cell, coords, facet)
                for qp in 1:getnquadpoints(fv)
                    n = getnormal(fv, qp)
                    t = Vec(-n[2], n[1]) # counterclockwise tangent
                    circ += (shape_value(fv, qp, i) ⋅ t) * getdetJdV(fv, qp)
                end
            end
            @test curlint ≈ circ atol = 1.0e-11
        end
    end
end
