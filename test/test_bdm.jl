using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "BDM" begin
    ip = BDM{RefTriangle, 2}()

    # (i) Interpolation-level tests
    test_interpolation_properties(ip)
    test_type_genericity(ip)

    # Live symfem cross-check: the Ferrite basis is a signed permutation of
    # symfem's (derived exactly with sympy from the Ferrite DOF conventions,
    # see src/bdm.jl): Ferrite edge 1 = -(symfem edge 2, weight order
    # {s=0, s=1, s=1/2}), edge 2 = +(symfem edge 1 reversed), edge 3 =
    # -(symfem edge 0), interior unchanged.
    test_symfem_reference_vector(
        ip, "BDM", 2,
        [(-1, 6), (-1, 8), (-1, 7), (1, 4), (1, 5), (1, 3), (-1, 0), (-1, 2), (-1, 1), (1, 9), (1, 10), (1, 11)],
    )

    # Full vector P2 span: every monomial vector (x^a y^b) e_d with
    # a + b <= 2 lies in the span (least-squares fit at random points).
    @testset "P2 vector span" begin
        pts = [sample_reference_point(RefTriangle) for _ in 1:40]
        A = zeros(2 * length(pts), 12)
        for (p, ξ) in pairs(pts), i in 1:12
            v = reference_shape_value(ip, ξ, i)
            A[2p - 1, i] = v[1]
            A[2p, i] = v[2]
        end
        for e in monomial_exponents(:P, 2, 2), d in 1:2
            b = zeros(2 * length(pts))
            for (p, ξ) in pairs(pts)
                b[2p - 2 + d] = prod(ξ .^ e)
            end
            c = A \ b
            @test norm(A * c - b) < 1.0e-10 * max(1, norm(b))
        end
    end

    # (ii) Integration tests on a two-cell mesh whose shared edge has
    # opposite local orientation (cell 1 edge (2,3), cell 2 edge (3,2)).
    nodes = [
        Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
        Node(Vec((0.0, 1.0))), Node(Vec((1.2, 1.1))),
    ]
    grid = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes)
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)

    @testset "DofHandler distribution" begin
        @test ndofs(dh) == 2 * 12 - 3
        @test length(intersect(celldofs(dh, 1), celldofs(dh, 2))) == 3
    end

    # H(div) conformity: the NORMAL component of the interpolated field is
    # continuous across the shared edge (the tangential part may jump).
    @testset "two-cell normal continuity" begin
        u = rand(ndofs(dh))
        fqr = FacetQuadratureRule{RefTriangle}(4)
        fv1 = FacetValues(fqr, ip, Lagrange{RefTriangle, 1}())
        fv2 = FacetValues(fqr, ip, Lagrange{RefTriangle, 1}())
        # Shared edge: facet 2 of cell 1, facet 3 of cell 2
        coords1 = getcoordinates(grid, 1)
        coords2 = getcoordinates(grid, 2)
        reinit!(fv1, getcells(grid, 1), coords1, 2)
        reinit!(fv2, getcells(grid, 2), coords2, 3)
        u1 = u[celldofs(dh, 1)]
        u2 = u[celldofs(dh, 2)]
        for qp1 in 1:getnquadpoints(fv1)
            x1 = spatial_coordinate(fv1, qp1, coords1)
            qp2 = findfirst(qp -> norm(spatial_coordinate(fv2, qp, coords2) - x1) < 1.0e-12, 1:getnquadpoints(fv2))
            @test qp2 !== nothing
            n1 = getnormal(fv1, qp1)
            v1 = function_value(fv1, qp1, u1)
            v2 = function_value(fv2, qp2, u2)
            # Same physical point, opposite outward normals
            @test v1 ⋅ n1 ≈ v2 ⋅ (-getnormal(fv2, qp2)) rtol = 1.0e-10 atol = 1.0e-12
            @test getnormal(fv2, qp2) ≈ -n1 atol = 1.0e-12
        end
    end

    # Divergence theorem on a single (affine) cell: int div(v) dV equals the
    # total boundary flux, for every basis function. This exercises the
    # contravariant Piola mapping of values and gradients consistently.
    @testset "divergence theorem" begin
        coords = getcoordinates(grid, 2)
        cell = getcells(grid, 2)
        cv = CellValues(QuadratureRule{RefTriangle}(4), ip, Lagrange{RefTriangle, 1}())
        reinit!(cv, cell, coords)
        fv = FacetValues(FacetQuadratureRule{RefTriangle}(4), ip, Lagrange{RefTriangle, 1}())
        for i in 1:12
            vol = sum(shape_divergence(cv, qp, i) * getdetJdV(cv, qp) for qp in 1:getnquadpoints(cv))
            flux = 0.0
            for facet in 1:3
                reinit!(fv, cell, coords, facet)
                flux += sum(
                    (shape_value(fv, qp, i) ⋅ getnormal(fv, qp)) * getdetJdV(fv, qp)
                        for qp in 1:getnquadpoints(fv)
                )
            end
            @test vol ≈ flux atol = 1.0e-11
        end
    end
end
