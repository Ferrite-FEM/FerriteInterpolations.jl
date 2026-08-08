using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "ConformingCrouzeixRaviart" begin
    ips = (
        ConformingCrouzeixRaviart{RefTriangle, 2}(),
        ConformingCrouzeixRaviart{RefTriangle, 3}(),
        ConformingCrouzeixRaviart{RefTriangle, 4}(),
    )

    # (i) Interpolation-level tests
    for ip in ips
        test_interpolation_properties(ip)
        test_partition_of_unity(ip)
        test_polynomial_reproduction(ip, Ferrite.getorder(ip); space = :P)
        test_type_genericity(ip)
        test_kronecker_delta(ip) # all DOFs are point evaluations
    end

    # Live symfem cross-check. The permutations (Ferrite DOF -> 0-based symfem
    # DOF) were derived geometrically from the DOF evaluation points: Ferrite
    # vertices (1, 2, 3) are symfem vertices (1, 2, 0); Ferrite edge 1
    # coincides with symfem edge 2 (same direction), Ferrite edge 2 is symfem
    # edge 1 reversed, Ferrite edge 3 is symfem edge 0 (same direction);
    # interior DOFs keep symfem's order.
    for (ip, perm) in (
            (ips[1], [1, 2, 0, 5, 4, 3, 6]),
            (ips[2], [1, 2, 0, 7, 8, 6, 5, 3, 4, 9, 10, 11]),
            (ips[3], [1, 2, 0, 9, 10, 11, 8, 7, 6, 3, 4, 5, 12, 13, 14, 15, 16, 17]),
        )
        test_symfem_reference(ip, "conforming Crouzeix-Raviart", Ferrite.getorder(ip), perm)
    end

    # (ii) Integration tests: nodal interpolation of a full degree-k
    # polynomial is exact on an affine cell.
    @testset "CellValues P$k reproduction" for (k, ip) in zip(2:4, ips)
        coords = [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))]
        geo = Lagrange{RefTriangle, 1}()
        f(x) = (1 + x[1] - x[2] / 2)^k
        ∇f(x) = k * (1 + x[1] - x[2] / 2)^(k - 1) * Vec((1.0, -1 / 2))
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:3)
        ue = [f(spatial(ξ)) for ξ in Ferrite.reference_coordinates(ip)]
        cv = CellValues(QuadratureRule{RefTriangle}(2k), ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) atol = 1.0e-11
            @test function_gradient(cv, qp, ue) ≈ ∇f(x) atol = 1.0e-11
        end
    end

    @testset "FacetValues partition of unity: $ip" for ip in ips
        coords = [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))]
        fv = FacetValues(FacetQuadratureRule{RefTriangle}(3), ip, Lagrange{RefTriangle, 1}())
        ue = ones(Ferrite.getnbasefunctions(ip))
        for facet in 1:3
            reinit!(fv, coords, facet)
            for qp in 1:getnquadpoints(fv)
                @test function_value(fv, qp, ue) ≈ 1.0 atol = 1.0e-13
            end
        end
    end

    # (iii) DofHandler + H1 conformity: two triangles whose shared edge has
    # opposite local orientation; the interpolated function must agree along
    # the WHOLE shared edge (the element is conforming: edge traces have
    # degree <= k and are determined by the k+1 shared point values).
    # Exercises dof distribution and, for k >= 3, the
    # adjust_dofs_during_distribution edge reversal.
    @testset "two-cell continuity: $ip" for (k, ip) in zip(2:4, ips)
        nodes = [
            Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
            Node(Vec((0.0, 1.0))), Node(Vec((1.0, 1.0))),
        ]
        # Shared edge: (2, 3) in cell 1, (3, 2) in cell 2 -> opposite orientation
        cells = [Triangle((1, 2, 3)), Triangle((2, 4, 3))]
        grid = Grid(cells, nodes)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        # Shared: 2 vertices + (k - 1) edge-interior DOFs
        N = Ferrite.getnbasefunctions(ip)
        @test ndofs(dh) == 2N - (2 + (k - 1))
        u = rand(ndofs(dh))
        dofs1 = celldofs(dh, 1)
        dofs2 = celldofs(dh, 2)
        eval_ref(ξ, dofs) = sum(u[dofs[i]] * Ferrite.reference_shape_value(ip, ξ, i) for i in 1:N)
        for s in (0.1, 0.15, 0.4, 0.5, 0.8, 0.95)
            # Physical point (1-s, s): cell 1 edge (2,3) param s -> ξ = (0, 1-s),
            # cell 2 (local vertices n2, n4, n3) edge (3,1) param 1-s -> ξ = (1-s, 0)
            v1 = eval_ref(Vec((0.0, 1 - s)), dofs1)
            v2 = eval_ref(Vec((1 - s, 0.0)), dofs2)
            @test v1 ≈ v2 atol = 1.0e-12
        end
    end
end
