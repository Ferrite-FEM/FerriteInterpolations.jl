using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "FortinSoulie" begin
    ip = FortinSoulie{RefTriangle, 2}()

    # (i) Interpolation-level tests
    test_interpolation_properties(ip)
    test_partition_of_unity(ip)
    test_polynomial_reproduction(ip, 2; space = :P)
    test_type_genericity(ip)
    test_kronecker_delta(ip) # all DOFs are point evaluations

    # Live symfem cross-check. The permutation (Ferrite DOF -> 0-based symfem
    # DOF) was derived geometrically from the DOF evaluation points: symfem
    # DOFs 2 and 3 sit at (0, 1/3) and (0, 2/3) on the Ferrite edge (2, 3),
    # which runs in the opposite direction, hence the swap.
    test_symfem_reference(ip, "Fortin-Soulie", 2, [0, 1, 3, 2, 4, 5])

    # (ii) Integration tests: nodal interpolation of a full quadratic is exact
    # on an affine cell.
    @testset "CellValues quadratic reproduction" begin
        coords = [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))]
        geo = Lagrange{RefTriangle, 1}()
        f(x) = 1 + 2x[1] - x[2] + x[1]^2 / 2 + x[1] * x[2] - 3x[2]^2 / 10
        ∇f(x) = Vec((2 + x[1] + x[2], -1 + x[1] - 3x[2] / 5))
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:3)
        ue = [f(spatial(ξ)) for ξ in Ferrite.reference_coordinates(ip)]
        cv = CellValues(QuadratureRule{RefTriangle}(3), ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) atol = 1.0e-12
            @test function_gradient(cv, qp, ue) ≈ ∇f(x) atol = 1.0e-12
        end
    end

    @testset "FacetValues partition of unity" begin
        coords = [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))]
        fv = FacetValues(FacetQuadratureRule{RefTriangle}(2), ip, Lagrange{RefTriangle, 1}())
        ue = ones(6)
        for facet in 1:3
            reinit!(fv, coords, facet)
            for qp in 1:getnquadpoints(fv)
                @test function_value(fv, qp, ue) ≈ 1.0 atol = 1.0e-13
            end
        end
    end

    # Two-cell test: all DOFs are cell DOFs (like DiscontinuousLagrange), so
    # nothing is shared between neighboring cells.
    @testset "two-cell independence" begin
        nodes = [
            Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
            Node(Vec((0.0, 1.0))), Node(Vec((1.0, 1.0))),
        ]
        grid = Grid([Triangle((1, 2, 3)), Triangle((3, 2, 4))], nodes)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        @test ndofs(dh) == 12
        @test isempty(intersect(celldofs(dh, 1), celldofs(dh, 2)))
    end
end
