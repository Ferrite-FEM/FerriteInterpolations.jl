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

    # Reference tabulation generated with symfem 2025.8.0
    # (dev/generate_fortin_soulie_table.py). Points in Ferrite reference
    # coordinates (identical to symfem's for the triangle) and values in
    # Ferrite DOF order; the generator assigns symfem DOFs to Ferrite entities
    # geometrically by their evaluation points (symfem DOF order [0, 1, 3, 2, 4, 5]).
    table = [
        (Vec((0.3, 0.2)), [-0.15, 0.06, -0.18, 0.27, 0.28, 0.72]),
        (Vec((0.1, 0.6)), [-0.075, 0.12, 0.6, -0.015, -0.08, 0.45]),
        (Vec((0.25, 0.25)), [-0.15234375, 0.046875, -0.140625, 0.31640625, 0.15625, 0.7734375]),
        (Vec((0.05, 0.9)), [-0.286875, 1.1475, 1.1475, -0.286875, 0.595, -1.31625]),
    ]
    @testset "symfem reference values" begin
        for (ξ, expected) in table
            vals = [Ferrite.reference_shape_value(ip, ξ, i) for i in 1:6]
            @test vals ≈ expected atol = 1.0e-13
        end
    end

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
