using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "Radau" begin
    ip = Radau{RefLine, 2}()

    # (i) Interpolation-level tests
    test_interpolation_properties(ip)
    test_partition_of_unity(ip)
    test_polynomial_reproduction(ip, 2; space = :P)
    test_type_genericity(ip)
    test_kronecker_delta(ip) # all DOFs are point evaluations

    # Live symfem cross-check ("Lagrange variant=radau"): same DOF order
    # (both vertices, then the interior point).
    test_symfem_reference(ip, "Lagrange", 2, [0, 1, 2]; variant = "radau")

    # (ii) Integration test: nodal interpolation of a quadratic is exact on
    # an affinely-mapped (stretched) interval.
    @testset "CellValues quadratic reproduction" begin
        coords = [Vec((-0.3,)), Vec((1.7,))]
        geo = Lagrange{RefLine, 1}()
        f(x) = 1 - 2x[1] + x[1]^2 / 2
        df(x) = -2 + x[1]
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:2)
        ue = [f(spatial(ξ)) for ξ in Ferrite.reference_coordinates(ip)]
        cv = CellValues(QuadratureRule{RefLine}(3), ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) atol = 1.0e-12
            @test function_gradient(cv, qp, ue)[1] ≈ df(x) atol = 1.0e-12
        end
    end

    # (iii) DofHandler: two line cells share the middle vertex DOF; the
    # interpolated field is continuous there.
    @testset "two-cell vertex sharing" begin
        nodes = [Node(Vec((0.0,))), Node(Vec((1.0,))), Node(Vec((2.5,)))]
        grid = Grid([Line((1, 2)), Line((2, 3))], nodes)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        @test ndofs(dh) == 2 * 3 - 1
        dofs1, dofs2 = celldofs(dh, 1), celldofs(dh, 2)
        @test length(intersect(dofs1, dofs2)) == 1
        u = rand(ndofs(dh))
        # Shared vertex: local vertex 2 of cell 1 (ξ = 1), local vertex 1 of
        # cell 2 (ξ = -1)
        v1 = sum(u[dofs1[i]] * reference_shape_value(ip, Vec((1.0,)), i) for i in 1:3)
        v2 = sum(u[dofs2[i]] * reference_shape_value(ip, Vec((-1.0,)), i) for i in 1:3)
        @test v1 ≈ v2 atol = 1.0e-13
    end
end
