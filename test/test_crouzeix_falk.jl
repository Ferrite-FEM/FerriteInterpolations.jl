using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "CrouzeixFalk" begin
    ip = CrouzeixFalk{RefTriangle, 3}()

    # (i) Interpolation-level tests
    test_interpolation_properties(ip)
    test_partition_of_unity(ip)
    test_polynomial_reproduction(ip, 3; space = :P)
    test_type_genericity(ip)
    test_kronecker_delta(ip) # all DOFs are point evaluations

    # Live symfem cross-check. The permutation (Ferrite DOF -> 0-based symfem
    # DOF) was derived geometrically from the DOF evaluation points: symfem
    # edge DOFs run at 1/4, 1/2, 3/4 from the first towards the second vertex
    # of each symfem edge; Ferrite edge 1 coincides with symfem edge 2 (same
    # direction), Ferrite edge 2 is symfem edge 1 reversed, and Ferrite edge 3
    # is symfem edge 0 (same direction).
    test_symfem_reference(ip, "Crouzeix-Falk", 3, [6, 7, 8, 5, 4, 3, 0, 1, 2, 9])

    # (ii) Integration tests: nodal interpolation of a full cubic is exact on
    # an affine cell.
    @testset "CellValues cubic reproduction" begin
        coords = [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))]
        geo = Lagrange{RefTriangle, 1}()
        f(x) = 1 - 2x[1] + x[2] + x[1] * x[2] - x[1]^2 + x[1]^3 / 3 - x[1]^2 * x[2] + x[2]^3 / 2
        ∇f(x) = Vec((-2 + x[2] - 2x[1] + x[1]^2 - 2x[1] * x[2], 1 + x[1] - x[1]^2 + 3x[2]^2 / 2))
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:3)
        ue = [f(spatial(ξ)) for ξ in Ferrite.reference_coordinates(ip)]
        cv = CellValues(QuadratureRule{RefTriangle}(4), ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) atol = 1.0e-12
            @test function_gradient(cv, qp, ue) ≈ ∇f(x) atol = 1.0e-12
        end
    end

    @testset "FacetValues partition of unity" begin
        coords = [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))]
        fv = FacetValues(FacetQuadratureRule{RefTriangle}(3), ip, Lagrange{RefTriangle, 1}())
        ue = ones(10)
        for facet in 1:3
            reinit!(fv, coords, facet)
            for qp in 1:getnquadpoints(fv)
                @test function_value(fv, qp, ue) ≈ 1.0 atol = 1.0e-13
            end
        end
    end

    # (iii) DofHandler: the edge DOFs are shared between neighboring cells.
    # Two triangles whose shared edge has opposite local orientation; the
    # interpolated function from both sides must agree at the three shared
    # DOF points (only there -- the element is nonconforming, traces are not
    # matched along the whole edge). Exercises dof distribution and the
    # adjust_dofs_during_distribution edge reversal.
    @testset "two-cell shared edge DOFs" begin
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
        @test ndofs(dh) == 2 * 10 - 3
        @test length(intersect(celldofs(dh, 1), celldofs(dh, 2))) == 3
        u = rand(ndofs(dh))
        dofs1 = celldofs(dh, 1)
        dofs2 = celldofs(dh, 2)
        eval_ref(ξ, dofs) = sum(u[dofs[i]] * Ferrite.reference_shape_value(ip, ξ, i) for i in 1:10)
        for t in (1 / 4, 1 / 2, 3 / 4)
            # Shared DOF point at physical (1-t)*n2 + t*n3: cell 1 edge (2,3)
            # param t -> ξ = (0, 1-t); cell 2 (local vertices n2, n4, n3) edge
            # (3,1) param 1-t -> ξ = (1-t, 0)
            v1 = eval_ref(Vec((0.0, 1 - t)), dofs1)
            v2 = eval_ref(Vec((1 - t, 0.0)), dofs2)
            @test v1 ≈ v2 atol = 1.0e-13
        end
    end
end
