using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "DPC" begin
    ips = (
        DPC{RefQuadrilateral, 1}(), DPC{RefQuadrilateral, 2}(), DPC{RefQuadrilateral, 3}(),
        DPC{RefHexahedron, 1}(), DPC{RefHexahedron, 2}(),
    )

    # (i) Interpolation-level tests. The span is the TOTAL-degree space P_k
    # (not Q_k), which the polynomial-reproduction helper checks exactly.
    for ip in ips
        test_interpolation_properties(ip)
        test_partition_of_unity(ip)
        test_polynomial_reproduction(ip, Ferrite.getorder(ip); space = :P)
        test_type_genericity(ip)
        test_kronecker_delta(ip) # all DOFs are point evaluations
    end

    # Live symfem cross-check. All DOFs are cell DOFs kept in symfem's order
    # (the helper maps the hypercube coordinates), so the permutations are
    # trivial.
    for ip in ips
        N = getnbasefunctions(ip)
        test_symfem_reference(ip, "dPc", Ferrite.getorder(ip), collect(0:(N - 1)))
    end

    # (ii) Integration tests: nodal interpolation of a degree-k polynomial is
    # exact on affinely-mapped cells (parallelogram / parallelepiped -- for
    # non-affine mappings P_k is not reproduced, cf. the mapping note in
    # src/dpc.jl).
    quad_coords = [Vec((0.0, 0.0)), Vec((2.0, 0.3)), Vec((2.4, 1.8)), Vec((0.4, 1.5))]
    hex_coords = begin
        a, u, v, w = Vec((0.0, 0.0, 0.0)), Vec((1.1, 0.1, 0.0)), Vec((0.2, 1.3, 0.1)), Vec((0.0, 0.2, 0.9))
        [a, a + u, a + u + v, a + v, a + w, a + u + w, a + u + v + w, a + v + w]
    end
    @testset "CellValues P$(Ferrite.getorder(ip)) reproduction: $ip" for ip in ips
        shape = getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        k = Ferrite.getorder(ip)
        geo = Lagrange{shape, 1}()
        coords = dim == 2 ? quad_coords : hex_coords
        b = Vec{dim}(ntuple(i -> Float64(i), dim))
        f(x) = (1 + x ⋅ b / 4)^k
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:length(coords))
        ue = [f(spatial(ξ)) for ξ in Ferrite.reference_coordinates(ip)]
        cv = CellValues(QuadratureRule{shape}(k + 1), ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) atol = 1.0e-11
        end
    end

    # (iii) DofHandler: all DOFs are cell DOFs, nothing is shared.
    @testset "two-cell independence" begin
        nodes = [Node(Vec((x, y))) for y in (0.0, 1.0) for x in (0.0, 1.0, 2.0)]
        grid = Grid([Quadrilateral((1, 2, 5, 4)), Quadrilateral((2, 3, 6, 5))], nodes)
        for ip in (DPC{RefQuadrilateral, 2}(), DPC{RefQuadrilateral, 3}())
            dh = DofHandler(grid)
            add!(dh, :u, ip)
            close!(dh)
            @test ndofs(dh) == 2 * getnbasefunctions(ip)
            @test isempty(intersect(celldofs(dh, 1), celldofs(dh, 2)))
        end
    end

    # Vector dPc is DPC^vdim via Ferrite's VectorizedInterpolation: check DOF
    # counts, distribution, and that a constant vector field is represented.
    @testset "vector dPc: DPC^vdim" begin
        ip = DPC{RefQuadrilateral, 2}()^2
        @test getnbasefunctions(ip) == 12
        nodes = [Node(Vec((x, y))) for y in (0.0, 1.0) for x in (0.0, 1.0, 2.0)]
        grid = Grid([Quadrilateral((1, 2, 5, 4)), Quadrilateral((2, 3, 6, 5))], nodes)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        @test ndofs(dh) == 2 * 12
        cv = CellValues(QuadratureRule{RefQuadrilateral}(3), ip, Lagrange{RefQuadrilateral, 1}())
        reinit!(cv, quad_coords)
        ue = ones(12) # every node carries (1, 1)
        for qp in 1:getnquadpoints(cv)
            @test function_value(cv, qp, ue) ≈ Vec((1.0, 1.0)) atol = 1.0e-13
        end
    end
end
