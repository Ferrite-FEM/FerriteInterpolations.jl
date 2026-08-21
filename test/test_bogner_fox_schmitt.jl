using FerriteInterpolations
using Ferrite
using Ferrite: reference_shape_value, reference_shape_gradient
using Test
include("test_utils.jl")

@testset "BognerFoxSchmitt" begin
    ip = BognerFoxSchmitt()

    test_interpolation_properties(ip)
    test_polynomial_reproduction(ip, 3; space = :Q)
    test_type_genericity(ip)

    # Live symfem cross-check ("Bogner-Fox-Schmit"): the vertex permutation
    # maps Ferrite's counterclockwise vertex order to symfem's lexicographic
    # one; the scales 1, 2, 4 are the reference-map derivative factors of
    # X = (ξ + 1)/2 for the (value, dX, dY, dXdY) DOFs.
    test_symfem_reference(
        ip, "Bogner-Fox-Schmit", 3,
        [0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, 8, 9, 10, 11];
        scales = repeat([1.0, 2.0, 2.0, 4.0], 4)
    )

    verts = [Vec((-1.0, -1.0)), Vec((1.0, -1.0)), Vec((1.0, 1.0)), Vec((-1.0, 1.0))]
    N(j, ξ) = reference_shape_value(ip, ξ, j)

    # The dof functionals at each vertex are (N, ∂N/∂ξ₁, ∂N/∂ξ₂, ∂²N/∂ξ₁∂ξ₂),
    # so the nodal Kronecker test does not apply.
    @testset "generalized Kronecker property" begin
        for i in 1:16
            for (w, ξw) in enumerate(verts)
                g = reference_shape_gradient(ip, ξw, i)
                H = Ferrite.Tensors.hessian(ξ -> N(i, ξ), ξw)
                for (kk, val) in enumerate((N(i, ξw), g[1], g[2], H[1, 2]))
                    @test val ≈ (i == 4 * (w - 1) + kk ? 1.0 : 0.0) atol = 1.0e-13
                end
            end
        end
    end

    # Nonuniform 2-cell rectangle grid: cell 1 is 1×2, cell 2 is 2×2. Different
    # cell sizes make the tests sensitive to the Jacobian scaling of the
    # derivative basis functions.
    nodes2 = [
        Node(Vec(0.0, 0.0)), Node(Vec(1.0, 0.0)), Node(Vec(3.0, 0.0)),
        Node(Vec(0.0, 2.0)), Node(Vec(1.0, 2.0)), Node(Vec(3.0, 2.0)),
    ]
    grid2 = Grid([Quadrilateral((1, 2, 5, 4)), Quadrilateral((2, 3, 6, 5))], nodes2)
    dh2 = DofHandler(grid2)
    add!(dh2, :u, ip)
    close!(dh2)
    cd1, cd2 = celldofs(dh2, 1), celldofs(dh2, 2)

    @testset "dof distribution" begin
        @test ndofs(dh2) == 24
        @test cd1[5:8] == cd2[1:4]     # shared dof quadruple at grid node 2
        @test cd1[9:12] == cd2[13:16]  # shared dof quadruple at grid node 5
        @test length(unique(vcat(cd1, cd2))) == 24
    end

    # Exact interpolant of w(x,y) = x³y²: dofs are (w, ∂w/∂x, ∂w/∂y, ∂²w/∂x∂y) (physical)
    w_ex(x) = x[1]^3 * x[2]^2
    dw_ex(x) = (3x[1]^2 * x[2]^2, 2x[1]^3 * x[2], 6x[1]^2 * x[2])
    u2 = zeros(ndofs(dh2))
    for cellid in 1:2
        cd = celldofs(dh2, cellid)
        for (v, x) in enumerate(getcoordinates(grid2, cellid))
            wx, wy, wxy = dw_ex(x)
            u2[cd[4 * (v - 1) + 1]] = w_ex(x)
            u2[cd[4 * (v - 1) + 2]] = wx
            u2[cd[4 * (v - 1) + 3]] = wy
            u2[cd[4 * (v - 1) + 4]] = wxy
        end
    end

    @testset "mapping exactness on nonuniform mesh" begin
        cv = CellValues(QuadratureRule{RefQuadrilateral}(4), ip; update_hessians = true)
        for cellid in 1:2
            coords = getcoordinates(grid2, cellid)
            reinit!(cv, coords)
            ue = u2[celldofs(dh2, cellid)]
            for qp in 1:getnquadpoints(cv)
                x = spatial_coordinate(cv, qp, coords)
                wx, wy, wxy = dw_ex(x)
                @test function_value(cv, qp, ue) ≈ w_ex(x) atol = 1.0e-11
                g = function_gradient(cv, qp, ue)
                @test g[1] ≈ wx atol = 1.0e-11
                @test g[2] ≈ wy atol = 1.0e-11
                H = function_hessian(cv, qp, ue)
                @test H[1, 1] ≈ 6x[1] * x[2]^2 atol = 1.0e-10
                @test H[1, 2] ≈ wxy atol = 1.0e-10
                @test H[2, 1] ≈ wxy atol = 1.0e-10
                @test H[2, 2] ≈ 2x[1]^3 atol = 1.0e-10
            end
        end
    end

    @testset "C1 continuity across the shared edge" begin
        urand = rand(ndofs(dh2))
        for η in (-0.8, -0.3, 0.44, 0.9)
            # Right edge of cell 1 and left edge of cell 2 at the same physical point
            cvr = CellValues(QuadratureRule{RefQuadrilateral}([1.0], [Vec((1.0, η))]), ip)
            cvl = CellValues(QuadratureRule{RefQuadrilateral}([1.0], [Vec((-1.0, η))]), ip)
            reinit!(cvr, getcoordinates(grid2, 1))
            reinit!(cvl, getcoordinates(grid2, 2))
            @test function_value(cvr, 1, urand[cd1]) ≈ function_value(cvl, 1, urand[cd2])
            @test function_gradient(cvr, 1, urand[cd1]) ≈ function_gradient(cvl, 1, urand[cd2])
        end
    end

    @testset "FacetValues" begin
        fv = FacetValues(FacetQuadratureRule{RefQuadrilateral}(3), ip)
        for cellid in 1:2, facetnr in 1:4
            coords = getcoordinates(grid2, cellid)
            reinit!(fv, coords, facetnr)
            ue = u2[celldofs(dh2, cellid)]
            for qp in 1:getnquadpoints(fv)
                x = spatial_coordinate(fv, qp, coords)
                wx, wy, _ = dw_ex(x)
                @test function_value(fv, qp, ue) ≈ w_ex(x) atol = 1.0e-11
                g = function_gradient(fv, qp, ue)
                @test g[1] ≈ wx atol = 1.0e-11
                @test g[2] ≈ wy atol = 1.0e-11
            end
        end
    end

    @testset "evaluate_at_grid_nodes and PointEvalHandler" begin
        @test evaluate_at_grid_nodes(dh2, u2, :u) ≈ [w_ex(n.x) for n in nodes2]
        ph = PointEvalHandler(grid2, [Vec((0.3, 0.77)), Vec((2.1, 1.44))])
        @test evaluate_at_points(ph, dh2, u2, :u) ≈ [w_ex(Vec((0.3, 0.77))), w_ex(Vec((2.1, 1.44)))]
    end

    @testset "Dirichlet on value and derivative dof kinds" begin
        grid = generate_grid(Quadrilateral, (4, 4))
        dh = DofHandler(grid)
        add!(dh, :c, ip)
        close!(dh)
        lr = union(getfacetset(grid, "left"), getfacetset(grid, "right"))
        tb = union(getfacetset(grid, "top"), getfacetset(grid, "bottom"))
        ch = ConstraintHandler(dh)
        add!(ch, Dirichlet(:c, lr, (x, t) -> x[2] + t; kind = :derivative_x))
        add!(ch, Dirichlet(:c, tb, Returns(0.0); kind = :derivative_y))
        add!(ch, Dirichlet(:c, union(lr, tb), Returns(0.0); kind = :derivative_xy))
        close!(ch)
        update!(ch, 0.0)
        # 5×5 vertex grid: ∂x dofs on the 2×5 left/right vertices, ∂y dofs on
        # the 2×5 top/bottom vertices, ∂x∂y dofs on all 16 boundary vertices
        @test length(ch.prescribed_dofs) == 10 + 10 + 16
        a = zeros(ndofs(dh))
        apply!(a, ch)
        # ∂c/∂x dof at the top-left corner should have the value y + t = 1.0 + 0.0
        cellid = getncells(grid) - 3 # top-left cell; corner is its 4th vertex
        topleft = celldofs(dh, cellid)[4 * 3 + 2]
        @test a[topleft] ≈ 1.0
        update!(ch, 2.0)
        apply!(a, ch)
        @test a[topleft] ≈ 3.0
    end

    @testset "unsupported usage errors" begin
        @test_throws ArgumentError BognerFoxSchmitt()^2
        @test_throws ArgumentError CellValues(QuadratureRule{RefQuadrilateral}(2), ip, Lagrange{RefQuadrilateral, 2}()^2) # higher order geometry
        @test_throws ArgumentError CellValues(QuadratureRule{RefQuadrilateral}(2), ip, Lagrange{RefQuadrilateral, 1}()^3) # embedded
        @test_throws ErrorException apply_analytical!(zeros(ndofs(dh2)), dh2, :u, x -> 0.0)
        # Non-rectangular (but valid bilinear) cells are caught at reinit! time
        badgrid = Grid(
            [Quadrilateral((1, 2, 3, 4))],
            [Node(Vec(0.0, 0.0)), Node(Vec(1.0, 0.1)), Node(Vec(1.2, 1.0)), Node(Vec(0.0, 1.0))]
        )
        cvb = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
        @test_throws ArgumentError reinit!(cvb, getcoordinates(badgrid, 1))
        # Rotated rectangle: affine, but the dof set is not closed under rotation
        s, c = sincos(π / 6)
        R = Ferrite.Tensors.Tensor{2, 2}((c, s, -s, c))
        rotgrid = Grid(
            [Quadrilateral((1, 2, 3, 4))],
            [Node(R ⋅ Vec(0.0, 0.0)), Node(R ⋅ Vec(1.0, 0.0)), Node(R ⋅ Vec(1.0, 1.0)), Node(R ⋅ Vec(0.0, 1.0))]
        )
        @test_throws ArgumentError reinit!(cvb, getcoordinates(rotgrid, 1))
    end

    @testset "clamped Kirchhoff plate patch test" begin
        # Pure bending patch test: with w = x² on the boundary (value, normal
        # derivative) the Galerkin solution of the biharmonic equation Δ²w = 0
        # reproduces w = x² exactly since it lies in the trial space.
        grid = generate_grid(Quadrilateral, (3, 2), Vec((0.0, 0.0)), Vec((1.5, 1.0)))
        dh = DofHandler(grid)
        add!(dh, :w, ip)
        close!(dh)
        w_a(x) = x[1]^2
        ∂Ω = union((getfacetset(grid, name) for name in ("left", "right", "top", "bottom"))...)
        ch = ConstraintHandler(dh)
        add!(ch, Dirichlet(:w, ∂Ω, x -> w_a(x)))
        add!(ch, Dirichlet(:w, ∂Ω, x -> 2x[1]; kind = :derivative_x))
        add!(ch, Dirichlet(:w, ∂Ω, Returns(0.0); kind = :derivative_y))
        add!(ch, Dirichlet(:w, ∂Ω, Returns(0.0); kind = :derivative_xy))
        close!(ch)
        cv = CellValues(QuadratureRule{RefQuadrilateral}(4), ip; update_hessians = true)
        K = allocate_matrix(dh)
        f = zeros(ndofs(dh))
        assembler = start_assemble(K, f)
        ndof = getnbasefunctions(ip)
        Ke = zeros(ndof, ndof)
        for cell in CellIterator(dh)
            reinit!(cv, cell)
            fill!(Ke, 0.0)
            for qp in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv, qp)
                for i in 1:ndof
                    Hi = shape_hessian(cv, qp, i)
                    for j in 1:ndof
                        Ke[i, j] += (Hi ⊡ shape_hessian(cv, qp, j)) * dΩ
                    end
                end
            end
            assemble!(assembler, celldofs(cell), Ke)
        end
        apply!(K, f, ch)
        w = K \ f
        ph = PointEvalHandler(grid, [Vec((0.35, 0.6)), Vec((1.2, 0.15))])
        @test evaluate_at_points(ph, dh, w, :w) ≈ [0.35^2, 1.2^2] atol = 1.0e-10
    end
end
