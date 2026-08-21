using FerriteInterpolations
using Ferrite
using Ferrite: reference_shape_value, reference_shape_gradient
using Test
include("test_utils.jl")
# Both Ferrite (the interval/tensor-product element) and FerriteInterpolations
# export `Hermite`; be explicit about which one is under test here.
using FerriteInterpolations: Hermite, HermiteDofTransformation

@testset "Hermite (triangle)" begin
    ip = Hermite{RefTriangle, 3}()

    test_interpolation_properties(ip)
    test_polynomial_reproduction(ip, 3; space = :P)
    test_type_genericity(ip)

    # Live symfem cross-check ("Hermite"): identical reference coordinates,
    # only the vertex numbering differs (Ferrite (1,0),(0,1),(0,0) is symfem's
    # vertices 1, 2, 0).
    test_symfem_reference(ip, "Hermite", 3, [3, 4, 5, 6, 7, 8, 0, 1, 2, 9])

    verts = [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0))]
    N(j, ξ) = reference_shape_value(ip, ξ, j)

    # The dof functionals are (N, ∂N/∂ξ₁, ∂N/∂ξ₂) at the vertices plus the
    # value at the centroid, so the nodal Kronecker test does not apply.
    @testset "generalized Kronecker property" begin
        cen = Vec((1 / 3, 1 / 3))
        for i in 1:10
            for (w, ξw) in enumerate(verts)
                g = reference_shape_gradient(ip, ξw, i)
                for (kk, val) in enumerate((N(i, ξw), g[1], g[2]))
                    @test val ≈ (i == 3 * (w - 1) + kk ? 1.0 : 0.0) atol = 1.0e-13
                end
            end
            @test N(i, cen) ≈ (i == 10 ? 1.0 : 0.0) atol = 1.0e-13
        end
    end

    # Distorted (but affine) 2-cell mesh: nothing is axis-aligned, so the
    # tests are sensitive to the full 2×2 vertex-gradient blocks of the dof
    # transformation.
    nodes2 = [Node(Vec(0.0, 0.0)), Node(Vec(2.5, 0.3)), Node(Vec(0.4, 1.8)), Node(Vec(2.9, 2.4))]
    grid2 = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes2)
    dh2 = DofHandler(grid2)
    add!(dh2, :u, ip)
    close!(dh2)
    cd1, cd2 = celldofs(dh2, 1), celldofs(dh2, 2)

    @testset "dof distribution" begin
        @test ndofs(dh2) == 14
        @test cd1[4:6] == cd2[1:3] # shared dof triple at grid node 2
        @test cd1[7:9] == cd2[7:9] # shared dof triple at grid node 3
        @test length(unique(vcat(cd1, cd2))) == 14
    end

    # Exact interpolant of a full cubic: dofs are (u, ∂u/∂x, ∂u/∂y) (physical)
    # at the vertices plus the centroid value
    f_ex(x) = 1 + 2x[1] - x[2] + x[1]^2 * x[2] - 0.5x[2]^3 + x[1]^3
    u2 = zeros(ndofs(dh2))
    for cellid in 1:2
        cd = celldofs(dh2, cellid)
        coords = getcoordinates(grid2, cellid)
        for (v, x) in enumerate(coords)
            g = Ferrite.Tensors.gradient(f_ex, x)
            u2[cd[3 * (v - 1) + 1]] = f_ex(x)
            u2[cd[3 * (v - 1) + 2]] = g[1]
            u2[cd[3 * (v - 1) + 3]] = g[2]
        end
        u2[cd[10]] = f_ex(sum(coords) / 3)
    end

    @testset "mapping exactness on distorted mesh" begin
        cv = CellValues(QuadratureRule{RefTriangle}(4), ip; update_hessians = true)
        for cellid in 1:2
            coords = getcoordinates(grid2, cellid)
            reinit!(cv, coords)
            ue = u2[celldofs(dh2, cellid)]
            for qp in 1:getnquadpoints(cv)
                x = spatial_coordinate(cv, qp, coords)
                @test function_value(cv, qp, ue) ≈ f_ex(x) atol = 1.0e-11
                @test function_gradient(cv, qp, ue) ≈ Ferrite.Tensors.gradient(f_ex, x) atol = 1.0e-11
                @test function_hessian(cv, qp, ue) ≈ Ferrite.Tensors.hessian(f_ex, x) atol = 1.0e-10
            end
        end
    end

    @testset "FacetValues" begin
        fv = FacetValues(FacetQuadratureRule{RefTriangle}(3), ip)
        for cellid in 1:2, facetnr in 1:3
            coords = getcoordinates(grid2, cellid)
            reinit!(fv, coords, facetnr)
            ue = u2[celldofs(dh2, cellid)]
            for qp in 1:getnquadpoints(fv)
                x = spatial_coordinate(fv, qp, coords)
                @test function_value(fv, qp, ue) ≈ f_ex(x) atol = 1.0e-11
                @test function_gradient(fv, qp, ue) ≈ Ferrite.Tensors.gradient(f_ex, x) atol = 1.0e-11
            end
        end
    end

    @testset "continuity for arbitrary dof values" begin
        urand = rand(ndofs(dh2))
        # C0 along the shared edge (grid nodes 2 and 3): the physical point
        # p(t) = (1-t) n₂ + t n₃ has reference coordinates (0, 1-t) in cell 1
        # and (1-t, 0) in cell 2
        for t in (0.2, 0.5, 0.7)
            cv1 = CellValues(QuadratureRule{RefTriangle}([1.0], [Vec((0.0, 1 - t))]), ip)
            cv2 = CellValues(QuadratureRule{RefTriangle}([1.0], [Vec((1 - t, 0.0))]), ip)
            reinit!(cv1, getcoordinates(grid2, 1))
            reinit!(cv2, getcoordinates(grid2, 2))
            @test function_value(cv1, 1, urand[cd1]) ≈ function_value(cv2, 1, urand[cd2])
        end
        # C1 at the shared vertices: the full physical gradient agrees there
        # (but not along the rest of the edge)
        for (ξ1, ξ2) in ((Vec((0.0, 1.0)), Vec((1.0, 0.0))), (Vec((0.0, 0.0)), Vec((0.0, 0.0))))
            cv1 = CellValues(QuadratureRule{RefTriangle}([1.0], [ξ1]), ip)
            cv2 = CellValues(QuadratureRule{RefTriangle}([1.0], [ξ2]), ip)
            reinit!(cv1, getcoordinates(grid2, 1))
            reinit!(cv2, getcoordinates(grid2, 2))
            @test function_gradient(cv1, 1, urand[cd1]) ≈ function_gradient(cv2, 1, urand[cd2])
        end
    end

    @testset "evaluate_at_grid_nodes and PointEvalHandler" begin
        @test evaluate_at_grid_nodes(dh2, u2, :u) ≈ [f_ex(n.x) for n in nodes2]
        ph = PointEvalHandler(grid2, [Vec((0.9, 0.8)), Vec((2.2, 1.4))])
        @test evaluate_at_points(ph, dh2, u2, :u) ≈ [f_ex(Vec((0.9, 0.8))), f_ex(Vec((2.2, 1.4)))]
    end

    @testset "Poisson patch test with Dirichlet on all dof kinds" begin
        # The exact solution lies in the trial space and all boundary-vertex
        # dofs (value + gradient) are prescribed to it, which makes the trace
        # of every test function vanish on ∂Ω: the Galerkin solution of
        # -Δu = -Δf_ex reproduces f_ex exactly.
        grid = generate_grid(Triangle, (3, 2), Vec((0.0, 0.0)), Vec((1.5, 1.0)))
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        ∂Ω = union((getfacetset(grid, name) for name in ("left", "right", "top", "bottom"))...)
        ch = ConstraintHandler(dh)
        add!(ch, Dirichlet(:u, ∂Ω, x -> f_ex(x)))
        add!(ch, Dirichlet(:u, ∂Ω, x -> Ferrite.Tensors.gradient(f_ex, x)[1]; kind = :derivative_x))
        add!(ch, Dirichlet(:u, ∂Ω, x -> Ferrite.Tensors.gradient(f_ex, x)[2]; kind = :derivative_y))
        close!(ch)
        update!(ch, 0.0)
        # 4×3 vertex grid with 10 boundary vertices, three dofs each
        @test length(ch.prescribed_dofs) == 3 * 10
        Δf(x) = Ferrite.Tensors.laplace(f_ex, x)
        cv = CellValues(QuadratureRule{RefTriangle}(4), ip)
        K = allocate_matrix(dh)
        f = zeros(ndofs(dh))
        assembler = start_assemble(K, f)
        ndof = getnbasefunctions(ip)
        Ke = zeros(ndof, ndof)
        fe = zeros(ndof)
        for cell in CellIterator(dh)
            reinit!(cv, cell)
            fill!(Ke, 0.0)
            fill!(fe, 0.0)
            for qp in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv, qp)
                x = spatial_coordinate(cv, qp, getcoordinates(cell))
                for i in 1:ndof
                    gi = shape_gradient(cv, qp, i)
                    fe[i] += -Δf(x) * shape_value(cv, qp, i) * dΩ
                    for j in 1:ndof
                        Ke[i, j] += (gi ⋅ shape_gradient(cv, qp, j)) * dΩ
                    end
                end
            end
            assemble!(assembler, celldofs(cell), Ke, fe)
        end
        apply!(K, f, ch)
        u = K \ f
        ph = PointEvalHandler(grid, [Vec((0.35, 0.6)), Vec((1.2, 0.15))])
        @test evaluate_at_points(ph, dh, u, :u) ≈ [f_ex(Vec((0.35, 0.6))), f_ex(Vec((1.2, 0.15)))] atol = 1.0e-10
    end

    @testset "unsupported usage errors" begin
        @test_throws ArgumentError ip^2
        @test_throws ArgumentError CellValues(QuadratureRule{RefTriangle}(2), ip, Lagrange{RefTriangle, 2}()^2) # curved geometry
        @test_throws ArgumentError CellValues(QuadratureRule{RefTriangle}(2), ip, Lagrange{RefTriangle, 1}()^3) # embedded
        @test_throws ErrorException apply_analytical!(zeros(ndofs(dh2)), dh2, :u, x -> 0.0)
    end
end
