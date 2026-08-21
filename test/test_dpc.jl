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

    # (ii) Integration tests: dPc uses the L2 Piola mapping N(x) = N̂(ξ)/detJ
    # (src/l2_piola.jl). On affinely-mapped cells (parallelogram /
    # parallelepiped) detJ is constant, so nodal interpolation of a degree-k
    # polynomial is exact with the coefficients scaled by detJ.
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
        detJ = Ferrite.Tensors.det(Ferrite.Tensors.gradient(spatial, zero(Vec{dim})))
        ue = [detJ * f(spatial(ξ)) for ξ in Ferrite.reference_coordinates(ip)]
        cv = CellValues(QuadratureRule{shape}(k + 1), ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) atol = 1.0e-11
        end
    end

    # L2 Piola mapping on a genuinely NON-affine quadrilateral: check values
    # and gradients from CellValues against an independent computation of
    # N̂(ξ)/detJ(ξ) and its chain-rule gradient via nested automatic
    # differentiation of the bilinear geometry map (this validates the
    # ∇ₓ(detJ) term in the mapping, which vanishes on affine cells).
    @testset "L2 Piola on non-affine cell: $ip" for ip in (DPC{RefQuadrilateral, 2}(), DPC{RefQuadrilateral, 3}())
        coords = [Vec((0.0, 0.0)), Vec((2.1, 0.2)), Vec((1.7, 2.4)), Vec((-0.3, 1.1))]
        geo = Lagrange{RefQuadrilateral, 1}()
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:4)
        jac(ξ) = Ferrite.Tensors.gradient(spatial, ξ)
        manual_value(ξ, i) = reference_shape_value(ip, ξ, i) / Ferrite.Tensors.det(jac(ξ))
        manual_grad(ξ, i) = inv(jac(ξ))' ⋅ Ferrite.Tensors.gradient(η -> manual_value(η, i), ξ)
        qr = QuadratureRule{RefQuadrilateral}(3)
        cv = CellValues(qr, ip, geo)
        reinit!(cv, coords)
        for (qp, ξ) in pairs(Ferrite.getpoints(qr)), i in 1:getnbasefunctions(ip)
            @test shape_value(cv, qp, i) ≈ manual_value(ξ, i) atol = 1.0e-12
            @test shape_gradient(cv, qp, i) ≈ manual_grad(ξ, i) atol = 1.0e-11
        end
        # FacetValues use the same mapping: with Σᵢ N̂ᵢ = 1 (reference
        # partition of unity) the mapped sum must equal 1/detJ on facets too.
        fqr = FacetQuadratureRule{RefQuadrilateral}(2)
        fv = FacetValues(fqr, ip, geo)
        for facet in 1:4
            reinit!(fv, coords, facet)
            for (qp, ξ) in pairs(Ferrite.getpoints(fqr, facet))
                s = sum(shape_value(fv, qp, i) for i in 1:getnbasefunctions(ip))
                @test s ≈ 1 / Ferrite.Tensors.det(jac(ξ)) atol = 1.0e-12
            end
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

    # Vector dPc is DPC^vdim via Ferrite's VectorizedInterpolation (with the
    # identity mapping, matching DefElement's vector dPc declaration). BROKEN
    # on this branch: the pinned Ferrite branch (Ferrite-FEM/Ferrite.jl#1391)
    # added a `VectorizedInterpolation` constructor guard that refuses to
    # vectorize any scalar interpolation with a non-identity mapping -- correct
    # for elements whose vectorization must not drop the mapping/transformation
    # (e.g. Hermite), but too strict for dPc where identity-mapped vectorization
    # is the intended element. Needs an upstream opt-out (e.g. a separate trait
    # for the guard instead of `physical_basis_is_reference_basis`) before that
    # PR merges; the full testset is in the git history of this branch.
    @testset "vector dPc: DPC^vdim" begin
        @test_broken DPC{RefQuadrilateral, 2}()^2 isa Ferrite.VectorizedInterpolation
    end
end
