using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "EnrichedGalerkin" begin
    ips = (
        EnrichedGalerkin{RefLine, 1}(), EnrichedGalerkin{RefLine, 2}(),
        EnrichedGalerkin{RefTriangle, 1}(), EnrichedGalerkin{RefTriangle, 2}(), EnrichedGalerkin{RefTriangle, 3}(),
        EnrichedGalerkin{RefQuadrilateral, 1}(), EnrichedGalerkin{RefQuadrilateral, 2}(),
        EnrichedGalerkin{RefTetrahedron, 1}(), EnrichedGalerkin{RefTetrahedron, 2}(),
        EnrichedGalerkin{RefHexahedron, 1}(),
        EnrichedGalerkin{RefPrism, 1}(), EnrichedGalerkin{RefPrism, 2}(),
        EnrichedGalerkin{RefPyramid, 1}(), EnrichedGalerkin{RefPyramid, 2}(),
    )

    # (i) Interpolation-level tests. Not a partition of unity: the Lagrange
    # part sums to one and the enrichment constant adds one, so the sum is 2.
    for ip in ips
        test_interpolation_properties(ip)
        test_type_genericity(ip)
        @testset "sum of basis == 2: $ip" begin
            for _ in 1:5
                ξ = sample_reference_point(getrefshape(ip))
                @test sum(i -> reference_shape_value(ip, ξ, i), 1:getnbasefunctions(ip)) ≈ 2 atol = 1.0e-12
            end
        end
    end

    # P1 is contained in the span for every cell (higher degrees are covered
    # through the Lagrange delegation and the symfem cross-check).
    for ip in ips
        test_polynomial_reproduction(ip, 1; space = :P)
    end

    # Live symfem cross-check. symfem's enriched Galerkin DOFs are its
    # Lagrange DOFs followed by the cell constant, so the permutation is the
    # Ferrite <-> symfem Lagrange permutation (derived geometrically from the
    # Lagrange DOF points at runtime) with the constant appended. Excluded:
    # the prism and the pyramid, since symfem does not implement enriched
    # Galerkin on them ("element cannot be created on a prism/pyramid") --
    # and Ferrite's pyramid Lagrange is a polynomial variant of the rational
    # DefElement basis anyway. The delegation has no per-cell code, so the
    # remaining cells cover the implementation.
    function symfem_lagrange_perm(shape, degree)
        lag = Lagrange{shape, degree}()
        el = symfem.create_element(symfem_cellname(shape), "Lagrange", degree)
        dim = Ferrite.getrefdim(lag)
        spts = [Vec{dim}(ntuple(j -> pyconvert(Float64, pybuiltins.float(d.point[j - 1])), dim)) for d in el.dofs]
        fpts = [symfem_coords(shape, ξ) for ξ in Ferrite.reference_coordinates(lag)]
        return [findfirst(p -> norm(p - fp) < 1.0e-10, spts) - 1 for fp in fpts]
    end
    for ip in ips
        getrefshape(ip) in (RefPrism, RefPyramid) && continue
        shape, k = getrefshape(ip), Ferrite.getorder(ip)
        perm = vcat(symfem_lagrange_perm(shape, k), getnbasefunctions(ip) - 1)
        test_symfem_reference(ip, "enriched Galerkin", k, perm)
    end

    # (ii) Integration test: an affine function plus a cell constant is
    # represented exactly (Lagrange nodal coefficients + enrichment value).
    @testset "CellValues affine + constant: $ip" for ip in (
            EnrichedGalerkin{RefTriangle, 2}(), EnrichedGalerkin{RefQuadrilateral, 1}(),
        )
        shape = getrefshape(ip)
        coords = shape === RefTriangle ?
            [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))] :
            [Vec((0.0, 0.0)), Vec((2.0, 0.3)), Vec((2.4, 1.8)), Vec((0.4, 1.5))]
        lag = Ferrite.Lagrange{shape, Ferrite.getorder(ip)}()
        geo = Lagrange{shape, 1}()
        c = 0.75
        f(x) = 1 + 2x[1] - x[2] # affine part, interpolated through the Lagrange nodes
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:length(coords))
        ue = vcat([f(spatial(ξ)) for ξ in Ferrite.reference_coordinates(lag)], c)
        cv = CellValues(QuadratureRule{shape}(3), ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) + c atol = 1.0e-12
        end
    end

    # (iii) DofHandler: Lagrange DOFs shared as usual, one unshared constant
    # per cell (the last cell dof).
    @testset "two-cell sharing and jumps: $ip" for ip in (
            EnrichedGalerkin{RefTriangle, 2}(), EnrichedGalerkin{RefTriangle, 3}(),
        )
        k = Ferrite.getorder(ip)
        nodes = [
            Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
            Node(Vec((0.0, 1.0))), Node(Vec((1.0, 1.0))),
        ]
        # Shared edge: (2, 3) in cell 1, (3, 2) in cell 2 -> opposite orientation
        grid = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        # Continuous Lagrange on this mesh: 4 vertices + 5 edges * (k-1) + 2 * interior
        nlag_mesh = 4 + 5 * (k - 1) + 2 * ((k - 1) * (k - 2) ÷ 2)
        @test ndofs(dh) == nlag_mesh + 2
        dofs1, dofs2 = celldofs(dh, 1), celldofs(dh, 2)
        # Shared: the Lagrange DOFs of the common edge (2 vertices + k-1 interior)
        @test length(intersect(dofs1, dofs2)) == 2 + (k - 1)
        const1, const2 = dofs1[end], dofs2[end]
        @test const1 != const2 && const1 ∉ dofs2 && const2 ∉ dofs1
        N = getnbasefunctions(ip)
        eval_ref(ξ, dofs, u) = sum(u[dofs[i]] * reference_shape_value(ip, ξ, i) for i in 1:N)
        u = rand(ndofs(dh))
        for s in (0.2, 0.5, 0.85)
            # Physical point (1-s, s) on the shared edge, in each cell's reference coords
            ξ1, ξ2 = Vec((0.0, 1 - s)), Vec((1 - s, 0.0))
            # With equal constants the field is continuous across the edge...
            u[const2] = u[const1]
            @test eval_ref(ξ1, dofs1, u) ≈ eval_ref(ξ2, dofs2, u) atol = 1.0e-12
            # ...and differing constants jump by exactly their difference.
            u[const2] = u[const1] + 3.0
            @test eval_ref(ξ2, dofs2, u) - eval_ref(ξ1, dofs1, u) ≈ 3.0 atol = 1.0e-12
        end
    end
end
