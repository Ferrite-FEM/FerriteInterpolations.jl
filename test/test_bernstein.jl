@testset "Bernstein" begin
    line_ips = (Bernstein{RefLine, 1}(), Bernstein{RefLine, 2}(), Bernstein{RefLine, 3}())
    tri_ips = (Bernstein{RefTriangle, 1}(), Bernstein{RefTriangle, 2}(), Bernstein{RefTriangle, 3}())
    tet_ips = (Bernstein{RefTetrahedron, 1}(), Bernstein{RefTetrahedron, 2}())

    # (i) Interpolation-level tests
    for ip in (line_ips..., tri_ips..., tet_ips...)
        test_interpolation_properties(ip)
        test_partition_of_unity(ip)
        test_polynomial_reproduction(ip, Ferrite.getorder(ip); space = :P)
        test_type_genericity(ip)
    end

    # Reference tabulation generated with symfem 2025.8.0 (dev/generate_bernstein_table.py).
    # Points are in Ferrite reference coordinates and values in Ferrite DOF order;
    # the generator maps coordinates (interval [0,1] -> [-1,1]) and permutes
    # symfem's DOF order (differing vertex/edge numbering and edge orientation)
    # structurally via symfem's entity_dofs.
    reference_tables = [
        Bernstein{RefLine, 1}() => [
            (Vec((-0.8,)), [0.9, 0.1]),
            (Vec((-0.2,)), [0.6, 0.4]),
            (Vec((0.5,)), [0.25, 0.75]),
            (Vec((0.9,)), [0.05, 0.95]),
        ],
        Bernstein{RefLine, 2}() => [
            (Vec((-0.8,)), [0.81, 0.01, 0.18]),
            (Vec((-0.2,)), [0.36, 0.16, 0.48]),
            (Vec((0.5,)), [0.0625, 0.5625, 0.375]),
            (Vec((0.9,)), [0.0025, 0.9025, 0.095]),
        ],
        Bernstein{RefLine, 3}() => [
            (Vec((-0.8,)), [0.729, 0.001, 0.243, 0.027]),
            (Vec((-0.2,)), [0.216, 0.064, 0.432, 0.288]),
            (Vec((0.5,)), [0.015625, 0.421875, 0.140625, 0.421875]),
            (Vec((0.9,)), [0.000125, 0.857375, 0.007125, 0.135375]),
        ],
        Bernstein{RefTriangle, 1}() => [
            (Vec((0.3, 0.2)), [0.3, 0.2, 0.5]),
            (Vec((0.1, 0.6)), [0.1, 0.6, 0.3]),
            (Vec((0.25, 0.25)), [0.25, 0.25, 0.5]),
            (Vec((0.05, 0.9)), [0.05, 0.9, 0.05]),
        ],
        Bernstein{RefTriangle, 2}() => [
            (Vec((0.3, 0.2)), [0.09, 0.04, 0.25, 0.12, 0.2, 0.3]),
            (Vec((0.1, 0.6)), [0.01, 0.36, 0.09, 0.12, 0.36, 0.06]),
            (Vec((0.25, 0.25)), [0.0625, 0.0625, 0.25, 0.125, 0.25, 0.25]),
            (Vec((0.05, 0.9)), [0.0025, 0.81, 0.0025, 0.09, 0.09, 0.005]),
        ],
        Bernstein{RefTriangle, 3}() => [
            (Vec((0.3, 0.2)), [0.027, 0.008, 0.125, 0.054, 0.036, 0.06, 0.15, 0.225, 0.135, 0.18]),
            (Vec((0.1, 0.6)), [0.001, 0.216, 0.027, 0.018, 0.108, 0.324, 0.162, 0.027, 0.009, 0.108]),
            (Vec((0.25, 0.25)), [0.015625, 0.015625, 0.125, 0.046875, 0.046875, 0.09375, 0.1875, 0.1875, 0.09375, 0.1875]),
            (Vec((0.05, 0.9)), [0.000125, 0.729, 0.000125, 0.00675, 0.1215, 0.1215, 0.00675, 0.000375, 0.000375, 0.0135]),
        ],
        Bernstein{RefTetrahedron, 1}() => [
            (Vec((0.2, 0.3, 0.4)), [0.1, 0.2, 0.3, 0.4]),
            (Vec((0.1, 0.1, 0.1)), [0.7, 0.1, 0.1, 0.1]),
            (Vec((0.5, 0.2, 0.25)), [0.05, 0.5, 0.2, 0.25]),
            (Vec((0.05, 0.6, 0.3)), [0.05, 0.05, 0.6, 0.3]),
        ],
        Bernstein{RefTetrahedron, 2}() => [
            (Vec((0.2, 0.3, 0.4)), [0.01, 0.04, 0.09, 0.16, 0.04, 0.12, 0.06, 0.08, 0.16, 0.24]),
            (Vec((0.1, 0.1, 0.1)), [0.49, 0.01, 0.01, 0.01, 0.14, 0.02, 0.14, 0.14, 0.02, 0.02]),
            (Vec((0.5, 0.2, 0.25)), [0.0025, 0.25, 0.04, 0.0625, 0.05, 0.2, 0.02, 0.025, 0.25, 0.1]),
            (Vec((0.05, 0.6, 0.3)), [0.0025, 0.0025, 0.36, 0.09, 0.005, 0.06, 0.06, 0.03, 0.03, 0.36]),
        ],
    ]
    @testset "symfem reference values: $ip" for (ip, table) in reference_tables
        for (ξ, expected) in table
            vals = [Ferrite.reference_shape_value(ip, ξ, i) for i in 1:Ferrite.getnbasefunctions(ip)]
            @test vals ≈ expected atol = 1.0e-13
        end
    end

    # Greville point of DOF i: the Bernstein basis has linear precision,
    # f = Σᵢ f(gᵢ) Nᵢ exactly for affine f, with gᵢ = Σⱼ αⱼ vⱼ / order.
    function greville_points(ip::Bernstein{shape, order}) where {shape, order}
        verts = Ferrite.reference_coordinates(Lagrange{shape, 1}())
        return [sum(α .* verts) / order for (_, α) in FerriteInterpolations._bernstein_data(ip)]
    end

    # (ii) Integration tests: CellValues on a single distorted (affine) cell,
    # checking exact reproduction of an affine function via Greville coefficients.
    @testset "CellValues affine reproduction: $ip" for (ip, coords, qr) in (
            (Bernstein{RefLine, 3}(), [Vec((-0.3,)), Vec((1.7,))], QuadratureRule{RefLine}(3)),
            (Bernstein{RefTriangle, 2}(), [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))], QuadratureRule{RefTriangle}(3)),
            (Bernstein{RefTriangle, 3}(), [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))], QuadratureRule{RefTriangle}(3)),
            (Bernstein{RefTetrahedron, 2}(), [Vec((0.0, 0.0, 0.0)), Vec((1.2, 0.1, 0.2)), Vec((0.1, 1.4, 0.3)), Vec((0.2, 0.3, 1.1))], QuadratureRule{RefTetrahedron}(3)),
        )
        shape = Ferrite.getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        geo = Lagrange{shape, 1}()
        b = Vec{dim}(ntuple(i -> Float64(i), dim))
        f(x) = 2 + x ⋅ b # affine
        # Affine geometry: spatial Greville points are the map of the reference ones
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:length(coords))
        ue = [f(spatial(g)) for g in greville_points(ip)]
        cv = CellValues(qr, ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) atol = 1.0e-12
            @test function_gradient(cv, qp, ue) ≈ b atol = 1.0e-12
        end
    end

    # FacetValues: partition of unity (all-ones dof vector interpolates 1 on facets).
    @testset "FacetValues partition of unity: $ip" for ip in (tri_ips[2], tri_ips[3], tet_ips[2])
        shape = Ferrite.getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        geo = Lagrange{shape, 1}()
        coords = dim == 2 ?
            [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))] :
            [Vec((0.0, 0.0, 0.0)), Vec((1.2, 0.1, 0.2)), Vec((0.1, 1.4, 0.3)), Vec((0.2, 0.3, 1.1))]
        fqr = FacetQuadratureRule{shape}(2)
        fv = FacetValues(fqr, ip, geo)
        ue = ones(Ferrite.getnbasefunctions(ip))
        for facet in 1:Ferrite.nfacets(ip)
            reinit!(fv, coords, facet)
            for qp in 1:getnquadpoints(fv)
                @test function_value(fv, qp, ue) ≈ 1.0 atol = 1.0e-13
            end
        end
    end

    # Two-cell continuity: two triangles whose shared edge has opposite local
    # orientation; the interpolated function must agree along the edge from
    # both sides (exercises dof distribution and, for order 3, the
    # adjust_dofs_during_distribution edge reversal).
    @testset "two-cell continuity: $ip" for ip in (tri_ips[2], tri_ips[3])
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
        u = rand(ndofs(dh))
        dofs1 = celldofs(dh, 1)
        dofs2 = celldofs(dh, 2)
        eval_ref(ξ, dofs) = sum(u[dofs[i]] * Ferrite.reference_shape_value(ip, ξ, i) for i in 1:Ferrite.getnbasefunctions(ip))
        for s in (0.15, 0.4, 0.5, 0.8)
            # Physical point (1-s, s): cell 1 edge (2,3) param -> ξ = (0, 1-s),
            # cell 2 (local vertices n2, n4, n3) edge (3,1) param -> ξ = (1-s, 0)
            v1 = eval_ref(Vec((0.0, 1 - s)), dofs1)
            v2 = eval_ref(Vec((1 - s, 0.0)), dofs2)
            @test v1 ≈ v2 atol = 1.0e-13
        end
    end
end
