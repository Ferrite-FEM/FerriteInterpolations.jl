using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "Bubble" begin
    ips = (
        Bubble{RefLine, 2}(), Bubble{RefLine, 3}(),
        Bubble{RefTriangle, 3}(), Bubble{RefTriangle, 4}(),
        Bubble{RefTetrahedron, 4}(),
    )

    # (i) Interpolation-level tests (no partition of unity: bubbles do not
    # sum to one, and the span is bubble * P_j, not a full P_k).
    for ip in ips
        test_interpolation_properties(ip)
        test_type_genericity(ip)
        test_kronecker_delta(ip) # all DOFs are point evaluations
    end

    # The defining property: every basis function vanishes on the whole cell
    # boundary.
    boundary_points(::Type{RefLine}) = [Vec((-1.0,)), Vec((1.0,))]
    function boundary_points(::Type{RefTriangle})
        # Vertices (1,0), (0,1), (0,0); points on all three edges
        return [Vec((1 - t, t)) for t in 0:0.25:1] ∪
            [Vec((0.0, t)) for t in 0:0.25:1] ∪
            [Vec((t, 0.0)) for t in 0:0.25:1]
    end
    function boundary_points(::Type{RefTetrahedron})
        # Points on all four faces of the reference tetrahedron
        ts = ((0.2, 0.3), (0.1, 0.6), (0.5, 0.25))
        return vcat(
            [Vec((a, b, 0.0)) for (a, b) in ts],       # z = 0
            [Vec((a, 0.0, b)) for (a, b) in ts],       # y = 0
            [Vec((0.0, a, b)) for (a, b) in ts],       # x = 0
            [Vec((a, b, 1 - a - b)) for (a, b) in ts], # x + y + z = 1
        )
    end
    @testset "vanishes on boundary: $ip" for ip in ips
        for ξ in boundary_points(getrefshape(ip)), i in 1:getnbasefunctions(ip)
            @test reference_shape_value(ip, ξ, i) ≈ 0 atol = 1.0e-13
        end
    end

    # Live symfem cross-check. Interior DOFs keep symfem's order and the
    # simplex coordinates agree (the RefLine map is handled by the helper),
    # so the permutations are trivial.
    for (ip, perm) in (
            (Bubble{RefLine, 2}(), [0]),
            (Bubble{RefLine, 3}(), [0, 1]),
            (Bubble{RefTriangle, 3}(), [0]),
            (Bubble{RefTriangle, 4}(), [0, 1, 2]),
            (Bubble{RefTetrahedron, 4}(), [0]),
        )
        test_symfem_reference(ip, "bubble", Ferrite.getorder(ip), perm)
    end

    # (ii) Integration tests: exact integrals of the single-bubble elements
    # over the reference cell (identity geometry mapping), computed with
    # CellValues. For the unit simplex, int x^a y^b (1-x-y)^c = a!b!c!/(a+b+c+2)!.
    @testset "CellValues bubble integral: $shape" for (shape, ip, coords, qr, exact) in (
            (RefLine, Bubble{RefLine, 2}(), [Vec((-1.0,)), Vec((1.0,))], QuadratureRule{RefLine}(3), 4 / 3),
            (RefTriangle, Bubble{RefTriangle, 3}(), [Vec((1.0, 0.0)), Vec((0.0, 1.0)), Vec((0.0, 0.0))], QuadratureRule{RefTriangle}(4), 27 / 120),
            (RefTetrahedron, Bubble{RefTetrahedron, 4}(), [Vec((0.0, 0.0, 0.0)), Vec((1.0, 0.0, 0.0)), Vec((0.0, 1.0, 0.0)), Vec((0.0, 0.0, 1.0))], QuadratureRule{RefTetrahedron}(5), 256 / 5040),
        )
        cv = CellValues(qr, ip, Lagrange{shape, 1}())
        reinit!(cv, coords)
        integral = sum(function_value(cv, qp, [1.0]) * getdetJdV(cv, qp) for qp in 1:getnquadpoints(cv))
        @test integral ≈ exact atol = 1.0e-12
    end

    # FacetValues: the interpolated field vanishes on every facet whatever
    # the coefficients are.
    @testset "FacetValues boundary zero: $ip" for ip in (Bubble{RefTriangle, 3}(), Bubble{RefTriangle, 4}(), Bubble{RefTetrahedron, 4}())
        shape = getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        coords = dim == 2 ?
            [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))] :
            [Vec((0.0, 0.0, 0.0)), Vec((1.2, 0.1, 0.2)), Vec((0.1, 1.4, 0.3)), Vec((0.2, 0.3, 1.1))]
        fv = FacetValues(FacetQuadratureRule{shape}(3), ip, Lagrange{shape, 1}())
        ue = rand(getnbasefunctions(ip))
        for facet in 1:Ferrite.nfacets(ip)
            reinit!(fv, coords, facet)
            for qp in 1:getnquadpoints(fv)
                @test function_value(fv, qp, ue) ≈ 0 atol = 1.0e-13
            end
        end
    end

    # (iii) DofHandler: all DOFs are cell-interior, nothing is shared.
    @testset "two-cell independence: $ip" for ip in (Bubble{RefTriangle, 3}(), Bubble{RefTriangle, 4}())
        nodes = [
            Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
            Node(Vec((0.0, 1.0))), Node(Vec((1.0, 1.0))),
        ]
        grid = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        @test ndofs(dh) == 2 * getnbasefunctions(ip)
        @test isempty(intersect(celldofs(dh, 1), celldofs(dh, 2)))
    end
end
