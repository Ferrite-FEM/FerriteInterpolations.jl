using FerriteInterpolations
using Ferrite
using Test
using LinearAlgebra: I
include("test_utils.jl")

@testset "GaussLegendre" begin
    ips = (
        GaussLegendre{RefLine, 1}(), GaussLegendre{RefLine, 2}(), GaussLegendre{RefLine, 3}(),
        GaussLegendre{RefQuadrilateral, 1}(), GaussLegendre{RefQuadrilateral, 2}(), GaussLegendre{RefQuadrilateral, 3}(),
        GaussLegendre{RefHexahedron, 1}(), GaussLegendre{RefHexahedron, 2}(),
    )

    # (i) Interpolation-level tests. Modal element: no partition of unity, no
    # Kronecker property. The span is the full Q_k space.
    for ip in ips
        test_interpolation_properties(ip)
        test_polynomial_reproduction(ip, Ferrite.getorder(ip); space = :Q)
        test_type_genericity(ip)
    end

    # Live symfem cross-check ("Lagrange variant=legendre"); all DOFs are
    # cell DOFs in symfem's order, so the permutations are trivial.
    for ip in (ips[1], ips[2], ips[3], ips[4], ips[5], ips[7])
        N = getnbasefunctions(ip)
        test_symfem_reference(ip, "Lagrange", Ferrite.getorder(ip), collect(0:(N - 1)); variant = "legendre")
    end

    # The defining property: the basis is orthonormal with respect to the
    # normalized reference measure, i.e. the reference mass matrix is
    # 2^dim * I. Computed with CellValues on the reference-shaped cell.
    refcoords(::Type{RefLine}) = [Vec((-1.0,)), Vec((1.0,))]
    refcoords(::Type{RefQuadrilateral}) = [Vec((-1.0, -1.0)), Vec((1.0, -1.0)), Vec((1.0, 1.0)), Vec((-1.0, 1.0))]
    function refcoords(::Type{RefHexahedron})
        return [
            Vec((-1.0, -1.0, -1.0)), Vec((1.0, -1.0, -1.0)), Vec((1.0, 1.0, -1.0)), Vec((-1.0, 1.0, -1.0)),
            Vec((-1.0, -1.0, 1.0)), Vec((1.0, -1.0, 1.0)), Vec((1.0, 1.0, 1.0)), Vec((-1.0, 1.0, 1.0)),
        ]
    end
    @testset "orthonormality: $ip" for ip in ips
        shape = getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        N = getnbasefunctions(ip)
        cv = CellValues(QuadratureRule{shape}(Ferrite.getorder(ip) + 1), ip, Lagrange{shape, 1}())
        reinit!(cv, refcoords(shape))
        M = [
            sum(shape_value(cv, qp, i) * shape_value(cv, qp, j) * getdetJdV(cv, qp) for qp in 1:getnquadpoints(cv))
                for i in 1:N, j in 1:N
        ]
        @test M ≈ 2^dim * I atol = 1.0e-11
    end

    # (ii) DofHandler: all DOFs are cell DOFs, nothing is shared.
    @testset "two-cell independence" begin
        nodes = [Node(Vec((x, y))) for y in (0.0, 1.0) for x in (0.0, 1.0, 2.0)]
        grid = Grid([Quadrilateral((1, 2, 5, 4)), Quadrilateral((2, 3, 6, 5))], nodes)
        for ip in (GaussLegendre{RefQuadrilateral, 1}(), GaussLegendre{RefQuadrilateral, 2}())
            dh = DofHandler(grid)
            add!(dh, :u, ip)
            close!(dh)
            @test ndofs(dh) == 2 * getnbasefunctions(ip)
            @test isempty(intersect(celldofs(dh, 1), celldofs(dh, 2)))
        end
    end
end
