using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "Taylor" begin
    ips = (
        Taylor{RefLine, 1}(), Taylor{RefLine, 2}(),
        Taylor{RefTriangle, 1}(), Taylor{RefTriangle, 2}(),
    )

    # (i) Interpolation-level tests (not nodal, no partition of unity; the
    # span is the full P_k).
    for ip in ips
        test_interpolation_properties(ip)
        test_polynomial_reproduction(ip, Ferrite.getorder(ip); space = :P)
        test_type_genericity(ip)
    end

    # Live symfem cross-check ("Taylor"): scaled permutations. The interval
    # scalings 1/2, 2, 4 are the reference-map factors for the integral and
    # derivative DOFs of x = (ξ+1)/2; on the triangle (same reference cell)
    # only symfem's lexicographic derivative order differs.
    test_symfem_reference(ips[1], "Taylor", 1, [0, 1]; scales = [1 / 2, 2.0])
    test_symfem_reference(ips[2], "Taylor", 2, [0, 1, 2]; scales = [1 / 2, 2.0, 4.0])
    test_symfem_reference(ips[3], "Taylor", 1, [0, 2, 1])
    test_symfem_reference(ips[4], "Taylor", 2, [0, 3, 1, 5, 4, 2])

    # (ii) The defining property of the TaylorMapping: on a distorted affine
    # cell, the PHYSICAL DOFs (cell integral + midpoint derivatives of the
    # physical function) are dual to the mapped basis, so setting the
    # coefficients to those DOF values of a degree-k polynomial reproduces it
    # exactly.
    @testset "Kirby mapping duality: $ip" for ip in ips
        shape = getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        k = Ferrite.getorder(ip)
        coords = dim == 1 ? [Vec((-0.3,)), Vec((1.7,))] :
            [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))]
        geo = Lagrange{shape, 1}()
        b = Vec{dim}(ntuple(i -> Float64(i), dim))
        f(x) = k == 1 ? 1 + x ⋅ b : 1 + x ⋅ b + (x ⋅ b)^2 / 4
        # Physical DOFs of f
        cv = CellValues(QuadratureRule{shape}(4), ip, geo)
        reinit!(cv, coords)
        vol_int = sum(f(spatial_coordinate(cv, qp, coords)) * getdetJdV(cv, qp) for qp in 1:getnquadpoints(cv))
        mid = sum(coords) / length(coords) # affine: midpoint = mapped reference centroid
        g = Ferrite.Tensors.gradient(f, mid)
        ue = if dim == 1
            k == 1 ? [vol_int, g[1]] : [vol_int, g[1], Ferrite.Tensors.hessian(f, mid)[1, 1]]
        else
            H = Ferrite.Tensors.hessian(f, mid)
            k == 1 ? [vol_int, g[1], g[2]] : [vol_int, g[1], g[2], H[1, 1], H[1, 2], H[2, 2]]
        end
        for qp in 1:getnquadpoints(cv)
            x = spatial_coordinate(cv, qp, coords)
            @test function_value(cv, qp, ue) ≈ f(x) atol = 1.0e-11
            @test function_gradient(cv, qp, ue) ≈ Ferrite.Tensors.gradient(f, x) atol = 1.0e-10
        end
    end

    # (iii) DofHandler: all DOFs are cell DOFs, nothing is shared.
    @testset "two-cell independence" begin
        nodes = [
            Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
            Node(Vec((0.0, 1.0))), Node(Vec((1.0, 1.0))),
        ]
        grid = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes)
        for ip in (ips[3], ips[4])
            dh = DofHandler(grid)
            add!(dh, :u, ip)
            close!(dh)
            @test ndofs(dh) == 2 * getnbasefunctions(ip)
            @test isempty(intersect(celldofs(dh, 1), celldofs(dh, 2)))
        end
    end
end
