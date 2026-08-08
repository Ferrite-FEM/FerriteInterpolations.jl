using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "Transition" begin
    combos = ((1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2))
    ips = map(eo -> Transition{RefTriangle, 2, eo}(), combos)

    # (i) Interpolation-level tests. P1 is always contained; the trace on an
    # order-1 edge must be linear.
    for ip in ips
        test_interpolation_properties(ip)
        test_partition_of_unity(ip)
        test_polynomial_reproduction(ip, 1; space = :P)
        test_type_genericity(ip)
        test_kronecker_delta(ip)
    end
    edge_param(::Val{1}, s) = Vec((1 - s, s))
    edge_param(::Val{2}, s) = Vec((0.0, 1 - s))
    edge_param(::Val{3}, s) = Vec((s, 0.0))
    @testset "order-1 edge traces are linear: $eo" for (eo, ip) in zip(combos, ips)
        for e in 1:3
            eo[e] == 1 || continue
            for i in 1:getnbasefunctions(ip)
                f(s) = reference_shape_value(ip, edge_param(Val(e), s), i)
                # linear <=> midpoint value is the average of the endpoints
                @test f(0.5) ≈ (f(0.0) + f(1.0)) / 2 atol = 1.0e-12
                @test f(0.25) ≈ (3f(0.0) + f(1.0)) / 4 atol = 1.0e-12
            end
        end
    end

    # Live symfem cross-check (edge orders permuted to symfem's numbering,
    # perms recorded in src/transition.jl).
    for (ip, perm) in (
            (ips[1], [1, 2, 0, 3]), (ips[2], [1, 2, 0, 3]), (ips[3], [1, 2, 0, 4, 3]),
            (ips[4], [1, 2, 0, 3]), (ips[5], [1, 2, 0, 4, 3]), (ips[6], [1, 2, 0, 4, 3]),
            (ips[7], [1, 2, 0, 5, 4, 3]),
        )
        eo = combos[findfirst(==(ip), ips)]
        test_symfem_reference(ip, "transition", 2, perm; edge_orders = pylist([eo[3], eo[2], eo[1]]))
    end

    # (ii) Two-cell conforming coupling: shared edge is Ferrite edge 2 of
    # cell 1 and edge 3 of cell 2 (opposite orientation). Both an order-2/2
    # and an order-1/1 pairing must be continuous along the WHOLE edge -- for
    # 1/1 because both traces are linear through the shared vertex DOFs.
    nodes = [
        Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
        Node(Vec((0.0, 1.0))), Node(Vec((1.0, 1.0))),
    ]
    grid = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes)
    @testset "two-cell continuity: shared order $shared_order" for (shared_order, eo1, eo2, nshared) in (
            (2, (1, 2, 1), (1, 1, 2), 3), # shared edge quadratic on both sides
            (1, (2, 1, 2), (2, 2, 1), 2), # shared edge linear on both sides
        )
        ip1 = Transition{RefTriangle, 2, eo1}()
        ip2 = Transition{RefTriangle, 2, eo2}()
        dh = DofHandler(grid)
        sdh1 = SubDofHandler(dh, Set([1]))
        add!(sdh1, :u, ip1)
        sdh2 = SubDofHandler(dh, Set([2]))
        add!(sdh2, :u, ip2)
        close!(dh)
        dofs1, dofs2 = celldofs(dh, 1), celldofs(dh, 2)
        @test length(intersect(dofs1, dofs2)) == nshared
        u = rand(ndofs(dh))
        eval_ref(ip, ξ, dofs) = sum(u[dofs[i]] * reference_shape_value(ip, ξ, i) for i in 1:getnbasefunctions(ip))
        for s in (0.15, 0.4, 0.5, 0.8)
            v1 = eval_ref(ip1, Vec((0.0, 1 - s)), dofs1)
            v2 = eval_ref(ip2, Vec((1 - s, 0.0)), dofs2)
            @test v1 ≈ v2 atol = 1.0e-13
        end
    end
end
