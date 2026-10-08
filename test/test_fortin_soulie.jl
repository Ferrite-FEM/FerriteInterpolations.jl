using FerriteInterpolations
using Ferrite
using LinearAlgebra: dot, norm, rank
using Test
include("test_utils.jl")

# The two Gauss-Legendre points of the edge from a to b
gauss_points(a, b) = (g = 1 / 2 - sqrt(3) / 6; (a + g * (b - a), a + (1 - g) * (b - a)))

@testset "FortinSoulie" begin
    ip = FortinSoulie{RefTriangle, 2}()
    lag = Lagrange{RefTriangle, 2}()
    verts = Ferrite.reference_coordinates(Lagrange{RefTriangle, 1}())

    # (i) Interpolation-level tests
    test_interpolation_properties(ip)
    test_type_genericity(ip)
    test_polynomial_reproduction(ip, 2; space = :P)

    @testset "Lagrange part and neutral function" begin
        for _ in 1:20
            ξ = sample_reference_point(RefTriangle)
            for i in 1:6
                @test reference_shape_value(ip, ξ, i) ≈ reference_shape_value(lag, ξ, i)
            end
            # The neutral function lies in P2: -1 at the vertices, 1/2 at the edge
            # midpoints (this is the one local linear dependency)
            φ₀ = sum(-reference_shape_value(lag, ξ, i) for i in 1:3) +
                sum(reference_shape_value(lag, ξ, i) / 2 for i in 4:6)
            @test reference_shape_value(ip, ξ, 7) ≈ φ₀ atol = 1.0e-12
        end
        # ...vanishes at the six Gauss-Legendre points and is 1 at the centroid
        for (a, b) in Ferrite.reference_edges(RefTriangle), p in gauss_points(verts[a], verts[b])
            @test reference_shape_value(ip, p, 7) ≈ 0 atol = 1.0e-14
        end
        @test reference_shape_value(ip, Vec((1 / 3, 1 / 3)), 7) ≈ 1
        # Seven functions spanning the six-dimensional P2
        pts = [sample_reference_point(RefTriangle) for _ in 1:30]
        @test rank([reference_shape_value(ip, ξ, i) for ξ in pts, i in 1:7]) == 6
    end

    # (ii) Integration test: continuous P2 part plus a neutral function on a
    # physical cell (the neutral function is mapped with the identity mapping)
    @testset "CellValues P2 + neutral function" begin
        coords = [Vec((0.0, 0.0)), Vec((2.5, 0.3)), Vec((0.4, 1.8))]
        geo = Lagrange{RefTriangle, 1}()
        f(x) = 1 + 2x[1] - x[2] + x[1]^2 / 2 + x[1] * x[2] - 3x[2]^2 / 10
        spatial(ξ) = sum(Ferrite.reference_shape_value(geo, ξ, j) * coords[j] for j in 1:3)
        c = 0.7
        ue = [[f(spatial(ξ)) for ξ in Ferrite.reference_coordinates(lag)]; c]
        cv = CellValues(QuadratureRule{RefTriangle}(4), ip, geo)
        reinit!(cv, coords)
        for qp in 1:getnquadpoints(cv)
            ξ = Ferrite.getpoints(QuadratureRule{RefTriangle}(4))[qp]
            @test function_value(cv, qp, ue) ≈ f(spatial_coordinate(cv, qp, coords)) + c * reference_shape_value(ip, ξ, 7)
        end
    end

    # (iii) DofHandler: the Lagrange DOFs are shared, the neutral functions are
    # not, and the field is continuous exactly at the Gauss-Legendre points
    @testset "two-cell Gauss-Legendre continuity" begin
        nodes = [
            Node(Vec((0.0, 0.0))), Node(Vec((1.0, 0.0))),
            Node(Vec((0.0, 1.0))), Node(Vec((1.1, 0.9))),
        ]
        # Shared edge (2, 3), with opposite orientation in the two cells
        grid = Grid([Triangle((1, 2, 3)), Triangle((2, 4, 3))], nodes)
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        @test ndofs(dh) == 4 + 5 + 2
        dofs1, dofs2 = celldofs(dh, 1), celldofs(dh, 2)
        @test length(intersect(dofs1, dofs2)) == 3 # two vertices and the midpoint
        @test dofs1[7] ∉ dofs2 && dofs2[7] ∉ dofs1

        # Evaluate the field at a physical point of a cell (affine map inverted exactly)
        function value(cell, dofs, u, x)
            X = [nodes[n].x for n in grid.cells[cell].nodes]
            J = hcat(X[1] - X[3], X[2] - X[3])
            ξ = Vec{2}(Tuple(J \ (x - X[3])))
            return sum(u[dofs[i]] * reference_shape_value(ip, ξ, i) for i in 1:7)
        end
        a, b = nodes[2].x, nodes[3].x
        for _ in 1:5
            u = rand(ndofs(dh))
            for p in gauss_points(a, b)
                @test value(1, dofs1, u, p) ≈ value(2, dofs2, u, p) atol = 1.0e-12
            end
            # Elsewhere on the edge the traces differ by the neutral functions
            p = a + 0.3 * (b - a)
            @test !isapprox(value(1, dofs1, u, p), value(2, dofs2, u, p); atol = 1.0e-6)
        end
    end

    # (iv) The global dependency: all neutral functions with coefficient 1, all
    # vertex DOFs 1 and all midpoint DOFs -1/2 represent the zero function
    @testset "global dependency (paper, Proposition 1)" begin
        grid = generate_grid(Triangle, (3, 3))
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        u = zeros(ndofs(dh))
        for cell in CellIterator(dh)
            d = celldofs(cell)
            u[d[1:3]] .= 1
            u[d[4:6]] .= -1 / 2
            u[d[7]] = 1
        end
        for cell in CellIterator(dh), _ in 1:5
            ξ = sample_reference_point(RefTriangle)
            d = celldofs(cell)
            @test sum(u[d[i]] * reference_shape_value(ip, ξ, i) for i in 1:7) ≈ 0 atol = 1.0e-12
        end
        # On this simply connected mesh, the global space has dimension
        # 2 * (number of edges), one less than ndofs
        nedges = length(Set(minmax(e...) for c in getcells(grid) for e in Ferrite.edges(c)))
        @test ndofs(dh) == 2 * nedges + 1
    end

    @testset "global space and domain topology" begin
        for hole in (false, true)
            grid = generate_grid(Triangle, (3, 3), Vec((0.0, 0.0)), Vec((3.0, 3.0)))
            # Remove the central square to obtain a connected mesh with one hole.
            cells = filter(getcells(grid)) do cell
                x = sum(grid.nodes[i].x for i in cell.nodes) / 3
                return !hole || !(1 < x[1] < 2 && 1 < x[2] < 2)
            end
            grid = Grid(cells, grid.nodes)
            dh = DofHandler(grid)
            add!(dh, :u, ip)
            close!(dh)

            # Map global coefficients to broken P2, represented by six nodal
            # values on each cell. Its rank is the implemented space dimension.
            B = zeros(6 * length(cells), ndofs(dh))
            for k in eachindex(cells), j in 1:6, i in 1:7
                ξ = Ferrite.reference_coordinates(lag)[j]
                B[6 * (k - 1) + j, celldofs(dh, k)[i]] = reference_shape_value(ip, ξ, i)
            end

            # Independently impose Gauss continuity on broken P2. The kernel
            # of C is the full Gauss-continuous space, including any hole modes.
            edges = Dict{Tuple{Int, Int}, Vector{Int}}()
            for (k, cell) in enumerate(cells), edge in Ferrite.edges(cell)
                push!(get!(edges, minmax(edge...), Int[]), k)
            end
            ninterior = count(ks -> length(ks) == 2, values(edges))
            C = zeros(2 * ninterior, 6 * length(cells))
            row = 0
            for ((a, b), ks) in edges
                length(ks) == 2 || continue
                for x in gauss_points(grid.nodes[a].x, grid.nodes[b].x)
                    row += 1
                    for (sign, k) in zip((1, -1), ks)
                        X = [grid.nodes[i].x for i in cells[k].nodes]
                        ξ = Vec{2}(Tuple(hcat(X[1] - X[3], X[2] - X[3]) \ (x - X[3])))
                        for j in 1:6
                            C[row, 6 * (k - 1) + j] = sign * reference_shape_value(lag, ξ, j)
                        end
                    end
                end
            end
            @test norm(C * B) < 1.0e-12
            @test rank(B) == ndofs(dh) - 1
            @test size(C, 2) - rank(C) == rank(B) + Int(hole)
        end
    end

    # (v) Dirichlet conditions act on the vertex and midpoint DOFs only
    @testset "Dirichlet on the Lagrange DOFs" begin
        grid = generate_grid(Triangle, (3, 3))
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        ch = ConstraintHandler(dh)
        ∂Ω = union(getfacetset.((grid,), ["left", "right", "bottom", "top"])...)
        add!(ch, Dirichlet(:u, ∂Ω, x -> x[1] + 2x[2]))
        close!(ch)
        update!(ch, 0.0)
        # 12 boundary vertices + 12 boundary edge midpoints
        @test length(ch.prescribed_dofs) == 24
        celldofs7 = Set(celldofs(dh, i)[7] for i in 1:getncells(grid))
        @test isempty(intersect(Set(ch.prescribed_dofs), celldofs7))
    end

    # (vi) Poisson convergence: second order in the broken H1 seminorm and
    # third order in L2 (paper: the element is second-order accurate)
    @testset "Poisson convergence" begin
        uex(x) = sin(π * x[1]) * sin(π * x[2])
        ∇uex(x) = π * Vec((cos(π * x[1]) * sin(π * x[2]), sin(π * x[1]) * cos(π * x[2])))
        fsrc(x) = 2π^2 * uex(x)
        function solve(n)
            grid = generate_grid(Triangle, (n, n), Vec((0.0, 0.0)), Vec((1.0, 1.0)))
            dh = DofHandler(grid)
            add!(dh, :u, ip)
            close!(dh)
            ch = ConstraintHandler(dh)
            ∂Ω = union(getfacetset.((grid,), ["left", "right", "bottom", "top"])...)
            add!(ch, Dirichlet(:u, ∂Ω, x -> 0.0))
            close!(ch)
            cv = CellValues(QuadratureRule{RefTriangle}(6), ip, Lagrange{RefTriangle, 1}())
            K = allocate_matrix(dh)
            f = zeros(ndofs(dh))
            asm = start_assemble(K, f)
            Ke, fe = zeros(7, 7), zeros(7)
            for cell in CellIterator(dh)
                reinit!(cv, cell)
                fill!(Ke, 0)
                fill!(fe, 0)
                for qp in 1:getnquadpoints(cv)
                    dΩ = getdetJdV(cv, qp)
                    x = spatial_coordinate(cv, qp, getcoordinates(cell))
                    for i in 1:7
                        fe[i] += fsrc(x) * shape_value(cv, qp, i) * dΩ
                        for j in 1:7
                            Ke[i, j] += dot(shape_gradient(cv, qp, i), shape_gradient(cv, qp, j)) * dΩ
                        end
                    end
                end
                assemble!(asm, celldofs(cell), Ke, fe)
            end
            apply!(K, f, ch)
            u = K \ f
            apply!(u, ch)
            eL2, eH1 = 0.0, 0.0
            for cell in CellIterator(dh)
                reinit!(cv, cell)
                ue = u[celldofs(cell)]
                for qp in 1:getnquadpoints(cv)
                    dΩ = getdetJdV(cv, qp)
                    x = spatial_coordinate(cv, qp, getcoordinates(cell))
                    eL2 += (function_value(cv, qp, ue) - uex(x))^2 * dΩ
                    eH1 += norm(function_gradient(cv, qp, ue) - ∇uex(x))^2 * dΩ
                end
            end
            return sqrt(eL2), sqrt(eH1)
        end
        (l2a, h1a), (l2b, h1b) = solve(8), solve(16)
        @test log2(h1a / h1b) > 1.9
        @test log2(l2a / l2b) > 2.8
    end

    # (vii) Boundary conditions. Dirichlet data is imposed on the Lagrange DOFs,
    # but the neutral function does not vanish on the boundary, so the trace
    # matches the data only at the two Gauss-Legendre points of every boundary
    # edge (the nonconforming sense of the boundary condition), not at the
    # vertices and midpoints.
    @testset "boundary conditions" begin
        gauss = ((3 - sqrt(3)) / 6, (3 + sqrt(3)) / 6)
        test_dirichlet_bc(ip, bc_poly(2); trace_points = gauss)
        test_neumann_bc(ip, bc_poly(2))
    end
end
