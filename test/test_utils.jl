# Shared test helpers. Element-specific data (symfem reference tables etc.)
# lives in the per-element test files.

using Ferrite
using Ferrite: getrefshape, getnbasefunctions, reference_shape_value
using LinearAlgebra: norm
using PythonCall
using Test

# Ferrite's generic interpolation property checker (test_interpolation_properties).
include(joinpath(pkgdir(Ferrite), "test", "interpolation_test_utils.jl"))

const symfem = pyimport("symfem")

# Sample a (not necessarily uniformly distributed) random point strictly inside
# the reference cell.
function sample_reference_point(::Type{Ferrite.RefHypercube{dim}}) where {dim}
    return Vec{dim}(ntuple(_ -> 2 * rand() - 1, dim))
end
function sample_reference_point(::Type{Ferrite.RefSimplex{dim}}) where {dim}
    # Dirichlet(1, ..., 1) via normalized exponentials; the reference simplex is
    # {ξ ≥ 0, sum(ξ) ≤ 1} for both RefTriangle and RefTetrahedron.
    w = ntuple(_ -> -log(rand()), dim + 1)
    s = sum(w)
    return Vec{dim}(ntuple(i -> w[i] / s, dim))
end
function sample_reference_point(::Type{RefPrism})
    tri = sample_reference_point(RefTriangle)
    return Vec{3}((tri[1], tri[2], rand()))
end
function sample_reference_point(::Type{RefPyramid})
    z = rand()
    return Vec{3}((rand() * (1 - z), rand() * (1 - z), z))
end

# Σᵢ Nᵢ(ξ) == 1 at random reference points (any partition-of-unity element).
function test_partition_of_unity(ip; npoints = 20)
    @testset "partition of unity: $ip" begin
        for _ in 1:npoints
            ξ = sample_reference_point(getrefshape(ip))
            s = sum(i -> reference_shape_value(ip, ξ, i), 1:getnbasefunctions(ip))
            @test s ≈ 1 atol = 1.0e-12
        end
    end
end

# Nᵢ(ξⱼ) == δᵢⱼ at reference_coordinates (nodal elements only).
function test_kronecker_delta(ip)
    @testset "Kronecker delta: $ip" begin
        coords = Ferrite.reference_coordinates(ip)
        N = getnbasefunctions(ip)
        @test length(coords) == N
        for (j, ξ) in pairs(coords), i in 1:N
            @test reference_shape_value(ip, ξ, i) ≈ (i == j ? 1.0 : 0.0) atol = 1.0e-12
        end
    end
end

# Exponent tuples for the monomial basis of P_degree (:P, total degree) or
# Q_degree (:Q, per-variable degree) in `dim` variables.
function monomial_exponents(space::Symbol, dim::Int, degree::Int)
    ranges = ntuple(_ -> 0:degree, dim)
    exps = vec(collect(Iterators.product(ranges...)))
    space === :Q && return exps
    space === :P && return filter(e -> sum(e) <= degree, exps)
    return error("unknown polynomial space $space")
end

# Check that every monomial of P_degree/Q_degree lies in the span of the shape
# functions: least-squares fit at random points, residual ≈ 0. Catches basis
# transcription typos without needing the dual basis.
function test_polynomial_reproduction(ip, degree::Int; space::Symbol)
    @testset "$space$degree ⊆ span: $ip" begin
        shape = getrefshape(ip)
        N = getnbasefunctions(ip)
        npts = 3 * N + 10
        points = [sample_reference_point(shape) for _ in 1:npts]
        A = [reference_shape_value(ip, ξ, i) for ξ in points, i in 1:N]
        for e in monomial_exponents(space, Ferrite.getrefdim(ip), degree)
            b = [prod(ξ .^ e) for ξ in points]
            c = A \ b
            @test norm(A * c - b) < 1.0e-10 * max(1, norm(b))
        end
    end
end

# --- symfem cross-check (via PythonCall; test-only dependency) ---------------

# DefElement/symfem cell name for a Ferrite reference shape.
symfem_cellname(::Type{RefLine}) = "interval"
symfem_cellname(::Type{RefTriangle}) = "triangle"
symfem_cellname(::Type{RefQuadrilateral}) = "quadrilateral"
symfem_cellname(::Type{RefTetrahedron}) = "tetrahedron"
symfem_cellname(::Type{RefHexahedron}) = "hexahedron"
symfem_cellname(::Type{RefPrism}) = "prism"
symfem_cellname(::Type{RefPyramid}) = "pyramid"

# Map a point from the Ferrite reference cell to the symfem/DefElement one:
# hypercubes are [-1, 1]^d in Ferrite but [0, 1]^d in symfem; the remaining
# cells have identical coordinates (up to entity numbering, which the
# permutation handles).
symfem_coords(::Type{<:Ferrite.RefHypercube}, ξ::Vec) = (ξ .+ 1) ./ 2
symfem_coords(::Type{<:Ferrite.AbstractRefShape}, ξ::Vec) = ξ

"""
    test_symfem_reference(ip, family, degree, perm; npoints = 10)

Compare `reference_shape_value` of `ip` against the symfem element
`create_element(cell, family, degree)` at `npoints` random reference points.
`perm` maps Ferrite DOF `i` to the 0-based symfem DOF `perm[i]`; it accounts
for the differing entity numbering (and must be derived per element, e.g.
geometrically from the DOF points -- symfem's `entity_dofs` entity numbering
is not consistent across element families).
"""
function test_symfem_reference(ip, family::String, degree::Int, perm::Vector{Int}; npoints = 10, scales = nothing, kwargs...)
    @testset "symfem cross-check: $ip" begin
        shape = getrefshape(ip)
        N = getnbasefunctions(ip)
        el = symfem.create_element(symfem_cellname(shape), family, degree; kwargs...)
        @test pyconvert(Int, el.space_dim) == N
        @test sort(perm) == 0:(N - 1)
        basis = el.get_basis_functions()
        x = symfem.symbols.x
        for _ in 1:npoints
            ξ = sample_reference_point(shape)
            sp = pytuple(Tuple(symfem_coords(shape, ξ)))
            for i in 1:N
                expected = pyconvert(Float64, pybuiltins.float(basis[perm[i]].subs(x, sp).as_sympy()))
                scales !== nothing && (expected *= scales[i])
                @test reference_shape_value(ip, ξ, i) ≈ expected atol = 1.0e-12
            end
        end
    end
end

"""
    test_symfem_reference_vector(ip, family, degree, sperm; npoints = 10)

Vector-valued version of [`test_symfem_reference`](@ref): `sperm[i]` is a
`(sign, j)` tuple mapping Ferrite DOF `i` to `sign` times the 0-based symfem
DOF `j` (H(div)/H(curl) conventions differ by edge order, intra-edge weight
order and normal/tangent sign, all absorbed in the signed permutation).
"""
function test_symfem_reference_vector(ip, family::String, degree::Int, sperm::Vector{<:Tuple{Real, Int}}; npoints = 10, kwargs...)
    @testset "symfem cross-check: $ip" begin
        shape = getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        N = getnbasefunctions(ip)
        el = symfem.create_element(symfem_cellname(shape), family, degree; kwargs...)
        @test pyconvert(Int, el.space_dim) == N
        @test sort(last.(sperm)) == 0:(N - 1)
        basis = el.get_basis_functions()
        x = symfem.symbols.x
        for _ in 1:npoints
            ξ = sample_reference_point(shape)
            sp = pytuple(Tuple(symfem_coords(shape, ξ)))
            for i in 1:N
                sign, j = sperm[i]
                fj = basis[j].subs(x, sp)
                expected = Vec{dim}(c -> sign * pyconvert(Float64, pybuiltins.float(fj[c - 1].as_sympy())))
                @test reference_shape_value(ip, ξ, i) ≈ expected atol = 1.0e-12
            end
        end
    end
end

# --- H(div) helpers ----------------------------------------------------------

# Two-cell DofHandler + normal-continuity test for H(div) elements: cells 1
# and 2 of `grid` share the edge (facet1 of cell 1, facet2 of cell 2, with
# opposite orientation in the test grids); `nshared` edge DOFs are identified.
# The normal component of an arbitrary field must agree at matching physical
# points on the shared edge; if `tangential_moment`, the zeroth tangential
# moment must also agree (MTW/HZ-type elements).
function test_hdiv_two_cell(ip, grid, facet1::Int, facet2::Int, nshared::Int; tangential_moment = false)
    @testset "two-cell H(div): $ip" begin
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        N = getnbasefunctions(ip)
        @test ndofs(dh) == 2N - nshared
        @test length(intersect(celldofs(dh, 1), celldofs(dh, 2))) == nshared
        u = rand(ndofs(dh))
        shape = getrefshape(ip)
        fqr = FacetQuadratureRule{shape}(4)
        geo = Lagrange{shape, 1}()
        fv1 = FacetValues(fqr, ip, geo)
        fv2 = FacetValues(fqr, ip, geo)
        coords1 = getcoordinates(grid, 1)
        coords2 = getcoordinates(grid, 2)
        reinit!(fv1, getcells(grid, 1), coords1, facet1)
        reinit!(fv2, getcells(grid, 2), coords2, facet2)
        u1 = u[celldofs(dh, 1)]
        u2 = u[celldofs(dh, 2)]
        tmoment = 0.0
        for qp1 in 1:getnquadpoints(fv1)
            x1 = spatial_coordinate(fv1, qp1, coords1)
            qp2 = findfirst(qp -> norm(spatial_coordinate(fv2, qp, coords2) - x1) < 1.0e-12, 1:getnquadpoints(fv2))
            @test qp2 !== nothing
            n1 = getnormal(fv1, qp1)
            t = Vec(-n1[2], n1[1])
            v1 = function_value(fv1, qp1, u1)
            v2 = function_value(fv2, qp2, u2)
            @test v1 ⋅ n1 ≈ v2 ⋅ n1 rtol = 1.0e-10 atol = 1.0e-12
            tmoment += ((v1 - v2) ⋅ t) * getdetJdV(fv1, qp1)
        end
        if tangential_moment
            @test tmoment ≈ 0 atol = 1.0e-11
        end
    end
end

# Two-cell DofHandler + tangential-continuity test for H(curl) elements,
# mirroring `test_hdiv_two_cell`.
function test_hcurl_two_cell(ip, grid, facet1::Int, facet2::Int, nshared::Int)
    @testset "two-cell H(curl): $ip" begin
        dh = DofHandler(grid)
        add!(dh, :u, ip)
        close!(dh)
        N = getnbasefunctions(ip)
        @test ndofs(dh) == 2N - nshared
        @test length(intersect(celldofs(dh, 1), celldofs(dh, 2))) == nshared
        u = rand(ndofs(dh))
        shape = getrefshape(ip)
        fqr = FacetQuadratureRule{shape}(4)
        geo = Lagrange{shape, 1}()
        fv1 = FacetValues(fqr, ip, geo)
        fv2 = FacetValues(fqr, ip, geo)
        coords1 = getcoordinates(grid, 1)
        coords2 = getcoordinates(grid, 2)
        reinit!(fv1, getcells(grid, 1), coords1, facet1)
        reinit!(fv2, getcells(grid, 2), coords2, facet2)
        u1 = u[celldofs(dh, 1)]
        u2 = u[celldofs(dh, 2)]
        for qp1 in 1:getnquadpoints(fv1)
            x1 = spatial_coordinate(fv1, qp1, coords1)
            qp2 = findfirst(qp -> norm(spatial_coordinate(fv2, qp, coords2) - x1) < 1.0e-12, 1:getnquadpoints(fv2))
            @test qp2 !== nothing
            n1 = getnormal(fv1, qp1)
            t = Vec(-n1[2], n1[1])
            v1 = function_value(fv1, qp1, u1)
            v2 = function_value(fv2, qp2, u2)
            @test v1 ⋅ t ≈ v2 ⋅ t rtol = 1.0e-10 atol = 1.0e-12
        end
    end
end

# Curl/Stokes theorem per basis function on a single cell: int curl(N_i) dA
# equals the counterclockwise circulation (2D).
function test_curl_theorem(ip, cell, coords)
    @testset "curl theorem: $ip" begin
        shape = getrefshape(ip)
        geo = Lagrange{shape, 1}()
        cv = CellValues(QuadratureRule{shape}(4), ip, geo)
        reinit!(cv, cell, coords)
        fv = FacetValues(FacetQuadratureRule{shape}(4), ip, geo)
        for i in 1:getnbasefunctions(ip)
            curlint = sum(
                (g = shape_gradient(cv, qp, i); (g[2, 1] - g[1, 2]) * getdetJdV(cv, qp))
                    for qp in 1:getnquadpoints(cv)
            )
            circ = 0.0
            for facet in 1:Ferrite.nfacets(ip)
                reinit!(fv, cell, coords, facet)
                for qp in 1:getnquadpoints(fv)
                    n = getnormal(fv, qp)
                    t = Vec(-n[2], n[1])
                    circ += (shape_value(fv, qp, i) ⋅ t) * getdetJdV(fv, qp)
                end
            end
            @test curlint ≈ circ atol = 1.0e-11
        end
    end
end

# Divergence theorem per basis function on a single cell: int div(N_i) dV
# equals the total boundary flux (checks the contravariant Piola mapping of
# values and gradients consistently).
function test_divergence_theorem(ip, cell, coords)
    @testset "divergence theorem: $ip" begin
        shape = getrefshape(ip)
        geo = Lagrange{shape, 1}()
        cv = CellValues(QuadratureRule{shape}(4), ip, geo)
        reinit!(cv, cell, coords)
        fv = FacetValues(FacetQuadratureRule{shape}(4), ip, geo)
        for i in 1:getnbasefunctions(ip)
            vol = sum(shape_divergence(cv, qp, i) * getdetJdV(cv, qp) for qp in 1:getnquadpoints(cv))
            flux = 0.0
            for facet in 1:Ferrite.nfacets(ip)
                reinit!(fv, cell, coords, facet)
                flux += sum(
                    (shape_value(fv, qp, i) ⋅ getnormal(fv, qp)) * getdetJdV(fv, qp)
                        for qp in 1:getnquadpoints(fv)
                )
            end
            @test vol ≈ flux atol = 1.0e-11
        end
    end
end

# Evaluate value and AD gradient with different float types (smoke test for
# type-generic implementations).
function test_type_genericity(ip)
    @testset "type genericity: $ip" begin
        shape = getrefshape(ip)
        dim = Ferrite.getrefdim(ip)
        for T in (Float32, Float64)
            ξ = Vec{dim, T}(sample_reference_point(shape))
            for i in 1:getnbasefunctions(ip)
                v = reference_shape_value(ip, ξ, i)
                @test v isa Ferrite.shape_value_type(ip, T)
                g = Ferrite.reference_shape_gradient(ip, ξ, i)
                @test all(isfinite, g)
            end
        end
    end
end
