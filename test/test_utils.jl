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
const sympy = pyimport("sympy")

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
    return @testset "partition of unity: $ip" begin
        for _ in 1:npoints
            ξ = sample_reference_point(getrefshape(ip))
            s = sum(i -> reference_shape_value(ip, ξ, i), 1:getnbasefunctions(ip))
            @test s ≈ 1 atol = 1.0e-12
        end
    end
end

# Nᵢ(ξⱼ) == δᵢⱼ at reference_coordinates (nodal elements only).
function test_kronecker_delta(ip)
    return @testset "Kronecker delta: $ip" begin
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
    return @testset "$space$degree ⊆ span: $ip" begin
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

# `ξ` mapped to the symfem cell as a tuple of exact sympy rationals, for
# substitution into symfem basis functions. These are expanded polynomials with
# large rational coefficients, so substituting floats loses digits to
# cancellation (errors > 1e-12 for e.g. Lagrange{RefHexahedron, 3}); the exact
# evaluation is rounded to Float64 only at the end.
function symfem_point(shape, ξ::Vec)
    s = symfem_coords(shape, Vec(Rational{BigInt}.(ξ.data)))
    return pytuple(Tuple(sympy.Rational(string(numerator(c)), string(denominator(c))) for c in s))
end

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
    return @testset "symfem cross-check: $ip" begin
        shape = getrefshape(ip)
        N = getnbasefunctions(ip)
        el = symfem.create_element(symfem_cellname(shape), family, degree; kwargs...)
        @test pyconvert(Int, el.space_dim) == N
        @test sort(perm) == 0:(N - 1)
        basis = el.get_basis_functions()
        x = symfem.symbols.x
        for _ in 1:npoints
            ξ = sample_reference_point(shape)
            sp = symfem_point(shape, ξ)
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
    return @testset "symfem cross-check: $ip" begin
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
            sp = symfem_point(shape, ξ)
            for i in 1:N
                sign, j = sperm[i]
                fj = basis[j].subs(x, sp)
                expected = Vec{dim}(c -> sign * pyconvert(Float64, pybuiltins.float(fj[c - 1].as_sympy())))
                @test reference_shape_value(ip, ξ, i) ≈ expected atol = 1.0e-12
            end
        end
    end
end

"""
    symfem_point_perm(ip, family, degree; kwargs...)

Permutation for [`test_symfem_reference`](@ref) of a nodal element, derived
geometrically: Ferrite DOF `i` (at `reference_coordinates(ip)[i]`) maps to the
0-based symfem point-evaluation DOF at the same (mapped) reference point.
"""
function symfem_point_perm(ip, family::String, degree::Int; kwargs...)
    shape = getrefshape(ip)
    dim = Ferrite.getrefdim(ip)
    el = symfem.create_element(symfem_cellname(shape), family, degree; kwargs...)
    spts = [Vec{dim}(ntuple(j -> pyconvert(Float64, pybuiltins.float(d.point[j - 1])), dim)) for d in el.dofs]
    fpts = [symfem_coords(shape, ξ) for ξ in Ferrite.reference_coordinates(ip)]
    return [findfirst(p -> norm(p - fp) < 1.0e-10, spts) - 1 for fp in fpts]
end

"""
    test_symfem_span(ip, family, degree; kwargs...)

Check that the (scalar) basis of `ip` spans the same space as the symfem
element `create_element(cell, family, degree)`, for elements whose basis
differs from symfem's because the DOFs differ (least-squares fits both ways
at random points, residuals ≈ 0).
"""
function test_symfem_span(ip, family::String, degree::Int; kwargs...)
    return @testset "symfem span: $ip" begin
        shape = getrefshape(ip)
        N = getnbasefunctions(ip)
        el = symfem.create_element(symfem_cellname(shape), family, degree; kwargs...)
        @test pyconvert(Int, el.space_dim) == N
        basis = el.get_basis_functions()
        x = symfem.symbols.x
        points = [sample_reference_point(shape) for _ in 1:(3 * N + 10)]
        A = [reference_shape_value(ip, ξ, i) for ξ in points, i in 1:N]
        B = [
            pyconvert(Float64, pybuiltins.float(basis[j].subs(x, symfem_point(shape, ξ)).as_sympy()))
                for ξ in points, j in 0:(N - 1)
        ]
        @test norm(A - B * (B \ A)) < 1.0e-10 * norm(A)
        @test norm(B - A * (A \ B)) < 1.0e-10 * norm(B)
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
    return @testset "two-cell H(div): $ip" begin
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
    return @testset "two-cell H(curl): $ip" begin
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
    return @testset "curl theorem: $ip" begin
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
    return @testset "divergence theorem: $ip" begin
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
    return @testset "type genericity: $ip" begin
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

# --- Boundary conditions -----------------------------------------------------
#
# The elements are defined on a single reference cell; these helpers check that
# they also work in a mesh with boundary conditions: strong Dirichlet
# conditions through `Dirichlet` (nodal) or `ProjectedDirichlet` (H(div)/
# H(curl)), weak Dirichlet conditions (SIPG/Nitsche) for discontinuous elements,
# and natural (Neumann) conditions. The manufactured solution lies in the
# discrete space, so the FE solution must be exact and the boundary trace must
# match the prescribed data.

using LinearAlgebra: qr

const BC_FACETSETS = ("left", "right", "bottom", "top", "front", "back")

# Unstructured-ish test mesh on [-1, 1]^dim with facetset "∂Ω". Cells stay affine
# (so that the polynomial solutions below lie in the mapped spaces): interior
# nodes of simplex/line meshes are perturbed, tensor-product meshes are
# sheared into parallelograms.
function bc_test_grid(::Type{CT}; n = 3) where {CT}
    dim = Ferrite.getrefdim(CT)
    grid = generate_grid(CT, ntuple(_ -> n, dim))
    ∂Ω = union((getfacetset(grid, k) for k in BC_FACETSETS if haskey(grid.facetsets, k))...)
    addfacetset!(grid, "∂Ω", ∂Ω)
    if Ferrite.getrefshape(CT) <: Ferrite.RefHypercube && dim > 1
        transform_coordinates!(grid, x -> Vec{dim}(i -> x[i] + (i == 1 ? 0.3 * x[2] : 0.1 * x[1])))
    else
        bnodes = Set(n for (c, f) in ∂Ω for n in Ferrite.facets(getcells(grid, c))[f])
        for i in eachindex(grid.nodes)
            i in bnodes && continue
            grid.nodes[i] = Node(grid.nodes[i].x + Vec{dim}(d -> 0.15 / n * sin(3i + d)))
        end
    end
    return grid
end

function bc_test_dofhandler(ip, grid)
    return close!(add!(DofHandler(grid), :u, ip))
end

# Iterate over (sdh, ip) of field :u.
bc_subdofhandlers(dh) = ((sdh, Ferrite.getfieldinterpolation(sdh, :u)) for sdh in dh.subdofhandlers)

# Scalar model problem -Δu + u = f (broken gradient for nonconforming
# elements) with either Dirichlet data `uex` on ∂Ω (`dirichlet = true`) or
# Neumann data ∇uex ⋅ n on ∂Ω; f = -Δuex + uex unless given. Solved with a
# (rank-revealing) sparse QR since some elements (Fortin-Soulie) have a
# non-unique representation.
function solve_scalar_bc(dh, uex; dirichlet::Bool, qr_order::Int, f = x -> -tr(hessian(uex, x)) + uex(x))
    grid = Ferrite.get_grid(dh)
    ∂Ω = getfacetset(grid, "∂Ω")
    ch = ConstraintHandler(dh)
    dirichlet && add!(ch, Dirichlet(:u, ∂Ω, uex))
    close!(ch)
    K = allocate_matrix(dh)
    F = zeros(ndofs(dh))
    asm = start_assemble(K, F)
    for (sdh, ip) in bc_subdofhandlers(dh)
        shape = getrefshape(ip)
        geo = geometric_interpolation(getcelltype(sdh))
        cv = CellValues(QuadratureRule{shape}(qr_order), ip, geo)
        fv = FacetValues(FacetQuadratureRule{shape}(qr_order), ip, geo)
        ΓN = dirichlet ? FacetIndex[] : filter(fi -> fi[1] in sdh.cellset, ∂Ω)
        assemble_scalar_bc!(asm, F, sdh, cv, fv, ΓN, f, uex)
    end
    update!(ch, 0.0)
    apply!(K, F, ch)
    u = qr(K) \ F
    apply!(u, ch)
    return u
end

# Function barrier for solve_scalar_bc (concretely typed `cv`, `fv`); Neumann
# data ∇uex ⋅ n on the facets ΓN.
function assemble_scalar_bc!(asm, F, sdh, cv, fv, ΓN, f, uex)
    n = getnbasefunctions(cv)
    Ke = zeros(n, n)
    Fe = zeros(n)
    for cc in CellIterator(sdh)
        reinit!(cv, cc)
        fill!(Ke, 0)
        fill!(Fe, 0)
        for qp in 1:getnquadpoints(cv)
            dV = getdetJdV(cv, qp)
            fx = f(spatial_coordinate(cv, qp, getcoordinates(cc)))
            for i in 1:n
                Ni, ∇Ni = shape_value(cv, qp, i), shape_gradient(cv, qp, i)
                Fe[i] += Ni * fx * dV
                for j in 1:n
                    Ke[i, j] += (∇Ni ⋅ shape_gradient(cv, qp, j) + Ni * shape_value(cv, qp, j)) * dV
                end
            end
        end
        assemble!(asm, celldofs(cc), Ke, Fe)
    end
    isempty(ΓN) && return
    for fc in FacetIterator(sdh, ΓN)
        reinit!(fv, fc)
        dofs = celldofs(fc)
        for qp in 1:getnquadpoints(fv)
            x = spatial_coordinate(fv, qp, getcoordinates(fc))
            h = gradient(uex, x) ⋅ getnormal(fv, qp)
            for i in 1:n
                F[dofs[i]] += h * shape_value(fv, qp, i) * getdetJdV(fv, qp)
            end
        end
    end
    return
end

# Points on every boundary facet, given by the facet parameters `s` in [0, 1]
# along the reference edge (2D); the facet vertex itself in 1D.
function bc_facet_point_rules(shape, s)
    verts = Ferrite.reference_coordinates(Lagrange{shape, 1}())
    return map(Ferrite.reference_facets(shape)) do fvs
        pts = length(fvs) == 1 ? [verts[fvs[1]]] : [verts[fvs[1]] + si * (verts[fvs[2]] - verts[fvs[1]]) for si in s]
        QuadratureRule{shape}(ones(length(pts)), pts)
    end
end

# max |u_h(x) - g(x)| over the points `s` of every boundary facet. `s` may be a
# function of the facet's cell id (for meshes mixing elements).
function max_boundary_error(dh, u, g; s = range(0, 1, length = 9))
    grid = Ferrite.get_grid(dh)
    err = 0.0
    for (sdh, ip) in bc_subdofhandlers(dh)
        geo = geometric_interpolation(getcelltype(sdh))
        for (c, f) in getfacetset(grid, "∂Ω")
            c in sdh.cellset || continue
            sc = s isa Function ? s(c, f) : s
            isempty(sc) && continue
            cv = CellValues(bc_facet_point_rules(getrefshape(ip), sc)[f], ip, geo)
            coords = getcoordinates(grid, c)
            reinit!(cv, getcells(grid, c), coords)
            ue = u[celldofs(dh, c)]
            for qp in 1:getnquadpoints(cv)
                err = max(err, abs(function_value(cv, qp, ue) - g(spatial_coordinate(cv, qp, coords))))
            end
        end
    end
    return err
end

# max |u_h - uex| (scalar) or ‖u_h - uex‖ (vector) at interior quadrature points.
function max_domain_error(dh, u, uex)
    err = 0.0
    for (sdh, ip) in bc_subdofhandlers(dh)
        cv = CellValues(QuadratureRule{getrefshape(ip)}(4), ip, geometric_interpolation(getcelltype(sdh)))
        for cc in CellIterator(sdh)
            reinit!(cv, cc)
            ue = u[celldofs(cc)]
            for qp in 1:getnquadpoints(cv)
                err = max(err, norm(function_value(cv, qp, ue) - uex(spatial_coordinate(cv, qp, getcoordinates(cc)))))
            end
        end
    end
    return err
end

"""
    test_dirichlet_bc(ip_or_dh, uex; trace_points, exact = true, qr_order, n = 3)

Solve -Δu + u = f with the strong Dirichlet condition u = uex on ∂Ω (inhomogeneous)
and with u = 0 (homogeneous, f = 1). Check that the boundary trace of the solution
matches the data at the facet points `trace_points` (default: the whole facet,
see [`max_boundary_error`](@ref)) and, if `exact`, that u_h == uex in the domain
(`uex` must then lie in the discrete space).
"""
function test_dirichlet_bc(dh::DofHandler, uex; trace_points = range(0, 1, length = 9), exact = true, qr_order::Int)
    return @testset "Dirichlet BC" begin
        u = solve_scalar_bc(dh, uex; dirichlet = true, qr_order)
        @test max_boundary_error(dh, u, uex; s = trace_points) < 1.0e-12
        exact && @test max_domain_error(dh, u, uex) < 1.0e-11
        u0 = solve_scalar_bc(dh, x -> 0.0; dirichlet = true, qr_order, f = x -> 1.0)
        @test norm(u0) > 0.01
        @test max_boundary_error(dh, u0, x -> 0.0; s = trace_points) < 1.0e-12
    end
end
function test_dirichlet_bc(ip::Interpolation, uex; n = 3, qr_order = 2 * Ferrite.getorder(ip) + 2, kwargs...)
    grid = bc_test_grid(bc_celltype(getrefshape(ip)); n)
    return test_dirichlet_bc(bc_test_dofhandler(ip, grid), uex; qr_order, kwargs...)
end

"""
    test_neumann_bc(ip_or_dh, uex; qr_order, n = 3)

Solve -Δu + u = f with the natural condition ∇u ⋅ n = ∇uex ⋅ n on all of ∂Ω and
check u_h == uex (`uex` must lie in the discrete space and the element must
pass the patch test for it).
"""
function test_neumann_bc(dh::DofHandler, uex; qr_order::Int)
    return @testset "Neumann BC" begin
        u = solve_scalar_bc(dh, uex; dirichlet = false, qr_order)
        @test max_domain_error(dh, u, uex) < 1.0e-11
    end
end
function test_neumann_bc(ip::Interpolation, uex; n = 3, qr_order = 2 * Ferrite.getorder(ip) + 2)
    grid = bc_test_grid(bc_celltype(getrefshape(ip)); n)
    return test_neumann_bc(bc_test_dofhandler(ip, grid), uex; qr_order)
end

bc_celltype(::Type{RefLine}) = Line
bc_celltype(::Type{RefTriangle}) = Triangle
bc_celltype(::Type{RefQuadrilateral}) = Quadrilateral
bc_celltype(::Type{RefTetrahedron}) = Tetrahedron
bc_celltype(::Type{RefHexahedron}) = Hexahedron

# 2D scalar curl of a vector gradient.
bc_curl(G) = G[2, 1] - G[1, 2]

# Vector model problem (q, δq) + (D q, D δq) = (f, δq) [+ natural BC], with
# D = div (H(div)) or curl (H(curl)) and f = qex - ∇ div qex resp.
# f = qex + rot curl qex (rot s = (∂s/∂y, -∂s/∂x)). With `essential`, the normal
# (H(div)) or tangential (H(curl)) trace of `bc` is prescribed on ∂Ω via
# `ProjectedDirichlet`; otherwise the natural condition D q = D qex enters
# through ∫ D qex (δq ⋅ n) resp. ∫ curl qex (δq ⋅ t).
function solve_vector_bc(dh, qex; essential::Bool, qr_order::Int, bc = qex)
    grid = Ferrite.get_grid(dh)
    ∂Ω = getfacetset(grid, "∂Ω")
    ip = only(ip for (_, ip) in bc_subdofhandlers(dh))
    hdiv = Ferrite.conformity(ip) isa Ferrite.HdivConformity
    D(G) = hdiv ? tr(G) : bc_curl(G)
    trace(q, n) = hdiv ? q ⋅ n : q ⋅ Vec((-n[2], n[1]))
    ch = ConstraintHandler(dh)
    essential && add!(ch, ProjectedDirichlet(:u, ∂Ω, hdiv ? ((x, t, n) -> bc(x) ⋅ n) : ((x, t, n) -> bc(x) × n)))
    close!(ch)
    update!(ch, 0.0)
    shape = getrefshape(ip)
    geo = geometric_interpolation(getcelltype(grid))
    cv = CellValues(QuadratureRule{shape}(qr_order), ip, geo)
    fv = FacetValues(FacetQuadratureRule{shape}(qr_order), ip, geo)
    K = allocate_matrix(dh)
    F = zeros(ndofs(dh))
    Dqex(x) = D(gradient(qex, x))
    f(x) = (g = gradient(Dqex, x); hdiv ? qex(x) - g : qex(x) + Vec((g[2], -g[1])))
    assemble_vector_bc!(start_assemble(K, F), F, dh, cv, fv, essential ? FacetIndex[] : ∂Ω, f, Dqex, D, trace)
    apply!(K, F, ch)
    u = K \ F
    apply!(u, ch)
    return u
end

# Function barrier for solve_vector_bc; natural boundary term on the facets ΓN.
function assemble_vector_bc!(asm, F, dh, cv, fv, ΓN, f, Dqex, D, trace)
    n = getnbasefunctions(cv)
    Ke = zeros(n, n)
    Fe = zeros(n)
    for cc in CellIterator(dh)
        reinit!(cv, cc)
        fill!(Ke, 0)
        fill!(Fe, 0)
        for qp in 1:getnquadpoints(cv)
            dV = getdetJdV(cv, qp)
            fx = f(spatial_coordinate(cv, qp, getcoordinates(cc)))
            for i in 1:n
                Ni, DNi = shape_value(cv, qp, i), D(shape_gradient(cv, qp, i))
                Fe[i] += Ni ⋅ fx * dV
                for j in 1:n
                    Ke[i, j] += (Ni ⋅ shape_value(cv, qp, j) + DNi * D(shape_gradient(cv, qp, j))) * dV
                end
            end
        end
        assemble!(asm, celldofs(cc), Ke, Fe)
    end
    isempty(ΓN) && return
    for fc in FacetIterator(dh, ΓN)
        reinit!(fv, fc)
        dofs = celldofs(fc)
        for qp in 1:getnquadpoints(fv)
            Dq = Dqex(spatial_coordinate(fv, qp, getcoordinates(fc)))
            for i in 1:n
                F[dofs[i]] += Dq * trace(shape_value(fv, qp, i), getnormal(fv, qp)) * getdetJdV(fv, qp)
            end
        end
    end
    return
end

# max |trace(q_h - qex)| at facet quadrature points of ∂Ω (normal trace for
# H(div), tangential trace for H(curl)).
function max_trace_error(dh, u, qex; qr_order = 6)
    grid = Ferrite.get_grid(dh)
    ip = only(ip for (_, ip) in bc_subdofhandlers(dh))
    hdiv = Ferrite.conformity(ip) isa Ferrite.HdivConformity
    fv = FacetValues(FacetQuadratureRule{getrefshape(ip)}(qr_order), ip, geometric_interpolation(getcelltype(grid)))
    err = 0.0
    for fc in FacetIterator(dh, getfacetset(grid, "∂Ω"))
        reinit!(fv, fc)
        ue = u[celldofs(fc)]
        for qp in 1:getnquadpoints(fv)
            n = getnormal(fv, qp)
            d = function_value(fv, qp, ue) - qex(spatial_coordinate(fv, qp, getcoordinates(fc)))
            err = max(err, abs(hdiv ? d ⋅ n : d ⋅ Vec((-n[2], n[1]))))
        end
    end
    return err
end

"""
    test_vector_bcs(ip; qex, n = 3)

H(div)/H(curl) boundary conditions for a 2D vector element: `ProjectedDirichlet`
with the normal resp. tangential trace of `qex` (inhomogeneous) and of zero
(homogeneous), checked pointwise along ∂Ω, and the natural condition on div resp.
curl. `qex` must be linear and lie in the discrete space (on the affine test
mesh), in which case the solution must be exact.
"""
function test_vector_bcs(ip; qex = x -> Vec((1 + 2x[1] - x[2], -2 + x[1] + 3x[2])), n = 3, qr_order = 2 * Ferrite.getorder(ip) + 2)
    return @testset "BCs: $ip" begin
        grid = bc_test_grid(bc_celltype(getrefshape(ip)); n)
        dh = bc_test_dofhandler(ip, grid)
        @testset "ProjectedDirichlet" begin
            u = solve_vector_bc(dh, qex; essential = true, qr_order)
            @test max_trace_error(dh, u, qex) < 1.0e-12
            @test max_domain_error(dh, u, qex) < 1.0e-11
            zq(x) = zero(Vec{2, Float64})
            u0 = solve_vector_bc(dh, qex; essential = true, qr_order, bc = zq)
            @test norm(u0) > 0.01
            @test max_trace_error(dh, u0, zq) < 1.0e-12
        end
        @testset "natural BC" begin
            u = solve_vector_bc(dh, qex; essential = false, qr_order)
            @test max_domain_error(dh, u, qex) < 1.0e-11
        end
    end
end

# Weak Dirichlet boundary: the facets of ∂Ω whose cell centroid has x < 0.
function weak_dirichlet_facets(grid)
    center(c) = sum(getcoordinates(grid, c)) / length(getcoordinates(grid, c))
    return Set(fi for fi in getfacetset(grid, "∂Ω") if center(fi[1])[1] < 0)
end

# Symmetric interior penalty (SIPG) discretization of -Δu + u = f for
# discontinuous elements: Dirichlet data uex weakly (Nitsche) on ΓD, Neumann data
# ∇uex ⋅ n on the rest of ∂Ω.
function solve_weak_bc(dh, uex, ΓD; qr_order::Int, penalty = 20)
    grid = Ferrite.get_grid(dh)
    ip = only(ip for (_, ip) in bc_subdofhandlers(dh))
    shape = getrefshape(ip)
    geo = geometric_interpolation(getcelltype(grid))
    cv = CellValues(QuadratureRule{shape}(qr_order), ip, geo)
    fv = FacetValues(FacetQuadratureRule{shape}(qr_order), ip, geo)
    iv = InterfaceValues(FacetQuadratureRule{shape}(qr_order), ip, geo)
    topo = ExclusiveTopology(grid)
    K = allocate_matrix(dh; topology = topo, interface_coupling = trues(1, 1))
    F = zeros(ndofs(dh))
    γ = penalty * (Ferrite.getorder(ip) + 1)^2
    assemble_weak_bc!(start_assemble(K, F), dh, topo, cv, fv, iv, uex, ΓD, γ)
    return qr(K) \ F
end

# Function barrier for solve_weak_bc.
function assemble_weak_bc!(asm, dh, topo, cv, fv, iv, uex, ΓD, γ)
    f(x) = -tr(hessian(uex, x)) + uex(x)
    hK(coords) = maximum(norm(a - b) for a in coords, b in coords)
    nb = getnbasefunctions(cv)
    Ke = zeros(nb, nb)
    Fe = zeros(nb)
    for cc in CellIterator(dh)
        reinit!(cv, cc)
        fill!(Ke, 0)
        fill!(Fe, 0)
        for qp in 1:getnquadpoints(cv)
            dV = getdetJdV(cv, qp)
            fx = f(spatial_coordinate(cv, qp, getcoordinates(cc)))
            for i in 1:nb
                Ni, ∇Ni = shape_value(cv, qp, i), shape_gradient(cv, qp, i)
                Fe[i] += Ni * fx * dV
                for j in 1:nb
                    Ke[i, j] += (∇Ni ⋅ shape_gradient(cv, qp, j) + Ni * shape_value(cv, qp, j)) * dV
                end
            end
        end
        assemble!(asm, celldofs(cc), Ke, Fe)
    end
    Ki = zeros(2nb, 2nb)
    for ic in InterfaceIterator(dh, topo)
        reinit!(iv, ic)
        fill!(Ki, 0)
        μ = γ / hK(getcoordinates(ic.a))
        for qp in 1:getnquadpoints(iv)
            # Ferrite's jump is (there - here), the normal points out of "here".
            n = getnormal(iv, qp)
            dΓ = getdetJdV(iv, qp)
            for i in 1:(2nb)
                ji, ai = shape_value_jump(iv, qp, i), shape_gradient_average(iv, qp, i) ⋅ n
                for j in 1:(2nb)
                    jj, aj = shape_value_jump(iv, qp, j), shape_gradient_average(iv, qp, j) ⋅ n
                    Ki[i, j] += (ai * jj + ji * aj + μ * ji * jj) * dΓ
                end
            end
        end
        # Merge DOFs shared by both cells (EnrichedGalerkin's continuous part):
        # the assembler requires unique DOFs.
        idofs = interfacedofs(ic)
        udofs = unique(idofs)
        P = indexin(idofs, udofs)
        Ku = zeros(length(udofs), length(udofs))
        for i in eachindex(idofs), j in eachindex(idofs)
            Ku[P[i], P[j]] += Ki[i, j]
        end
        assemble!(asm, udofs, Ku)
    end
    for fc in FacetIterator(dh, getfacetset(Ferrite.get_grid(dh), "∂Ω"))
        reinit!(fv, fc)
        fill!(Ke, 0)
        fill!(Fe, 0)
        isD = FacetIndex(cellid(fc), Ferrite.getcurrentfacet(fv)) in ΓD
        μ = γ / hK(getcoordinates(fc))
        for qp in 1:getnquadpoints(fv)
            x = spatial_coordinate(fv, qp, getcoordinates(fc))
            n = getnormal(fv, qp)
            dΓ = getdetJdV(fv, qp)
            for i in 1:nb
                Ni, dNi = shape_value(fv, qp, i), shape_gradient(fv, qp, i) ⋅ n
                if isD
                    Fe[i] += (μ * Ni - dNi) * uex(x) * dΓ
                    for j in 1:nb
                        Nj, dNj = shape_value(fv, qp, j), shape_gradient(fv, qp, j) ⋅ n
                        Ke[i, j] += (μ * Ni * Nj - dNi * Nj - Ni * dNj) * dΓ
                    end
                else
                    Fe[i] += (gradient(uex, x) ⋅ n) * Ni * dΓ
                end
            end
        end
        assemble!(asm, celldofs(fc), Ke, Fe)
    end
    return
end

"""
    test_weak_bcs(ip, uex; n = 3, qr_order)

Weak boundary conditions for discontinuous elements, which have no boundary DOFs
for strong Dirichlet conditions (see [`solve_weak_bc`](@ref)): Dirichlet data
uex on the left half of ∂Ω, Neumann data on the rest. `uex` must lie in the
discrete space; the solution must then be exact.
"""
function test_weak_bcs(ip, uex; n = 3, qr_order = 2 * Ferrite.getorder(ip) + 2)
    return @testset "weak (SIPG) BCs: $ip" begin
        grid = bc_test_grid(bc_celltype(getrefshape(ip)); n)
        dh = bc_test_dofhandler(ip, grid)
        ΓD = weak_dirichlet_facets(grid)
        @test !isempty(ΓD) && length(ΓD) < length(getfacetset(grid, "∂Ω"))
        u = solve_weak_bc(dh, uex, ΓD; qr_order)
        @test max_domain_error(dh, u, uex) < 1.0e-10
        @test max_boundary_error(dh, u, uex; s = (c, f) -> FacetIndex(c, f) in ΓD ? range(0, 1, length = 9) : Float64[]) < 1.0e-10
    end
end

# Degree-k polynomial test solutions for the BC tests (in any dimension).
function bc_poly(k::Int)
    p1(x) = 1 + 2x[1] - x[end] / 2
    p2(x) = p1(x) + x[1] * x[end] - x[end]^2 / 3
    p3(x) = p2(x) + x[1]^3 / 3 - x[1] * x[end]^2
    return (p1, p2, p3)[k]
end

# ---------------------------------------------------------------------------
# Convergence
# ---------------------------------------------------------------------------
# The BC tests above use solutions in the discrete space. The convergence tests
# solve the same model problems for smooth solutions on two meshes and check
# the observed rates, which catches e.g. consistency errors of nonconforming
# elements that only show up for higher-degree solutions.

conv_u(x) = exp(x[1] / 2) * sin(2x[end] + 0.5) + (length(x) == 3 ? cos(x[2]) : 0.0)
conv_q(x) = Vec((sin(x[1] + 2x[2]), cos(2x[1] - x[2])))

# (‖u_h - uex‖, ‖D(u_h - uex)‖) in L2, with D the broken gradient (scalar) or
# div resp. curl (vector).
function error_norms(dh, u, uex; qr_order::Int)
    e0, e1 = 0.0, 0.0
    for (sdh, ip) in bc_subdofhandlers(dh)
        cv = CellValues(QuadratureRule{getrefshape(ip)}(qr_order), ip, geometric_interpolation(getcelltype(sdh)))
        D = if ip isa ScalarInterpolation
            identity
        elseif Ferrite.conformity(ip) isa Ferrite.HdivConformity
            tr
        else
            bc_curl
        end
        for cc in CellIterator(sdh)
            reinit!(cv, cc)
            ue = u[celldofs(cc)]
            for qp in 1:getnquadpoints(cv)
                x = spatial_coordinate(cv, qp, getcoordinates(cc))
                dV = getdetJdV(cv, qp)
                e0 += norm(function_value(cv, qp, ue) - uex(x))^2 * dV
                e1 += norm(D(function_gradient(cv, qp, ue)) - D(gradient(uex, x)))^2 * dV
            end
        end
    end
    return sqrt(e0), sqrt(e1)
end

"""
    test_convergence(ip_or_make_dh, rates; bc = :dirichlet, ns, qr_order, atol = 0.2)

Solve the model problem of the BC tests with a smooth exact solution on two
meshes with `ns = (n, 2n)` cells per direction and check that the observed
convergence rates of the L2 error and of the derivative error (broken H1
seminorm for scalar elements, div resp. curl for vector elements) are at least
`rates .- atol`. `bc` selects the problem: `:dirichlet` (strong, `Dirichlet`) or
`:neumann` ([`solve_scalar_bc`](@ref)), `:weak` ([`solve_weak_bc`](@ref), with
SIPG `penalty`), `:essential` (`ProjectedDirichlet`) or `:natural`
([`solve_vector_bc`](@ref)).
Instead of an interpolation, a function `n -> dh` can be passed (e.g. for mixed
meshes).
"""
function test_convergence(make_dh::Function, rates; bc::Symbol = :dirichlet, ns, qr_order::Int, atol = 0.2, penalty = 20, label = "")
    return @testset "convergence ($bc) $label" begin
        uex = bc in (:essential, :natural) ? conv_q : conv_u
        errs = map(ns) do n
            dh = make_dh(n)
            u = if bc in (:dirichlet, :neumann)
                solve_scalar_bc(dh, uex; dirichlet = bc == :dirichlet, qr_order)
            elseif bc == :weak
                solve_weak_bc(dh, uex, weak_dirichlet_facets(Ferrite.get_grid(dh)); qr_order, penalty)
            else
                solve_vector_bc(dh, uex; essential = bc == :essential, qr_order)
            end
            error_norms(dh, u, uex; qr_order)
        end
        observed = ntuple(i -> log2(errs[1][i] / errs[2][i]) / log2(ns[2] / ns[1]), 2)
        @test observed[1] > rates[1] - atol
        @test observed[2] > rates[2] - atol
    end
end
function test_convergence(ip::Interpolation, rates; ns = conv_meshes(getrefshape(ip)), qr_order = conv_qr_order(ip), kwargs...)
    CT = bc_celltype(getrefshape(ip))
    return test_convergence(n -> bc_test_dofhandler(ip, bc_test_grid(CT; n)), rates; ns, qr_order, label = string(ip), kwargs...)
end

conv_meshes(::Type{<:Ferrite.AbstractRefShape{1}}) = (8, 16)
conv_meshes(::Type{<:Ferrite.AbstractRefShape{2}}) = (4, 8)
conv_meshes(::Type{<:Ferrite.AbstractRefShape{3}}) = (2, 4)

# Ferrite's default quadrature rules go up to order 8 (triangles) resp. 5
# (tetrahedra).
conv_qr_order(ip) = min(2 * Ferrite.getorder(ip) + 2, getrefshape(ip) == RefTetrahedron ? 5 : 8)
