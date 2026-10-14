using Ferrite
using Test
include("test_utils.jl")

# symfem cross-check of Ferrite's own Lagrange family (Lagrange,
# DiscontinuousLagrange and their vectorizations). Ferrite's test suite covers
# the interpolation properties and integration; this file only checks the
# bases against the DefElement reference implementation.

@testset "Ferrite Lagrange (symfem)" begin
    lagrange_ips = (
        Lagrange{RefLine, 1}(), Lagrange{RefLine, 2}(),
        Lagrange{RefTriangle, 1}(), Lagrange{RefTriangle, 2}(), Lagrange{RefTriangle, 3}(),
        Lagrange{RefTriangle, 4}(), Lagrange{RefTriangle, 5}(),
        Lagrange{RefQuadrilateral, 1}(), Lagrange{RefQuadrilateral, 2}(), Lagrange{RefQuadrilateral, 3}(),
        Lagrange{RefTetrahedron, 1}(), Lagrange{RefTetrahedron, 2}(), Lagrange{RefTetrahedron, 3}(),
        Lagrange{RefTetrahedron, 4}(),
        Lagrange{RefHexahedron, 1}(), Lagrange{RefHexahedron, 2}(), Lagrange{RefHexahedron, 3}(),
        Lagrange{RefPrism, 1}(), Lagrange{RefPrism, 2}(),
        Lagrange{RefPyramid, 1}(), Lagrange{RefPyramid, 2}(),
    )

    # Ferrite >= 1.8 places the nodes of order >= 3 Lagrange at Gauss-Lobatto-Legendre
    # (GLL) points: tensor products of the 1D GLL points on hypercubes and the GLL-warped
    # lattice of Basix on simplices. symfem only implements equispaced nodes on simplices,
    # and its GLL variant on hypercubes is prohibitively slow (symbolic basis with radicals),
    # so for these elements check that the basis spans the same space as symfem's
    # (equispaced) Lagrange element and that it is nodal (Kronecker delta property) instead.
    # Together with the node checks in Ferrite's own tests this determines the basis.
    # The GLL nodes are detected from the loaded Ferrite rather than from the version so that
    # the tests also work against Ferrite master (which has the same version number as the
    # latest release).
    ferrite_gll_nodes = isdefined(Ferrite, :_gll_points)
    symfem_has_basis(ip) = !ferrite_gll_nodes || Ferrite.getorder(ip) <= 2

    # The Ferrite <-> symfem permutation is derived geometrically from the
    # DOF points (the nodal bases agree exactly, including on the pyramid).
    for ip in lagrange_ips
        k = Ferrite.getorder(ip)
        if symfem_has_basis(ip)
            test_symfem_reference(ip, "Lagrange", k, symfem_point_perm(ip, "Lagrange", k))
        else
            test_symfem_span(ip, "Lagrange", k)
            test_kronecker_delta(ip)
        end
    end

    # DiscontinuousLagrange: same basis as Lagrange, plus degree 0 on every cell.
    for lag in lagrange_ips
        shape, k = getrefshape(lag), Ferrite.getorder(lag)
        ip = DiscontinuousLagrange{shape, k}()
        if symfem_has_basis(lag)
            test_symfem_reference(ip, "discontinuous Lagrange", k, symfem_point_perm(ip, "discontinuous Lagrange", k))
        else
            test_symfem_span(ip, "discontinuous Lagrange", k)
            test_kronecker_delta(ip)
        end
    end
    for shape in (RefLine, RefTriangle, RefQuadrilateral, RefTetrahedron, RefHexahedron, RefPrism, RefPyramid)
        test_symfem_reference(DiscontinuousLagrange{shape, 0}(), "discontinuous Lagrange", 0, [0])
    end

    # Vectorized Lagrange ("vector Lagrange" on simplices, "vector Q" on
    # hypercubes). Both Ferrite and symfem order the DOFs node by node with
    # the components innermost, so the permutation is the scalar one blown up
    # per component.
    for (ip, family) in (
            (Lagrange{RefTriangle, 2}()^2, "vector Lagrange"),
            (Lagrange{RefTetrahedron, 2}()^3, "vector Lagrange"),
            (Lagrange{RefQuadrilateral, 2}()^2, "vector Q"),
            (Lagrange{RefHexahedron, 2}()^3, "vector Q"),
        )
        sip = ip.ip
        k = Ferrite.getorder(sip)
        vdim = Ferrite.n_components(ip)
        p = symfem_point_perm(sip, "Lagrange", k)
        sperm = [(1, vdim * p[n] + c - 1) for n in eachindex(p) for c in 1:vdim]
        test_symfem_reference_vector(ip, family, k, sperm)
    end
end
