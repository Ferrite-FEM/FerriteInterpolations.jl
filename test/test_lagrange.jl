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

    # The Ferrite <-> symfem permutation is derived geometrically from the
    # DOF points (the nodal bases agree exactly, including on the pyramid).
    for ip in lagrange_ips
        k = Ferrite.getorder(ip)
        test_symfem_reference(ip, "Lagrange", k, symfem_point_perm(ip, "Lagrange", k))
    end

    # DiscontinuousLagrange: same basis as Lagrange, plus degree 0 on every cell.
    for lag in lagrange_ips
        shape, k = getrefshape(lag), Ferrite.getorder(lag)
        ip = DiscontinuousLagrange{shape, k}()
        test_symfem_reference(ip, "discontinuous Lagrange", k, symfem_point_perm(ip, "discontinuous Lagrange", k))
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
