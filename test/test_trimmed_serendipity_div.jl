using FerriteInterpolations
using Ferrite
using Test
include("test_utils.jl")

@testset "TrimmedSerendipityDiv" begin
    ip = TrimmedSerendipityDiv{RefQuadrilateral, 1}()

    # (i) Interpolation-level tests
    test_interpolation_properties(ip)
    test_type_genericity(ip)

    # Live symfem cross-check ("TSdiv"): signed permutation derived with
    # sympy from the Ferrite DOF conventions; the factor 1/2 (contravariant)
    # is the Piola factor of the [0,1]^2 -> [-1,1]^2 reference map (see the
    # source file).
    test_symfem_reference_vector(
        ip, "TSdiv", 1,
        [(-1//2, 0), (-1//2, 1), (-1//2, 4), (-1//2, 5), (1//2, 7), (1//2, 6), (1//2, 3), (1//2, 2), (1//2, 8), (1//2, 9)],
    )

    # (ii) Integration tests on a two-cell quadrilateral mesh whose shared
    # edge (2,5)/(5,2) has opposite local orientation.
    nodes = [Node(Vec((x, y))) for y in (0.0, 1.1) for x in (0.0, 1.0, 2.1)]
    grid = Grid([Quadrilateral((1, 2, 5, 4)), Quadrilateral((2, 3, 6, 5))], nodes)
    test_hdiv_two_cell(ip, grid, 2, 4, 2)
    test_divergence_theorem(ip, getcells(grid, 1), getcoordinates(grid, 1))
end
