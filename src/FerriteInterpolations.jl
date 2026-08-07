module FerriteInterpolations

using Ferrite:
    Ferrite, ScalarInterpolation, VectorInterpolation,
    RefLine, RefTriangle, RefQuadrilateral, RefTetrahedron, RefHexahedron,
    RefPrism, RefPyramid
using Tensors: Vec

include("common.jl")

# Implemented elements (one file per element)
include("bernstein.jl")
export Bernstein

# Blocked elements: documentation stubs describing required upstream Ferrite
# changes (see each file). Included but not exported.

end # module FerriteInterpolations
