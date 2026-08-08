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
include("fortin_soulie.jl")
export FortinSoulie
include("crouzeix_falk.jl")
export CrouzeixFalk
include("conforming_crouzeix_raviart.jl")
export ConformingCrouzeixRaviart
include("bubble.jl")
export Bubble
include("dpc.jl")
export DPC
include("enriched_galerkin.jl")
export EnrichedGalerkin

# Blocked elements: documentation stubs describing required upstream Ferrite
# changes (see each file). Included but not exported.

end # module FerriteInterpolations
