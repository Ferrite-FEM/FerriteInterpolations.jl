module FerriteInterpolations

using Ferrite:
    Ferrite, ScalarInterpolation, VectorInterpolation,
    RefLine, RefTriangle, RefQuadrilateral, RefTetrahedron, RefHexahedron,
    RefPrism, RefPyramid
using Tensors: Vec

include("common.jl")
include("l2_piola.jl")
export L2PiolaMapping

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
include("gauss_legendre.jl")
export GaussLegendre
include("radau.jl")
export Radau
include("bdm.jl")
export BDM
include("nedelec2.jl")
export NedelecSecondKind
include("bdfm.jl")
export BDFM
include("trimmed_serendipity_div.jl")
export TrimmedSerendipityDiv
include("trimmed_serendipity_curl.jl")
export TrimmedSerendipityCurl
include("tnt_div.jl")
export TNTDiv
include("tnt_curl.jl")
export TNTCurl
include("taylor.jl")
export Taylor, TaylorMapping
include("transition.jl")
export Transition

# Blocked elements: documentation stubs describing required upstream Ferrite
# changes (see each file). Included but not exported.
include("mardal_tai_winther.jl")
include("huang_zhang.jl")

end # module FerriteInterpolations
