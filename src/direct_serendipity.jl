# Direct serendipity element (https://defelement.org/elements/direct-serendipity.html,
# DefElement: elements/direct-serendipity.def).
#
# Direct serendipity element.
#
# STATUS: BLOCKED -- cannot be expressed as a reference element
#
# The direct serendipity element's shape functions are RATIONAL functions constructed directly on the physical quadrilateral (they depend on the physical vertex positions to retain full approximation order on non-parallelogram cells). There is no factorization into a fixed reference basis composed
# with a (geometry-dependent) mapping, which is the only model Ferrite (and
# this package) supports; the element would need to be constructed per
# physical cell at reinit! time, essentially a different finite element
# framework feature.

"""
    DirectSerendipity

Placeholder: not implementable in Ferrite today. See the `STATUS: BLOCKED`
notes in `src/direct_serendipity.jl`.
"""
struct DirectSerendipity
    function DirectSerendipity()
        return error("DirectSerendipity is not implemented; see the STATUS: BLOCKED notes in src/direct_serendipity.jl.")
    end
end
