# Elements already implemented in Ferrite.jl

These DefElement entries are **not** re-implemented in this package because
Ferrite.jl (v1.6) already provides them, possibly for a subset of the degrees
and reference cells that DefElement defines. Extending their coverage belongs
upstream in Ferrite — except where a note below says a gap is filled here under
a different type name.

Degrees are given in *Ferrite's* numbering. For Raviart–Thomas and Nédélec
(first kind), DefElement degree `k` corresponds to Ferrite degree `k + 1`.

| DefElement entry | Ferrite type | Coverage in Ferrite v1.6 | Gaps |
|---|---|---|---|
| [lagrange](https://defelement.org/elements/lagrange.html) | `Lagrange` | line 1–2, triangle 1–5, quadrilateral 1–3, tetrahedron 1–4, hexahedron 1–3, prism 1–2, pyramid 1–2 | higher degrees on all cells (DefElement is unbounded) |
| [lagrange (discontinuous)](https://defelement.org/elements/lagrange.html) | `DiscontinuousLagrange` | as `Lagrange`, plus degree 0 on every cell | as `Lagrange` |
| [vector-lagrange](https://defelement.org/elements/vector-lagrange.html) | `Lagrange^vdim` | via vectorization of the scalar element | follows `Lagrange` |
| [vector-q](https://defelement.org/elements/vector-q.html) | `Lagrange^vdim` | via vectorization of the scalar element | follows `Lagrange` |
| [crouzeix-raviart](https://defelement.org/elements/crouzeix-raviart.html) | `CrouzeixRaviart` | triangle 1, tetrahedron 1 | higher odd degrees |
| [crouzeix-raviart](https://defelement.org/elements/crouzeix-raviart.html) (quad/hex) | `RannacherTurek` | quadrilateral 1, hexahedron 1 | — |
| [raviart-thomas](https://defelement.org/elements/raviart-thomas.html) | `RaviartThomas` | triangle 1–2, quadrilateral 1, tetrahedron 1, hexahedron 1 | triangle ≥3, quadrilateral ≥2; tet/hex ≥2 additionally blocked upstream (>1 DOF per face needs face-orientation DOF transformations) |
| [nedelec1](https://defelement.org/elements/nedelec1.html) | `Nedelec` | triangle 1–2, quadrilateral 1, tetrahedron 1, hexahedron 1 | same situation as Raviart–Thomas; prism missing entirely |
| [brezzi-douglas-marini](https://defelement.org/elements/brezzi-douglas-marini.html) | `BrezziDouglasMarini` | triangle 1 | triangle ≥2 (planned in this package as `BDM` to avoid the name clash); tetrahedron blocked upstream (face-orientation DOF transformations) |
| [serendipity](https://defelement.org/elements/serendipity.html) | `Serendipity` | quadrilateral 2, hexahedron 2 | degree ≥3 (edge DOFs become integral moments; edge-reversal permutation no longer applies — needs upstream DOF transformations). Note: Ferrite uses edge-midpoint point evaluations instead of DefElement's integral moments at degree 2 (same span, different basis). |
| [bubble-enriched-lagrange](https://defelement.org/elements/bubble-enriched-lagrange.html) | `BubbleEnrichedLagrange` | triangle 1 | DefElement degree 2 |
| [vector-bubble-enriched-lagrange](https://defelement.org/elements/vector-bubble-enriched-lagrange.html) | `BubbleEnrichedLagrange^2` | via vectorization | follows the scalar element |
