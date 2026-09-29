# FerriteInterpolations.jl

Finite element interpolations from the [DefElement encyclopedia](https://defelement.org)
for use with [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl).

The implementations are, to a large extent, automated/robot translations of
the element tabulations in [DefElement](https://defelement.org) (through its
reference implementation [symfem](https://github.com/mscroggs/symfem)) into
Ferrite's conventions: reference cells, entity numbering, DOF ordering, and
orientation handling. Every implemented basis is cross-checked against symfem
in the test suite.

This package implements elements that Ferrite does not (yet) provide itself.
Elements that already exist in Ferrite are deliberately **not** re-implemented
here — see [Implemented in Ferrite](#implemented-in-ferrite) for that list, including
notes on degrees/cells Ferrite covers and gaps.

## Elements

Element names link to their DefElement entries. The Type column gives the
interpolation type to use; "-" means the element is not available.

### Implemented

Interpolations implemented in this package.

| Element | Type | Notes |
|---|---|---|
| [bernstein](https://defelement.org/elements/bernstein.html) | `Bernstein` | line 1–3, triangle 1–3, tetrahedron 1–2 (tet ≥ 3 needs [#5](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/5)) |
| [brezzi-douglas-fortin-marini](https://defelement.org/elements/brezzi-douglas-fortin-marini.html) | `BrezziDouglasFortinMarini` | triangle, degrees 1–2 |
| [brezzi-douglas-marini](https://defelement.org/elements/brezzi-douglas-marini.html) | `BDM` | triangle, degree 2 (degree 1 is Ferrite's `BrezziDouglasMarini`, see [below](#implemented-in-ferrite)) |
| [bubble](https://defelement.org/elements/bubble.html) | `Bubble` | line 2–3, triangle 3–4, tetrahedron 4 |
| [conforming-crouzeix-raviart](https://defelement.org/elements/conforming-crouzeix-raviart.html) | `ConformingCrouzeixRaviart` | triangle, degrees 2–4 (degree 1 coincides with P1 Lagrange) |
| [crouzeix-falk](https://defelement.org/elements/crouzeix-falk.html) | `CrouzeixFalk` | triangle (degree 3, as defined) |
| [dPc](https://defelement.org/elements/dpc.html) | `DPC` | quadrilateral 1–3, hexahedron 1–2, with the L2 Piola mapping; on the interval it coincides with `DiscontinuousLagrange` |
| [enriched-galerkin](https://defelement.org/elements/enriched-galerkin.html) | `EnrichedGalerkin` | every cell/degree Ferrite's `Lagrange` supports |
| [fortin-soulie](https://defelement.org/elements/fortin-soulie.html) | `FortinSoulie` | triangle (degree 2), as in the original paper: continuous P2 plus one neutral function per cell, continuous at the Gauss–Legendre points of each edge. Differs from DefElement's definition (edge points 1/3 and 2/3, asymmetric local DOFs), see the file header and [#33](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/33) |
| [gauss-legendre](https://defelement.org/elements/gauss-legendre.html) | `GaussLegendre` | line/quadrilateral/hexahedron, any degree (modal basis, following symfem) |
| [nedelec2](https://defelement.org/elements/nedelec2.html) | `NedelecSecondKind` | triangle, degrees 1–2 |
| [radau](https://defelement.org/elements/radau.html) | `Radau` | interval degree 2 (degree ≥ 3 has no symfem reference; quad/hex need [#2](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/2)) |
| [taylor](https://defelement.org/elements/taylor.html) | `Taylor` | line/triangle, degrees 1–2, with the cell-dependent `TaylorMapping` |
| [tnt-curl](https://defelement.org/elements/tnt-curl.html) | `TNTCurl` | quadrilateral, degree 1 |
| [tnt-div](https://defelement.org/elements/tnt-div.html) | `TNTDiv` | quadrilateral, degree 1 |
| [transition](https://defelement.org/elements/transition.html) | `Transition` | triangle, interior order 2, per-edge orders 1–2 (tuple type parameter) |
| [trimmed-serendipity-curl](https://defelement.org/elements/trimmed-serendipity-curl.html) | `TrimmedSerendipityCurl` | quadrilateral, degree 1 |
| [trimmed-serendipity-div](https://defelement.org/elements/trimmed-serendipity-div.html) | `TrimmedSerendipityDiv` | quadrilateral, degree 1 |
| [vector dPc](https://defelement.org/elements/vector-dpc.html) | `DPC^vdim` | vectorized `DPC` (`VectorizedInterpolation`) |

### Blocked

Elements that cannot currently be implemented against Ferrite (each linked
issue describes the required upstream changes), or that are out of scope.

| Element | Type | Notes |
|---|---|---|
| [alfeld-sorokina](https://defelement.org/elements/alfeld-sorokina.html) | - | blocked (macro split + derivative DOFs): [#21](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/21) |
| [argyris](https://defelement.org/elements/argyris.html) | - | blocked (derivative DOFs): [#7](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/7) |
| [arnold-boffi-falk](https://defelement.org/elements/arnold-boffi-falk.html) | - | blocked (Kirby mapping): [#14](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/14) |
| [arnold-winther](https://defelement.org/elements/arnold-winther.html) | - | blocked (matrix-valued): [#25](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/25) |
| [bell](https://defelement.org/elements/bell.html) | - | blocked (derivative DOFs): [#8](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/8) |
| [bernardi-raugel](https://defelement.org/elements/bernardi-raugel.html) | - | blocked (geometry-dependent facet DOFs): [#13](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/13) |
| [bogner-fox-schmitt](https://defelement.org/elements/bogner-fox-schmitt.html) | - | blocked (derivative DOFs): [#12](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/12) |
| [brezzi-douglas-duran-fortin](https://defelement.org/elements/brezzi-douglas-duran-fortin.html) | - | blocked (3D face DOF orientation): [#28](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/28) |
| [buffa-christiansen](https://defelement.org/elements/buffa-christiansen.html) | - | out of scope (dual polygon element) |
| [direct-serendipity](https://defelement.org/elements/direct-serendipity.html) | - | blocked (basis built on the physical cell): [#29](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/29) |
| [dual](https://defelement.org/elements/dual.html) | - | out of scope (dual polygon element) |
| [gopalakrishnan-lederer-schoberl](https://defelement.org/elements/gopalakrishnan-lederer-schoberl.html) | - | blocked (matrix-valued): [#27](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/27) |
| [guzman-neilan (first kind)](https://defelement.org/elements/guzman-neilan.html) | - | blocked (macro split): [#22](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/22) |
| [guzman-neilan (second kind)](https://defelement.org/elements/guzman-neilan2.html) | - | blocked (macro split): [#22](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/22) |
| [hellan-herrmann-johnson](https://defelement.org/elements/hellan-herrmann-johnson.html) | - | blocked (matrix-valued): [#26](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/26) |
| [hermite](https://defelement.org/elements/hermite.html) | - | blocked (derivative DOFs): [#6](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/6) |
| [hsieh-clough-tocher](https://defelement.org/elements/hsieh-clough-tocher.html) | - | blocked (macro split): [#18](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/18) |
| [huang-zhang](https://defelement.org/elements/huang-zhang.html) | - | blocked (tangential edge moments): [#16](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/16) |
| [johnson-mercier](https://defelement.org/elements/johnson-mercier.html) | - | blocked (macro split, matrix-valued): [#23](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/23) |
| [mardal-tai-winther](https://defelement.org/elements/mardal-tai-winther.html) | - | blocked (tangential edge moments): [#15](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/15) |
| [morley](https://defelement.org/elements/morley.html) | - | blocked (normal-derivative DOFs): [#9](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/9) |
| [morley-wang-xu](https://defelement.org/elements/morley-wang-xu.html) | - | blocked (normal-derivative DOFs): [#10](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/10) |
| [nonconforming-arnold-winther](https://defelement.org/elements/nonconforming-arnold-winther.html) | - | blocked (matrix-valued): [#25](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/25) |
| [p1-iso-p2](https://defelement.org/elements/p1-iso-p2.html) | - | blocked (macro split): [#19](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/19) |
| [p1-macro](https://defelement.org/elements/p1-macro.html) | - | blocked (macro split): [#20](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/20) |
| [reduced-hsieh-clough-tocher](https://defelement.org/elements/reduced-hsieh-clough-tocher.html) | - | blocked (macro split): [#18](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/18) |
| [regge](https://defelement.org/elements/regge.html) | - | blocked (matrix-valued): [#24](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/24) |
| [rotated-buffa-christiansen](https://defelement.org/elements/rotated-buffa-christiansen.html) | - | out of scope (dual polygon element) |
| [tnt (scalar)](https://defelement.org/elements/tnt.html) | - | blocked (sign-flipping edge moment weights): [#17](https://github.com/Ferrite-FEM/FerriteInterpolations.jl/issues/17) |

### Implemented in Ferrite

These DefElement entries are **not** re-implemented in this package because
Ferrite.jl (v1.6) already provides them, possibly for a subset of the degrees
and reference cells that DefElement defines. Extending their coverage belongs
upstream in Ferrite — except where a note says a gap is filled here under a
different type name.

Degrees are given in *Ferrite's* numbering. For Raviart–Thomas and Nédélec
(first kind), DefElement degree `k` corresponds to Ferrite degree `k + 1`.

| Element | Type | Notes |
|---|---|---|
| [brezzi-douglas-marini](https://defelement.org/elements/brezzi-douglas-marini.html) | `BrezziDouglasMarini` | triangle 1. Degree 2 is `BDM` in this package (distinct name to avoid clashing with Ferrite's type); triangle ≥ 3 missing; tetrahedron blocked upstream (face-orientation DOF transformations) |
| [bubble-enriched-lagrange](https://defelement.org/elements/bubble-enriched-lagrange.html) | `BubbleEnrichedLagrange` | triangle 1; DefElement degree 2 missing |
| [crouzeix-raviart](https://defelement.org/elements/crouzeix-raviart.html) | `CrouzeixRaviart` | triangle 1, tetrahedron 1; higher odd degrees missing |
| [crouzeix-raviart](https://defelement.org/elements/crouzeix-raviart.html) (quad/hex) | `RannacherTurek` | quadrilateral 1, hexahedron 1 |
| [lagrange](https://defelement.org/elements/lagrange.html) | `Lagrange` | line 1–2, triangle 1–5, quadrilateral 1–3, tetrahedron 1–4, hexahedron 1–3, prism 1–2, pyramid 1–2; higher degrees missing on all cells (DefElement is unbounded) |
| [lagrange](https://defelement.org/elements/lagrange.html) (discontinuous) | `DiscontinuousLagrange` | as `Lagrange`, plus degree 0 on every cell |
| [nedelec1](https://defelement.org/elements/nedelec1.html) | `Nedelec` | triangle 1–2, quadrilateral 1, tetrahedron 1, hexahedron 1; same gaps as `RaviartThomas`, prism missing entirely |
| [raviart-thomas](https://defelement.org/elements/raviart-thomas.html) | `RaviartThomas` | triangle 1–2, quadrilateral 1, tetrahedron 1, hexahedron 1; triangle ≥ 3 and quadrilateral ≥ 2 missing; tet/hex ≥ 2 additionally blocked upstream (> 1 DOF per face needs face-orientation DOF transformations) |
| [serendipity](https://defelement.org/elements/serendipity.html) | `Serendipity` | quadrilateral 2, hexahedron 2. Ferrite uses edge-midpoint point evaluations instead of DefElement's integral moments (same span, different basis). Degree ≥ 3 missing (edge DOFs become integral moments; the edge-reversal permutation no longer applies — needs upstream DOF transformations) |
| [vector-bubble-enriched-lagrange](https://defelement.org/elements/vector-bubble-enriched-lagrange.html) | `BubbleEnrichedLagrange^2` | via vectorization; follows the scalar element |
| [vector-lagrange](https://defelement.org/elements/vector-lagrange.html) | `Lagrange^vdim` | via vectorization; follows `Lagrange` |
| [vector-q](https://defelement.org/elements/vector-q.html) | `Lagrange^vdim` | via vectorization; follows `Lagrange` |

## Structure

- One source file per element in `src/`, defining an interpolation type in
  Ferrite's style (e.g. `Bernstein{RefTriangle, 2}()`).
- Elements that cannot currently be implemented against Ferrite (they need
  upstream changes, e.g. general DOF transformations, matrix-valued
  interpolations, or macro-element support) also get one file each, containing
  a description of the required upstream changes and a placeholder type whose
  constructor throws. These are not exported.
- One test file per element in `test/`, with interpolation-level tests
  (properties, reference tabulation tables) and integration tests via
  `CellValues`/`FacetValues`.
- Elements that exist in Ferrite also get one test file each, with symfem
  cross-checks only (Ferrite's own test suite covers the rest). Their bases
  agree with symfem up to a signed permutation (scaled by the Piola factor of
  the reference-cell map for `RaviartThomas`/`Nedelec` on hypercubes), except
  `Serendipity`, where only the spans agree since Ferrite uses edge-midpoint
  point evaluations instead of integral moments.

## Development setup

The package is developed against a local Ferrite checkout (sibling directory
`../Ferrite`):

```
julia --project=. -e 'using Pkg; Pkg.develop(path="../Ferrite"); Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

`Manifest.toml` is intentionally not committed.

## Caveats

- Ferrite's interpolation interface is unexported and documented as not fully
  stable; this package implements `Ferrite.`-qualified internal methods and
  pins `Ferrite = "1.6"` in `[compat]` as a guard.
- The test suite cross-checks every basis against
  [symfem](https://github.com/mscroggs/symfem) (DefElement's reference
  implementation), imported live via PythonCall — a test-only dependency;
  `test/CondaPkg.toml` provisions Python and symfem automatically. Note that
  symfem/DefElement and Ferrite differ in reference cell coordinates and in
  entity numbering/ordering — `test_symfem_reference` maps the coordinates,
  and each element's test file documents its DOF permutation.
