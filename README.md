# FerriteInterpolations.jl

Finite element interpolations from the [DefElement encyclopedia](https://defelement.org)
for use with [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl).

This package implements elements that Ferrite does not (yet) provide itself.
Elements that already exist in Ferrite are deliberately **not** re-implemented
here — see [ExistingElements.md](ExistingElements.md) for that list, including
notes on degrees/cells Ferrite covers and gaps.

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
