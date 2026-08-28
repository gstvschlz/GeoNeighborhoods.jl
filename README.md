# GeoNeighborhoods.jl

[![CI](https://github.com/gstvschlz/GeoNeighborhoods.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/gstvschlz/GeoNeighborhoods.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/gstvschlz/GeoNeighborhoods.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gstvschlz/GeoNeighborhoods.jl)

Mining-style search neighborhoods for the
[GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) ecosystem: anisotropic
ellipsoids, angular sectors, sample quotas per sector and per drillhole,
half-space balancing, multipass searches, and explicit rules for refusing to
interpolate. The search is a `Meshes.BoundedNeighborSearchMethod`, so it
composes with the estimation machinery already in the ecosystem instead of
replacing it.

> [!NOTE]
> The package is tested but unregistered, and the public API is still settling.
> Every table below is generated output. The scripts in
> [`examples/`](examples/) are executed on each change, so the documented
> numbers and the tested behavior describe one implementation.

## Why a search neighborhood needs constraints

A block estimate is only as good as the samples that inform it. Down-hole
sampling is typically an order of magnitude denser than the spacing between
drillholes, so a plain "nearest N samples" search collapses onto whichever hole
passes closest to the block. What reaches the estimator is then a vertical
string of correlated composites standing in for a three-dimensional
neighborhood.

The cost of that collapse is measurable. Every example in this document uses
one synthetic dataset — a 5 × 5 grid of vertical holes on 45 m centers,
composited every 4 m — and estimates at the origin, which one hole passes
within 4 m of.

**Table 1.** Samples selected at the origin by two searches sharing an
ellipsoid of 95 × 95 × 62 m and a budget of twelve samples.

```
search              samples  holes  most from one hole  octants filled
------------------  -------  -----  ------------------  --------------
nearest 12          12       1      12                  1
sectors + hole cap  12       5      3                   8
```

The unconstrained search drew all twelve composites from a single hole and
filled one octant of the ellipsoid. Octant sectors combined with a cap of three
composites per hole spread the same budget of twelve samples over five holes
and eight octants. Neither the data, the ellipsoid, nor the sample budget
changed between the two rows.

Sample composition then propagates into the estimate itself.

**Table 2.** Ordinary kriging at the origin from each neighborhood of Table 1,
with a spherical variogram of range 120 m.

```
search              Au estimate
------------------  -----------
nearest 12          1.516
sectors + hole cap  1.501
```

Declustering moved the estimate by 0.015, or 1.0 % — away from the grade of the
one hole that had dominated it.
[Full example →](examples/01-why-sectors.md)

## Defining a neighborhood

A `SearchNeighborhood` collects the geometry and every constraint that governs
selection:

```julia
using GeoNeighborhoods
using GeoStatsBase: GslibAngles

search = SearchNeighborhood(
  (100, 50, 20),                     # semi-axes of the ellipsoid
  rotation = GslibAngles(30, 0, 0),  # or Datamine / Vulcan / Minesight / Rotations.jl

  sectors      = Octants(),
  minsamples   = 4,  maxsamples   = 24,
  minpersector = 1,  maxpersector = 3,
  minsectors   = 3,

  # refuse to interpolate rather than produce an indefensible estimate
  maxemptyconsecutive = 2,

  # force samples from both sides of a plane through the block center
  split = HalfSpace((0, 0, 1)),

  # at most 3 composites from any one hole, from at least 2 distinct holes
  category = CategoryRule(:BHID, maxper=3, mindistinct=2),
)
```

Sectors are assigned in the rotated, radius-normalized frame — the frame in
which the ellipsoid is a unit sphere — so anisotropy and sector geometry cannot
disagree. Rotating the search rotates the sectors with it. Two schemes are
available: `Octants()` splits the neighborhood by the sign of each rotated
coordinate, whereas `Azimuthal(n)` cuts `n` wedges about the vertical, halved
above and below on request. The distinction is critical for vertical holes,
because every composite in one such hole shares a single azimuthal wedge.
[How the schemes differ →](examples/02-sector-schemes.md)

A `CategoryRule` ties selection to a categorical column, typically the
drillhole identifier (BHID) or a domain code. Caps per value, minimum distinct
values, hard domain matching, and explicit per-value quotas are available, and
they compose:

```julia
category = [
  CategoryRule(:BHID, maxper=3, mindistinct=2),
  CategoryRule(:ROCK, match=:block),                  # only samples in the block's own domain
  CategoryRule(:ZONE, quotas=Dict("A" => 2:6, "B" => 0:4)),
]
```

[What each rule does →](examples/03-category-rules.md)

## Estimating

`interpolate` fits the model to the samples the search returns and asks it to
predict, leaving the model itself untouched:

```julia
using GeoStatsModels: Kriging
using GeoStatsFunctions: SphericalVariogram
using Unitful: @u_str

model = Kriging(SphericalVariogram(range=120.0u"m"))
estimate = interpolate(samples, blocks, model; search, vars=(:Au,))
```

A location no neighborhood can serve is left `missing` rather than filled with
a number nobody can defend. Setting `diagnostics=true` appends the search
outcome to every location, so a gap in the model becomes explainable rather
than merely visible.

**Table 3.** Audit trail for the first two blocks of a 144-block model
estimated under octant sectors and a per-hole cap, reporting the samples,
sectors, and holes each estimate drew on.

```
block  Au     samples  sectors  holes  reason
-----  -----  -------  -------  -----  --------
1      1.438  12       8        4      Accepted
2      1.539  12       8        4      Accepted
```

**Table 4.** Search outcome across the same 144 blocks under a stricter
neighborhood, which demands eight samples from at least six octants and three
distinct holes.

```
outcome        blocks
-------------  ------
TooFewSectors  100
Accepted       44
```

The 100 rejected blocks are not failures of the estimator. Each one names the
constraint that could not be met, which is the information a resource
estimator needs to decide whether to relax the neighborhood or to leave the
ground uninformed.

## Multipass

A single neighborhood forces one compromise on the whole model: tight enough to
be defensible near data, and every block further out goes uninterpolated. A
`MultiPass` states the compromise explicitly instead. Each pass is a complete
neighborhood with its own quotas and rejection rules rather than merely a wider
radius, and the first pass whose constraints hold serves the block.

```julia
search = MultiPass(
  SearchNeighborhood(( 35,  35, 35), sectors=Octants(), minsamples=8, minsectors=4),
  SearchNeighborhood(( 70,  70, 55), sectors=Octants(), minsamples=5, minsectors=2),
  SearchNeighborhood((140, 140, 90), minsamples=2),
)
```

**Table 5.** Blocks served by each pass of that search over a 676-block model.

```
pass  blocks  % of model
----  ------  ----------
1     36      5.3
2     196     29.0
3     236     34.9

estimated:     468 of 676
not estimated: 208
```

The output records which pass served each block, so the model carries its own
confidence ordering.
[Full example →](examples/04-multipass.md)

## Verification against the ecosystem

Kriging never sees the neighborhood: it receives a subset of samples and
solves. This package changes which samples go in, and nothing else. Switching
the sectors and quotas off therefore reduces the search to "every sample inside
the ellipsoid, nearest first" — precisely what `Meshes.KBallSearch` supplies to
`GeoStatsModels.fitpredict` — and the two must then return the same estimate.

Two checks confirm that they do. The test suite compares the searcher against
`KBallSearch` on the unconstrained case, and `interpolate` reproduces
`fitpredict` across a 144-block model:

```
blocks:              144
estimated by both:   144
same missing blocks: true
largest difference:  8.881784197001252e-16
```

A largest difference of 9 × 10⁻¹⁶ is floating-point noise rather than agreement
to within a tolerance. Once the neighborhood does something, the estimates
diverge, which is the entire purpose of the package.
[Full example →](examples/05-kriging.md)

## Limitations

- The package is not registered, and the public API may still change.
- `interpolate` runs its own prediction loop over the public `fit` and
  `predict` interface, because `fitpredict` constructs its searcher internally
  from a `MetricBall` and accepts no searcher of its own. A `searcher=` keyword
  upstream would remove the duplication.
- Unlike `fitpredict`, `interpolate` does not rescale the data, domain, and
  model before solving, so the results depend only on the coordinates supplied.
- Every number documented here comes from one synthetic dataset. The behavior
  is tested, but no production estimation run is reported.

## Installation

The package is not registered, so install it from the repository:

```julia
using Pkg; Pkg.add(url="https://github.com/gstvschlz/GeoNeighborhoods.jl")
```

## Development

The project uses [mise](https://mise.jdx.dev) for tooling:

```bash
mise install         # julia 1.11.9
mise run instantiate
mise run test
```

Regenerate the examples after changing the library:

```bash
mise run examples
```

## License

MIT
