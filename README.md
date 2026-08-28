# GeoNeighborhoods.jl

[![CI](https://github.com/gstvschlz/GeoNeighborhoods.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/gstvschlz/GeoNeighborhoods.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/gstvschlz/GeoNeighborhoods.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gstvschlz/GeoNeighborhoods.jl)

Mining-style search neighbourhoods for the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl)
ecosystem — anisotropic ellipsoids, angular sectors, sample quotas per sector
and per drillhole, half-space balancing, and multipass searches.

> [!NOTE]
> Working but unregistered, and the API may still move. Everything shown below
> is executed on every change: the tables are real output, produced by the
> scripts in [`examples/`](examples/).

## Why

A block estimate is only as good as the samples that inform it. Down-hole
sampling is typically an order of magnitude denser than the spacing between
holes, so a plain "nearest N samples" search collapses onto whichever hole
happens to pass closest — a vertical string of correlated composites standing in
for a three-dimensional neighbourhood.

Twelve samples, one ellipsoid of 95 × 95 × 62 m, a 5 × 5 grid of vertical holes
on 45 m centres:

```
search              samples  holes  most from one hole  octants filled
------------------  -------  -----  ------------------  --------------
nearest 12          12       1      12                  1
sectors + hole cap  12       5      3                   8
```

Same data, same ellipsoid, same sample budget — a completely different
neighbourhood, and a different answer from the same kriging model:

```
search              Au estimate
------------------  -----------
nearest 12          1.516
sectors + hole cap  1.501
```

[Full example →](examples/01-why-sectors.md)

## It works with the kriging you already have

Kriging never sees the neighbourhood; it receives a subset of samples and
solves. This package changes *which samples go in*, and nothing else — so with
the constraints switched off it reproduces the ecosystem's own `fitpredict`:

```
blocks:              144
estimated by both:   144
same missing blocks: true
largest difference:  8.881784197001252e-16
```

[Full example →](examples/05-kriging.md)

`NeighborhoodSearch` is a `Meshes.BoundedNeighborSearchMethod`, so it also drops
into anything else in the ecosystem that accepts a searcher.

## Defining a neighbourhood

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

  # force samples from both sides of a plane through the block centre
  split = HalfSpace((0, 0, 1)),

  # at most 3 composites from any one hole, from at least 2 distinct holes
  category = CategoryRule(:BHID, maxper=3, mindistinct=2),
)
```

Sectors are measured in the **rotated, radius-normalised frame**, so anisotropy
and sector geometry can never disagree. Two conventions are available —
`Octants()` splits by the sign of each rotated coordinate, `Azimuthal(n)` cuts
`n` wedges about the vertical, optionally halved above and below.
[How they differ →](examples/02-sector-schemes.md)

Category rules cover caps per hole, minimum distinct holes, hard domain
boundaries and explicit per-value quotas:

```julia
category = [
  CategoryRule(:BHID, maxper=3, mindistinct=2),
  CategoryRule(:ROCK, match=:block),                  # only samples in the block's own domain
  CategoryRule(:ZONE, quotas=Dict("A" => 2:6, "B" => 0:4)),
]
```

[What each one does →](examples/03-category-rules.md)

## Estimating

```julia
using GeoStatsModels: Kriging
using GeoStatsFunctions: SphericalVariogram

model = Kriging(SphericalVariogram(range=120.0u"m"))
estimate = interpolate(samples, blocks, model; search, vars=(:Au,))
```

A location no neighbourhood can serve is left `missing` rather than filled with
a number nobody can defend. With `diagnostics=true` every location also carries
the reason:

```
block  Au     samples  sectors  holes  reason
-----  -----  -------  -------  -----  --------
1      1.438  12       8        4      Accepted
2      1.539  12       8        4      Accepted
```

```
outcome        blocks
-------------  ------
TooFewSectors  100
Accepted       44
```

## Multipass

Each pass is a complete neighbourhood, not merely a wider radius. The first
whose constraints hold wins, and the output records which:

```julia
search = MultiPass(
  SearchNeighborhood(( 35,  35, 35), sectors=Octants(), minsamples=8, minsectors=4),
  SearchNeighborhood(( 70,  70, 55), sectors=Octants(), minsamples=5, minsectors=2),
  SearchNeighborhood((140, 140, 90), minsamples=2),
)
```

```
pass  blocks  % of model
----  ------  ----------
1     36      5.3
2     196     29.0
3     236     34.9

estimated:     468 of 676
not estimated: 208
```

[Full example →](examples/04-multipass.md)

## Install

Not yet registered.

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

To regenerate the examples after changing the library:

```bash
julia --project=examples examples/generate.jl
```

## License

MIT
