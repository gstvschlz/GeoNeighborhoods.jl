# Why sector search

A block estimate is only as good as the samples that inform it. Down-hole
sampling is typically an order of magnitude denser than the spacing between
holes, so a plain "nearest N samples" search collapses onto whichever hole
happens to pass closest — a vertical string of correlated composites
standing in for a three-dimensional neighbourhood.

The data here is a 5 × 5 grid of vertical holes on 45 m centres, composited
every 4 m, defined in [`drillholes.jl`](drillholes.jl). One hole passes
within 4 m of the origin, which is where we estimate.

```julia
using Neighborhoods
using Meshes: Point

samples = drillholes()
block = Point(0.0, 0.0, 0.0)

println("composites: ", length(getproperty(values(samples), :Au)))
println("holes:      ", length(unique(getproperty(values(samples), :BHID))))
```

```
composites: 975
holes:      25
```

## Nearest twelve

The ellipsoid is 95 × 95 × 62 m, so it reaches four holes in every
direction. Nothing about that shape is wrong — but with no further
constraint, the twelve nearest composites all come from one hole.

```julia
ball = (95, 95, 62)

plain = NeighborhoodSearch(samples, SearchNeighborhood(ball, maxsamples=12))
flat = searchreport(block, plain)

println("accepted: ", length(flat.indices))
println("by hole:  ", sort(collect(tally(samples, flat.indices, :BHID))))
```

```
accepted: 12
by hole:  ["DH13" => 12]
```

## The same ellipsoid, with sectors and a per-hole cap

`Octants()` splits the neighbourhood by the sign of each rotated
coordinate, and the two quotas below say "at most two samples from any one
octant, and at most three from any one hole". The search budget is
unchanged at twelve samples.

```julia
sectored = NeighborhoodSearch(
  samples,
  SearchNeighborhood(
    ball,
    sectors=Octants(),
    maxsamples=12,
    maxpersector=2,
    category=CategoryRule(:BHID, maxper=3)
  )
)
spread = searchreport(block, sectored)

println("accepted: ", length(spread.indices))
println("by hole:  ", sort(collect(tally(samples, spread.indices, :BHID))))
```

```
accepted: 12
by hole:  ["DH07" => 2, "DH08" => 3, "DH12" => 3, "DH13" => 3, "DH14" => 1]
```

## Side by side

Same data, same ellipsoid, same number of samples — a completely different
neighbourhood.

```julia
printtable(
  ["search", "samples", "holes", "most from one hole", "octants filled"],
  [
    ["nearest 12", length(flat.indices), length(tally(samples, flat.indices, :BHID)),
      topshare(samples, flat.indices, :BHID), flat.nsectors],
    ["sectors + hole cap", length(spread.indices), length(tally(samples, spread.indices, :BHID)),
      topshare(samples, spread.indices, :BHID), spread.nsectors]
  ]
)
```

```
search              samples  holes  most from one hole  octants filled
------------------  -------  -----  ------------------  --------------
nearest 12          12       1      12                  1
sectors + hole cap  12       5      3                   8
```

## It changes the estimate, not just the sample list

Both neighbourhoods are handed to the same ordinary kriging model. The
declustered neighbourhood pulls the estimate away from the single hole that
dominated it.

```julia
using GeoStatsModels: Kriging
using GeoStatsFunctions: SphericalVariogram
using Meshes: PointSet
using Unitful: @u_str

model = Kriging(SphericalVariogram(range=120.0u"m"))
target = PointSet([block])

estimate(search) = first(getproperty(values(
  interpolate(samples, target, model; search, vars=(:Au,))
), :Au))

flatest = estimate(SearchNeighborhood(ball, maxsamples=12))
spreadest = estimate(SearchNeighborhood(
  ball, sectors=Octants(), maxsamples=12, maxpersector=2,
  category=CategoryRule(:BHID, maxper=3)
))

printtable(
  ["search", "Au estimate"],
  [["nearest 12", fmt(flatest)], ["sectors + hole cap", fmt(spreadest)]]
)
println()
println("difference: ", fmt(abs(spreadest - flatest)), " (",
  fmt(100 * abs(spreadest - flatest) / flatest; digits=1), "%)")
```

```
search              Au estimate
------------------  -----------
nearest 12          1.516
sectors + hole cap  1.501

difference: 0.015 (1.0%)
```
