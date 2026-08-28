#md # Why sector search
#md
#md A block estimate is only as good as the samples that inform it. Down-hole
#md sampling is typically an order of magnitude denser than the spacing between
#md holes, so a plain "nearest N samples" search collapses onto whichever hole
#md happens to pass closest — a vertical string of correlated composites
#md standing in for a three-dimensional neighbourhood.
#md
#md The data here is a 5 × 5 grid of vertical holes on 45 m centres, composited
#md every 4 m, defined in [`drillholes.jl`](drillholes.jl). One hole passes
#md within 4 m of the origin, which is where we estimate.

using GeoNeighborhoods
using Meshes: Point

samples = drillholes()
block = Point(0.0, 0.0, 0.0)

println("composites: ", length(getproperty(values(samples), :Au)))
println("holes:      ", length(unique(getproperty(values(samples), :BHID))))

#-
#md ## Nearest twelve
#md
#md The ellipsoid is 95 × 95 × 62 m, so it reaches four holes in every
#md direction. Nothing about that shape is wrong — but with no further
#md constraint, the twelve nearest composites all come from one hole.

ball = (95, 95, 62)

plain = NeighborhoodSearch(samples, SearchNeighborhood(ball, maxsamples=12))
flat = searchreport(block, plain)

println("accepted: ", length(flat.indices))
println("by hole:  ", sort(collect(tally(samples, flat.indices, :BHID))))

#-
#md ## The same ellipsoid, with sectors and a per-hole cap
#md
#md `Octants()` splits the neighbourhood by the sign of each rotated
#md coordinate, and the two quotas below say "at most two samples from any one
#md octant, and at most three from any one hole". The search budget is
#md unchanged at twelve samples.

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

#-
#md ## Side by side
#md
#md Same data, same ellipsoid, same number of samples — a completely different
#md neighbourhood.

printtable(
  ["search", "samples", "holes", "most from one hole", "octants filled"],
  [
    ["nearest 12", length(flat.indices), length(tally(samples, flat.indices, :BHID)),
      topshare(samples, flat.indices, :BHID), flat.nsectors],
    ["sectors + hole cap", length(spread.indices), length(tally(samples, spread.indices, :BHID)),
      topshare(samples, spread.indices, :BHID), spread.nsectors]
  ]
)

#-
#md ## It changes the estimate, not just the sample list
#md
#md Both neighbourhoods are handed to the same ordinary kriging model. The
#md declustered neighbourhood pulls the estimate away from the single hole that
#md dominated it.

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
