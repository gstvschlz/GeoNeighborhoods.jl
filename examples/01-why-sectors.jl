#md # Why sector search matters
#md
#md A block estimate is only as good as the samples that inform it. Down-hole
#md sampling is typically an order of magnitude denser than the spacing between
#md drillholes, so a plain "nearest N samples" search collapses onto whichever
#md hole passes closest to the block. What reaches the estimator is then a
#md vertical string of correlated composites standing in for a
#md three-dimensional neighborhood.
#md
#md The data below is synthetic: a 5 × 5 grid of vertical holes on 45 m
#md centers, composited every 4 m, defined in [`drillholes.jl`](drillholes.jl).
#md One hole passes within 4 m of the origin, and the origin is where every
#md estimate in this document is made.

using GeoNeighborhoods
using Meshes: Point

samples = drillholes()
block = Point(0.0, 0.0, 0.0)

println("composites: ", length(getproperty(values(samples), :Au)))
println("holes:      ", length(unique(getproperty(values(samples), :BHID))))

#-
#md ## Nearest twelve
#md
#md The ellipsoid measures 95 × 95 × 62 m, so it reaches four holes in every
#md direction. Nothing about that shape is wrong. Under no further constraint,
#md however, the twelve nearest composites all come from one hole.

ball = (95, 95, 62)

plain = NeighborhoodSearch(samples, SearchNeighborhood(ball, maxsamples=12))
flat = searchreport(block, plain)

println("accepted: ", length(flat.indices))
println("by hole:  ", sort(collect(tally(samples, flat.indices, :BHID))))

#-
#md ## The same ellipsoid, with sectors and a per-hole cap
#md
#md `Octants()` splits the neighborhood by the sign of each rotated coordinate.
#md The two quotas below admit at most two samples from any one octant and at
#md most three from any one hole. The search budget stays at twelve samples.

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
#md ## The two neighborhoods side by side
#md
#md The table below reports the composition of each selection. Neither the
#md data, the ellipsoid, nor the number of samples differs between the two
#md rows — only the constraints do, and the neighborhoods they produce have
#md nothing in common.

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
#md ## The neighborhood changes the estimate, not merely the sample list
#md
#md Both neighborhoods are handed to the same ordinary kriging model, with the
#md same variogram and the same data. The declustered neighborhood pulls the
#md estimate away from the one hole that had dominated it.

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
