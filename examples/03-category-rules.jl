#md # Category rules
#md
#md A `CategoryRule` ties sample selection to a categorical column — a drillhole
#md id, a domain code, a weathering zone. Four kinds of constraint are
#md available, and they compose.

using GeoNeighborhoods
using Meshes: Point

samples = drillholes()
block = Point(0.0, 0.0, 0.0)
ball = (95, 95, 62)

search(rule) = NeighborhoodSearch(samples, SearchNeighborhood(ball, maxsamples=12, category=rule))

#-
#md ## A cap per value
#md
#md The classic rule: no more than N composites from the same hole. Without it
#md the twelve nearest samples come from one hole; with it they cannot.

for cap in (nothing, 4, 2)
  r = if isnothing(cap)
    searchreport(block, NeighborhoodSearch(samples, SearchNeighborhood(ball, maxsamples=12)))
  else
    searchreport(block, search(CategoryRule(:BHID, maxper=cap)))
  end
  println(rpad(isnothing(cap) ? "no cap" : "maxper=$cap", 10),
    " samples=", lpad(length(r.indices), 2),
    "  holes=", lpad(length(tally(samples, r.indices, :BHID)), 2),
    "  most from one hole=", topshare(samples, r.indices, :BHID))
end

#-
#md Note that `maxper` alone does not spread the search across space — it
#md spreads it across holes. With a cap of 2 the search still takes the nearest
#md twelve samples it is allowed, which is six different holes.
#md
#md ## A minimum number of distinct values
#md
#md `mindistinct` is a rejection rule: if the accepted samples come from too few
#md holes, the location is not interpolated at all.
#md
#md On its own it is close to useless, and the reason is worth seeing. The
#md greedy fill takes the nearest samples first, so without a cap it empties the
#md closest hole before reaching a second one — the rule then rejects a location
#md that had plenty of holes within range:

for radius in (12, 30, 95)
  r = searchreport(block, NeighborhoodSearch(
    samples, SearchNeighborhood((radius, radius, radius), category=CategoryRule(:BHID, mindistinct=3))
  ))
  println("radius ", lpad(radius, 2), " m: ", rpad(string(r.reason), 14),
    " samples=", lpad(length(r.indices), 2), "  holes=", r.ndistinct)
end

#-
#md Pair it with a cap and the two rules do what was intended: spread the
#md selection, then insist the spread was actually achieved.

for radius in (12, 30, 95)
  r = searchreport(block, NeighborhoodSearch(
    samples,
    SearchNeighborhood((radius, radius, radius), category=CategoryRule(:BHID, maxper=3, mindistinct=3))
  ))
  println("radius ", lpad(radius, 2), " m: ", rpad(string(r.reason), 14),
    " samples=", lpad(length(r.indices), 2), "  holes=", r.ndistinct)
end

#-
#md ## Hard domain matching
#md
#md `match=:block` accepts only samples whose value equals the estimated
#md location's own. The value has to come from somewhere, so the target must
#md carry the column — [`interpolate`](../src/estimate.jl) reads it from the
#md target geotable, and a direct search takes it as `blockvals`.

hard = NeighborhoodSearch(samples, SearchNeighborhood(ball, maxsamples=12, category=CategoryRule(:ROCK, match=:block)))

for rock in ("WEST", "EAST")
  r = searchreport(block, hard; blockvals=(; ROCK=rock))
  println("block in ", rpad(rock, 5), ": ", collect(keys(tally(samples, r.indices, :ROCK))),
    "  samples=", length(r.indices))
end

#-
#md Without `blockvals` the search cannot answer the question, and says so
#md rather than quietly ignoring the rule:

try
  searchreport(block, hard)
catch e
  println(sprint(showerror, e))
end

#-
#md ## Explicit per-value quotas
#md
#md `quotas` sets a floor and a ceiling per value. The ceiling is enforced while
#md selecting; the floor is checked afterwards and rejects the location when it
#md cannot be met. Where `maxper` and a quota ceiling both apply, the tighter
#md one wins.

quota = NeighborhoodSearch(samples, SearchNeighborhood(
  ball, maxsamples=12,
  category=CategoryRule(:ROCK, quotas=Dict("WEST" => 4:6, "EAST" => 0:6))
))
r = searchreport(block, quota)
println("with quotas WEST 4:6, EAST 0:6")
println("  reason:  ", r.reason)
println("  by rock: ", sort(collect(tally(samples, r.indices, :ROCK))))

impossible = NeighborhoodSearch(samples, SearchNeighborhood(
  (30, 30, 30), category=CategoryRule(:ROCK, quotas=Dict("WEST" => 5:9))
))
r = searchreport(block, impossible)
println("\nwith a floor the data cannot meet nearby")
println("  reason:  ", r.reason, "  samples=", length(r.indices))
