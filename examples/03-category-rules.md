# Category rules

A `CategoryRule` ties sample selection to a categorical column, such as a
drillhole identifier, a domain code, or a weathering zone. Four kinds of
constraint are available, and they compose.

```julia
using GeoNeighborhoods
using Meshes: Point

samples = drillholes()
block = Point(0.0, 0.0, 0.0)
ball = (95, 95, 62)

search(rule) = NeighborhoodSearch(samples, SearchNeighborhood(ball, maxsamples=12, category=rule))
```

## A cap per value

The classic rule admits no more than N composites from the same hole.
Without it the twelve nearest samples come from a single hole; with it they
cannot.

```julia
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
```

```
no cap     samples=12  holes= 1  most from one hole=12
maxper=4   samples=12  holes= 3  most from one hole=4
maxper=2   samples=12  holes= 6  most from one hole=2
```

A cap does not spread the search across space. It spreads the search across
holes: under `maxper=2` the selection remains the nearest twelve samples
the rule allows, drawn from six different holes.

## A minimum number of distinct values

`mindistinct` is a rejection rule. When the accepted samples come from too
few holes, the location is not interpolated at all.

On its own the rule is close to useless, and the reason repays attention.
The greedy fill takes the nearest samples first, so without a cap it empties
the closest hole before reaching a second one. The rule then rejects a
location that had ample holes within range.

```julia
for radius in (12, 30, 95)
  r = searchreport(block, NeighborhoodSearch(
    samples, SearchNeighborhood((radius, radius, radius), category=CategoryRule(:BHID, mindistinct=3))
  ))
  println("radius ", lpad(radius, 2), " m: ", rpad(string(r.reason), 14),
    " samples=", lpad(length(r.indices), 2), "  holes=", r.ndistinct)
end
```

```
radius 12 m: CategoryUnmet  samples= 0  holes=1
radius 30 m: CategoryUnmet  samples= 0  holes=1
radius 95 m: CategoryUnmet  samples= 0  holes=2
```

Paired with a cap, the two rules do what was intended: they spread the
selection, then insist that the spread was actually achieved.

```julia
for radius in (12, 30, 95)
  r = searchreport(block, NeighborhoodSearch(
    samples,
    SearchNeighborhood((radius, radius, radius), category=CategoryRule(:BHID, maxper=3, mindistinct=3))
  ))
  println("radius ", lpad(radius, 2), " m: ", rpad(string(r.reason), 14),
    " samples=", lpad(length(r.indices), 2), "  holes=", r.ndistinct)
end
```

```
radius 12 m: CategoryUnmet  samples= 0  holes=1
radius 30 m: CategoryUnmet  samples= 0  holes=1
radius 95 m: Accepted       samples=24  holes=8
```

## Hard domain matching

`match=:block` accepts only samples whose value equals the estimated
location's own. That value has to come from somewhere, so the target must
carry the column: [`interpolate`](../src/estimate.jl) reads it from the
target geotable, and a direct search takes it as `blockvals`.

```julia
hard = NeighborhoodSearch(samples, SearchNeighborhood(ball, maxsamples=12, category=CategoryRule(:ROCK, match=:block)))

for rock in ("WEST", "EAST")
  r = searchreport(block, hard; blockvals=(; ROCK=rock))
  println("block in ", rpad(rock, 5), ": ", collect(keys(tally(samples, r.indices, :ROCK))),
    "  samples=", length(r.indices))
end
```

```
block in WEST : ["WEST"]  samples=12
block in EAST : ["EAST"]  samples=12
```

Without `blockvals` the search cannot answer the question, and it says so
rather than quietly ignoring the rule.

```julia
try
  searchreport(block, hard)
catch e
  println(sprint(showerror, e))
end
```

```
ArgumentError: CategoryRule on :ROCK uses match=:block, which needs the estimated location's own value; pass blockvals=(; ROCK=…), or use interpolate with a geotable target
```

## Explicit per-value quotas

`quotas` sets a floor and a ceiling per value. The ceiling is enforced while
selecting, whereas the floor is checked afterwards and rejects the location
when it cannot be met. Where `maxper` and a quota ceiling both apply, the
tighter of the two wins.

```julia
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
```

```
with quotas WEST 4:6, EAST 0:6
  reason:  Accepted
  by rock: ["EAST" => 6, "WEST" => 6]

with a floor the data cannot meet nearby
  reason:  CategoryUnmet  samples=0
```
