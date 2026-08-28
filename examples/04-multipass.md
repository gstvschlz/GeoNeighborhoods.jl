# Multipass searches

A single neighborhood forces one compromise on the whole model: tight enough
to be defensible near data, and every block further out goes uninterpolated.
A `MultiPass` states the compromise explicitly instead. It holds an ordered
list of complete neighborhoods, each with its own quotas and rejection
rules, and the first one whose constraints hold serves the block. The output
records which pass did so.

```julia
using GeoNeighborhoods
using Meshes: CartesianGrid, nelements

samples = drillholes()
grid = CartesianGrid((-260.0, -260.0, -60.0), (260.0, 260.0, 60.0), dims=(13, 13, 4))

println("blocks: ", nelements(grid))
```

```
blocks: 676
```

## One neighborhood at a time

Each neighborhood on its own trades coverage against how well informed the
estimate is. The tight pass demands eight samples from at least four
octants, whereas the wide pass asks for almost nothing.

```julia
tight = SearchNeighborhood((35, 35, 35), sectors=Octants(), minsamples=8, minsectors=4)
medium = SearchNeighborhood((70, 70, 55), sectors=Octants(), minsamples=5, minsectors=2)
wide = SearchNeighborhood((140, 140, 90), minsamples=2)

using GeoStatsModels: Kriging
using GeoStatsFunctions: SphericalVariogram
using Unitful: @u_str

model = Kriging(SphericalVariogram(range=120.0u"m"))

function coverage(search)
  est = interpolate(samples, grid, model; search, vars=(:Au,), diagnostics=true)
  tab = values(est)
  au = getproperty(tab, :Au)
  filled = count(!ismissing, au)
  Any[filled, round(100 * filled / length(au), digits=1), round(sum(getproperty(tab, :nsamples)) / max(filled, 1), digits=1)]
end

printtable(
  ["neighborhood", "blocks estimated", "% of model", "mean samples used"],
  [
    ["tight", coverage(tight)...],
    ["medium", coverage(medium)...],
    ["wide", coverage(wide)...]
  ]
)
```

```
neighborhood  blocks estimated  % of model  mean samples used
------------  ----------------  ----------  -----------------
tight         36                5.3         21.6
medium        232               34.3        21.1
wide          468               69.2        23.4
```

## All three, in order

Every block now receives the tightest neighborhood able to serve it.

```julia
passes = MultiPass(tight, medium, wide)
est = interpolate(samples, grid, model; search=passes, vars=(:Au,), diagnostics=true)
tab = values(est)

au = getproperty(tab, :Au)
pass = getproperty(tab, :pass)
used = .!ismissing.(au)

printtable(
  ["pass", "blocks", "% of model"],
  [Any[p, count(used .& (pass .== p)), round(100 * count(used .& (pass .== p)) / length(au), digits=1)]
   for p in 1:length(passes)]
)
println()
println("estimated:     ", count(used), " of ", length(au))
println("not estimated: ", count(ismissing, au))
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

## Why the remaining blocks stayed missing

With `diagnostics=true` every location carries the reason its search ended,
so an unestimated block is explainable rather than mysterious.

```julia
reasons = getproperty(tab, :reject)
printtable(
  ["reason", "blocks"],
  sort([Any[string(r), count(==(r), reasons)] for r in unique(reasons)], by=x -> -x[2])
)
```

```
reason        blocks
------------  ------
Accepted      468
NoCandidates  208
```

## Relaxing one neighborhood outwards

When the passes differ only in size, `MultiPass(base, factors)` builds them
and, on request, relaxes the sample floor alongside the radii.

```julia
expanded = MultiPass(tight, (1, 2, 4), minsamples=(8, 5, 2))
for (i, p) in enumerate(expanded)
  println("pass ", i, ": ", p)
end
```

```
pass 1: SearchNeighborhood((35.0 m, 35.0 m, 35.0 m), Octants, 8-24 samples, ≥4 sectors)
pass 2: SearchNeighborhood((70.0 m, 70.0 m, 70.0 m), Octants, 5-24 samples, ≥4 sectors)
pass 3: SearchNeighborhood((140.0 m, 140.0 m, 140.0 m), Octants, 2-24 samples, ≥4 sectors)
```
