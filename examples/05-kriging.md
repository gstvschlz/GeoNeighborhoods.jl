# Working with the existing kriging

This package does not implement estimation. It decides which samples reach
the model and hands them to GeoStatsModels unchanged, so any model from that
package works, and the numbers agree with the ecosystem's own `fitpredict`
whenever the neighborhood does nothing extra.

```julia
using GeoNeighborhoods
using GeoStatsModels: GeoStatsModels, Kriging
using GeoStatsFunctions: SphericalVariogram
using GeoTables: georef, domain
using Meshes: CartesianGrid, nelements
using Unitful: @u_str

samples = drillholes()
grid = CartesianGrid((-75.0, -75.0, -40.0), (75.0, 75.0, 40.0), dims=(6, 6, 4))
model = Kriging(SphericalVariogram(range=120.0u"m"))
```

## The same answer as `fitpredict`

With no sectors and no quotas, the search reduces to "every sample inside
the ellipsoid, nearest first", which is exactly what `KBallSearch` supplies
to `fitpredict`. The ball below is small enough that everything inside it
fits within the sample budget, so neither search has to break a tie, and the
two must agree.

```julia
spec = SearchNeighborhood((30, 30, 30), maxsamples=60)

ours = interpolate(samples, grid, model; search=spec, vars=(:Au,))
theirs = GeoStatsModels.fitpredict(
  model,
  georef((Au=getproperty(values(samples), :Au),), domain(samples)),
  grid;
  neighbors=true, minneighbors=1, maxneighbors=60, neighborhood=metricball(spec)
)

a = getproperty(values(ours), :Au)
b = getproperty(values(theirs), :Au)
both = .!ismissing.(a) .& .!ismissing.(b)

println("blocks:              ", length(a))
println("estimated by both:   ", count(both))
println("same missing blocks: ", ismissing.(a) == ismissing.(b))
println("largest difference:  ", maximum(abs.(a[both] .- b[both])))
```

```
blocks:              144
estimated by both:   144
same missing blocks: true
largest difference:  8.881784197001252e-16
```

## What the neighborhood is for

Once the neighborhood does something, the estimates diverge, which is the
point. The model, the variogram, and the data are identical in both rows
below.

```julia
function estimates(search)
  au = getproperty(values(interpolate(samples, grid, model; search, vars=(:Au,))), :Au)
  au, count(!ismissing, au)
end

wide = SearchNeighborhood((95, 95, 62), maxsamples=12)
declustered = SearchNeighborhood(
  (95, 95, 62), sectors=Octants(), maxsamples=12, maxpersector=2,
  category=CategoryRule(:BHID, maxper=3)
)

x, nx = estimates(wide)
y, ny = estimates(declustered)
shared = .!ismissing.(x) .& .!ismissing.(y)

printtable(
  ["search", "blocks estimated", "mean Au"],
  [
    ["nearest 12", nx, fmt(sum(skipmissing(x)) / nx)],
    ["octants + hole cap", ny, fmt(sum(skipmissing(y)) / ny)]
  ]
)
println()
println("largest difference on shared blocks: ", fmt(maximum(abs.(x[shared] .- y[shared]))))
```

```
search              blocks estimated  mean Au
------------------  ----------------  -------
nearest 12          144               1.494
octants + hole cap  144               1.507

largest difference on shared blocks: 0.068
```

## An audit trail per block

`diagnostics=true` appends the search outcome to every location: which pass
served it, how many samples and sectors it used, how many distinct values of
the first category column it drew on, and why the search stopped.

```julia
est = interpolate(samples, grid, model; search=declustered, vars=(:Au,), diagnostics=true)
tab = values(est)

using Tables
println("columns: ", Tables.columnnames(Tables.columns(tab)))
println()

printtable(
  ["block", "Au", "samples", "sectors", "holes", "reason"],
  [[i, fmt(getproperty(tab, :Au)[i]), getproperty(tab, :nsamples)[i],
    getproperty(tab, :nsectors)[i], getproperty(tab, :ndistinct)[i],
    string(getproperty(tab, :reject)[i])] for i in 1:6]
)
```

```
columns: (:Au, :pass, :nsamples, :nsectors, :ndistinct, :reject)

block  Au     samples  sectors  holes  reason
-----  -----  -------  -------  -----  --------
1      1.438  12       8        4      Accepted
2      1.539  12       8        4      Accepted
3      1.64   12       8        4      Accepted
4      1.736  12       8        4      Accepted
5      1.84   12       8        4      Accepted
6      1.936  12       8        4      Accepted
```

## Refusing to interpolate

A neighborhood that cannot be satisfied leaves the block `missing` rather
than estimating it from whatever happened to be nearby. The reason is
recorded, so the gaps in a model can be explained.

```julia
strict = SearchNeighborhood(
  (40, 40, 40), sectors=Octants(), minsamples=8, minsectors=6,
  category=CategoryRule(:BHID, mindistinct=3)
)
est = interpolate(samples, grid, model; search=strict, vars=(:Au,), diagnostics=true)
tab = values(est)
reasons = getproperty(tab, :reject)

printtable(
  ["outcome", "blocks"],
  sort([[string(r), count(==(r), reasons)] for r in unique(reasons)], by=x -> -x[2])
)
println()
println("estimated: ", count(!ismissing, getproperty(tab, :Au)), " of ", nelements(grid))
```

```
outcome        blocks
-------------  ------
TooFewSectors  100
Accepted       44

estimated: 44 of 144
```
