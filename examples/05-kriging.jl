#md # Working with the existing kriging
#md
#md This package does not implement estimation. It decides which samples reach
#md the model and hands them to GeoStatsModels unchanged, so any model from that
#md package works, and the numbers agree with the ecosystem's own `fitpredict`
#md whenever the neighborhood does nothing extra.

using GeoNeighborhoods
using GeoStatsModels: GeoStatsModels, Kriging
using GeoStatsFunctions: SphericalVariogram
using GeoTables: georef, domain
using Meshes: CartesianGrid, nelements
using Unitful: @u_str

samples = drillholes()
grid = CartesianGrid((-75.0, -75.0, -40.0), (75.0, 75.0, 40.0), dims=(6, 6, 4))
model = Kriging(SphericalVariogram(range=120.0u"m"))

#-
#md ## The same answer as `fitpredict`
#md
#md With no sectors and no quotas, the search reduces to "every sample inside
#md the ellipsoid, nearest first", which is exactly what `KBallSearch` supplies
#md to `fitpredict`. The ball below is small enough that everything inside it
#md fits within the sample budget, so neither search has to break a tie, and the
#md two must agree.

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

#-
#md ## What the neighborhood is for
#md
#md Once the neighborhood does something, the estimates diverge, which is the
#md point. The model, the variogram, and the data are identical in both rows
#md below.

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

#-
#md ## An audit trail per block
#md
#md `diagnostics=true` appends the search outcome to every location: which pass
#md served it, how many samples and sectors it used, how many distinct values of
#md the first category column it drew on, and why the search stopped.

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

#-
#md ## Refusing to interpolate
#md
#md A neighborhood that cannot be satisfied leaves the block `missing` rather
#md than estimating it from whatever happened to be nearby. The reason is
#md recorded, so the gaps in a model can be explained.

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
