#md # Multipass searches
#md
#md A single neighbourhood forces one compromise for the whole model: tight
#md enough to be defensible near data, and every block further out goes
#md uninterpolated. A `MultiPass` states the compromise explicitly instead —
#md an ordered list of complete neighbourhoods, each with its own quotas and
#md rejection rules. The first one whose constraints hold wins, and the output
#md records which.

using Neighborhoods
using Meshes: CartesianGrid, nelements

samples = drillholes()
grid = CartesianGrid((-260.0, -260.0, -60.0), (260.0, 260.0, 60.0), dims=(13, 13, 4))

println("blocks: ", nelements(grid))

#-
#md ## One neighbourhood at a time
#md
#md Each of these on its own trades coverage against how well informed the
#md estimate is. The tight pass demands eight samples from at least four
#md octants; the wide one asks for almost nothing.

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
  ["neighbourhood", "blocks estimated", "% of model", "mean samples used"],
  [
    ["tight", coverage(tight)...],
    ["medium", coverage(medium)...],
    ["wide", coverage(wide)...]
  ]
)

#-
#md ## All three, in order
#md
#md Every block gets the tightest neighbourhood that can actually serve it.

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

#-
#md ## Why the blocks that stayed missing did
#md
#md With `diagnostics=true` every location carries the reason its search ended,
#md so an unestimated block is explainable rather than mysterious.

reasons = getproperty(tab, :reject)
printtable(
  ["reason", "blocks"],
  sort([Any[string(r), count(==(r), reasons)] for r in unique(reasons)], by=x -> -x[2])
)

#-
#md ## Relaxing one neighbourhood outwards
#md
#md When the passes differ only in size, `MultiPass(base, factors)` builds them,
#md optionally relaxing the sample floor alongside the radii.

expanded = MultiPass(tight, (1, 2, 4), minsamples=(8, 5, 2))
for (i, p) in enumerate(expanded)
  println("pass ", i, ": ", p)
end
