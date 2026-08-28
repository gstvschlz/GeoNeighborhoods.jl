# The synthetic drillhole data is defined alongside the examples so that the
# documented numbers and the tested behaviour come from one dataset.
include(joinpath(@__DIR__, "..", "examples", "drillholes.jl"))

"count of accepted samples per value of `column`"
function tally(data, indices, column)
  vals = Tables.getcolumn(Tables.columns(values(data)), column)
  counts = Dict{eltype(vals),Int}()
  for i in indices
    counts[vals[i]] = get(counts, vals[i], 0) + 1
  end
  counts
end

"raw coordinates of element `i`"
_xyz(data, i) = CoordRefSystems.raw(coords(centroid(domain(data), i)))

"count occurrences of each element"
function _counts(xs)
  c = Dict{eltype(xs),Int}()
  for x in xs
    c[x] = get(c, x, 0) + 1
  end
  c
end
