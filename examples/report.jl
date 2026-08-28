using Printf
using Tables
using CoordRefSystems: CoordRefSystems
using GeoTables: domain
using Meshes: centroid, coords

"counts of accepted samples per value of `column`"
function tally(data, indices, column)
  vals = Tables.getcolumn(Tables.columns(values(data)), column)
  counts = Dict{eltype(vals),Int}()
  for i in indices
    counts[vals[i]] = get(counts, vals[i], 0) + 1
  end
  counts
end

"how many samples came from the most heavily used value of `column`"
topshare(data, indices, column) = isempty(indices) ? 0 : maximum(values(tally(data, indices, column)))

function printtable(headers, rows)
  cells = Vector{String}[[string(h) for h in headers]]
  for r in rows
    push!(cells, [string(c) for c in r])
  end
  widths = [maximum(length(row[i]) for row in cells) for i in eachindex(headers)]
  line(row) = rstrip(join((rpad(row[i], widths[i]) for i in eachindex(widths)), "  "))
  println(line(cells[1]))
  println(join(("-"^w for w in widths), "  "))
  for row in cells[2:end]
    println(line(row))
  end
end

fmt(x::Real; digits=3) = string(round(x; digits))
fmt(::Missing; digits=3) = "missing"

"z coordinate of element `i`, for describing a selection above/below a block"
_z(data, i) = CoordRefSystems.raw(coords(centroid(domain(data), i)))[3]
