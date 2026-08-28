# Synthetic drillhole data shared by the tests and the examples.
#
# A 5 x 5 grid of vertical holes on 45 m centres, composited every 4 m: the
# ordinary case where down-hole sampling is an order of magnitude denser than
# the spacing between holes.
#
# The grid is deliberately offset by a few metres so that no collar sits
# exactly on x = 0 or y = 0. Collars on the axes would make sign-based octants
# degenerate, since a coordinate of exactly zero always counts as positive.
# One hole passes within 4 m of the origin, so a plain nearest-N search there
# collapses onto that single hole.

"""
    drillholes(; ngrid=5, spacing=45.0, composite=4.0, halfheight=76.0)

A geotable of composites with columns `Au` (grade), `BHID` (hole id) and
`ROCK` (a domain code split by the plane x = 0).
"""
function drillholes(; ngrid=5, spacing=45.0, composite=4.0, halfheight=76.0)
  points = Point[]
  bhid = String[]
  au = Float64[]
  rock = String[]

  half = (ngrid - 1) ÷ 2
  h = 0
  for k in (-half):half, j in (-half):half
    h += 1
    x = spacing * k + 3.0
    y = spacing * j + 2.0
    for z in (-halfheight):composite:halfheight
      push!(points, Point(x, y, z))
      push!(bhid, "DH" * lpad(h, 2, '0'))
      # a smooth trend, so kriging has something sensible to reproduce
      push!(au, 1.5 + 0.004x + 0.002y - 0.010z)
      push!(rock, x < 0 ? "WEST" : "EAST")
    end
  end

  georef((Au=au, BHID=bhid, ROCK=rock), PointSet(points))
end

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
