using GeoTables: georef
using Meshes: Point, PointSet

"""
    drillholes(; ngrid=5, spacing=45.0, composite=4.0, halfheight=76.0)

Synthetic composites from a `ngrid × ngrid` grid of vertical holes on `spacing`
centers, sampled every `composite` meters. Down-hole sampling is an order of
magnitude denser than the spacing between holes, which is the ordinary case and
the reason a plain nearest-N search degenerates.

The grid is offset a few meters off the axes so that no collar sits exactly on
`x = 0` or `y = 0`: a coordinate of exactly zero always counts as positive, so
collars on an axis would make sign-based octants degenerate. One hole passes
within 4 m of the origin.

Columns are `Au` (grade), `BHID` (hole id) and `ROCK` (a domain code split by
the plane `x = 0`).
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
      push!(au, 1.5 + 0.004x + 0.002y - 0.010z)
      push!(rock, x < 0 ? "WEST" : "EAST")
    end
  end

  georef((Au=au, BHID=bhid, ROCK=rock), PointSet(points))
end
