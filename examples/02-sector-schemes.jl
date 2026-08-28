#md # Sector schemes
#md
#md Sectors are assigned in the **rotated, radius-normalised frame** — the frame
#md in which the ellipsoid is a unit sphere. Anisotropy and sector geometry
#md therefore cannot disagree: rotating the search rotates the sectors with it.
#md
#md Three schemes are available, and the choice matters more than it looks.

using GeoNeighborhoods
using Meshes: Point

samples = drillholes()
block = Point(0.0, 0.0, 0.0)

#-
#md ## The schemes compared
#md
#md The same ellipsoid and the same budget of 24 samples, capped at two per
#md sector. What differs is only how the neighbourhood is carved up.

ball = (95, 95, 62)

function compare(scheme)
  s = SearchNeighborhood(ball, sectors=scheme, maxsamples=24, maxpersector=2)
  r = searchreport(block, NeighborhoodSearch(samples, s))
  [
    string(scheme),
    nsectors(s),
    length(r.indices),
    length(tally(samples, r.indices, :BHID)),
    r.nsectors
  ]
end

printtable(
  ["scheme", "sectors", "samples", "holes", "sectors filled"],
  [compare(s) for s in (NoSectors(), Octants(), Azimuthal(8), Azimuthal(8, halves=true))]
)

#-
#md With `NoSectors()` there is only one sector, so `maxpersector=2` throttles
#md the entire search to two samples. Per-sector quotas are only meaningful
#md alongside a scheme that actually subdivides the neighbourhood.
#md
#md ## Azimuthal measures azimuth about the vertical
#md
#md This is the trap worth knowing. `Azimuthal(n)` divides the horizontal plane
#md into wedges, so **every composite in a vertical hole shares one wedge**. Cap
#md a sector and you cap the whole hole. Narrow the search onto a single hole
#md and the effect is unmistakable:

tall = (30, 30, 60)   # reaches only the hole through the block

function onehole(scheme)
  s = SearchNeighborhood(tall, sectors=scheme, maxsamples=24, maxpersector=2)
  r = searchreport(block, NeighborhoodSearch(samples, s))
  zs = [_z(samples, i) for i in r.indices]
  [string(scheme), length(r.indices), count(≥(0), zs), count(<(0), zs)]
end

printtable(
  ["scheme", "samples", "above the block", "below"],
  [onehole(s) for s in (Azimuthal(8), Azimuthal(8, halves=true), Octants())]
)

#-
#md `halves=true` splits each wedge by the sign of the third axis, and `Octants`
#md splits by sign in all three. Either recovers the vertical resolution that
#md plain `Azimuthal` cannot express.
#md
#md ## Rotation carries the sectors with it
#md
#md A sample sitting on the rotated major axis lands in the same sector
#md regardless of the rotation, because the sector is measured after rotating
#md into the ellipsoid frame.

using Rotations: RotZ
using Unitful: @u_str
using StaticArrays: SVector
using LinearAlgebra: norm

for angle in (0, 30, 60)
  spec = SearchNeighborhood((100, 40, 20), rotation=RotZ(deg2rad(angle)), sectors=Octants())
  P = GeoNeighborhoods.localprojection(spec, u"m")
  along = RotZ(deg2rad(angle)) * SVector(60.0, 0.0, 0.0)
  println("rotation ", lpad(angle, 2), "°: octant ", sectorid(Octants(), P * along),
    ", normalised distance ", fmt(norm(P * along)))
end
