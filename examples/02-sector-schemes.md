# Sector schemes

Sectors are assigned in the rotated, radius-normalized frame — the frame in
which the ellipsoid is a unit sphere. Anisotropy and sector geometry
therefore cannot disagree, and rotating the search rotates the sectors with
it.

Three schemes are available, and the choice between them carries more
consequence than it appears to.

```julia
using GeoNeighborhoods
using Meshes: Point

samples = drillholes()
block = Point(0.0, 0.0, 0.0)
```

## The schemes compared

The comparison below holds the ellipsoid fixed, along with a budget of 24
samples capped at two per sector. Only the way the neighborhood is
subdivided changes.

```julia
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
```

```
scheme                     sectors  samples  holes  sectors filled
-------------------------  -------  -------  -----  --------------
NoSectors()                1        2        1      1
Octants()                  8        16       4      8
Azimuthal(8)               8        14       7      7
Azimuthal(8, halves=true)  16       24       6      12
```

`NoSectors()` defines a single sector, so `maxpersector=2` throttled the
entire search to two samples. Per-sector quotas are meaningful only
alongside a scheme that actually subdivides the neighborhood.

## Azimuthal sectors measure azimuth about the vertical

This behavior is the trap worth knowing. `Azimuthal(n)` divides the
horizontal plane into wedges, so every composite in a vertical hole shares
one wedge, and capping a sector caps the whole hole. Narrowing the search
onto a single hole makes the effect unmistakable.

```julia
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
```

```
scheme                     samples  above the block  below
-------------------------  -------  ---------------  -----
Azimuthal(8)               2        1                1
Azimuthal(8, halves=true)  4        2                2
Octants()                  4        2                2
```

`halves=true` splits each wedge by the sign of the third axis, and
`Octants` splits by sign in all three. Either scheme recovers the vertical
resolution that plain `Azimuthal` cannot express.

## Rotation carries the sectors with it

A sample sitting on the rotated major axis lands in the same sector
whatever the rotation, because the sector is measured after rotating into
the ellipsoid frame.

```julia
using Rotations: RotZ
using Unitful: @u_str
using StaticArrays: SVector
using LinearAlgebra: norm

for angle in (0, 30, 60)
  spec = SearchNeighborhood((100, 40, 20), rotation=RotZ(deg2rad(angle)), sectors=Octants())
  P = GeoNeighborhoods.localprojection(spec, u"m")
  along = RotZ(deg2rad(angle)) * SVector(60.0, 0.0, 0.0)
  println("rotation ", lpad(angle, 2), "°: octant ", sectorid(Octants(), P * along),
    ", normalized distance ", fmt(norm(P * along)))
end
```

```
rotation  0°: octant 1, normalized distance 0.6
rotation 30°: octant 1, normalized distance 0.6
rotation 60°: octant 1, normalized distance 0.6
```
