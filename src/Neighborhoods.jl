# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Neighborhoods

Mining-style search neighborhoods for the GeoStats.jl ecosystem:
anisotropic ellipsoids with angular sectors, per-sector and
per-category sample quotas, half-space balancing, and multipass
searches.

The search is exposed as a `Meshes.BoundedNeighborSearchMethod`,
so it composes with the existing estimation machinery instead of
replacing it.

!!! warning "Pre-release"
    Under active implementation; the public API is still settling.
"""
module Neighborhoods

include("sectors.jl")

export
  # sector schemes
  SectorScheme,
  NoSectors,
  Octants,
  Azimuthal,
  nsectors,
  sectorid

end
