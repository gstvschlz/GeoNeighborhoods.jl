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
    The design is settled (see `docs/superpowers/specs/`) but the
    implementation has not landed yet. The API below is the target.
"""
module Neighborhoods

# Design settled 2026-08-28. Implementation follows the plan in
# docs/superpowers/specs/2026-08-28-neighborhoods-jl-design.md
#
#   sectors.jl      Octants / Azimuthal sector assignment
#   constraints.jl  category quotas, half-space split, rejection rules
#   spec.jl         SearchNeighborhood
#   searcher.jl     NeighborhoodSearch <: BoundedNeighborSearchMethod
#   multipass.jl    MultiPass driver
#   estimate.jl     prediction loop + TableTransform

end
