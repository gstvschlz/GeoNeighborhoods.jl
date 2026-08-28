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

using Meshes: MetricBall
using StaticArrays: SVector, SMatrix
using LinearAlgebra: I, UniformScaling, Diagonal
using Unitful: Unitful, @u_str, ustrip

include("sectors.jl")
include("constraints.jl")
include("spec.jl")

export
  # sector schemes
  SectorScheme,
  NoSectors,
  Octants,
  Azimuthal,
  nsectors,
  sectorid,

  # constraints
  HalfSpace,
  CategoryRule,
  RejectReason,

  # neighbourhood specification
  SearchNeighborhood,
  metricball

end
