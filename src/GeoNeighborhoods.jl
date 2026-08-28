"""
    GeoNeighborhoods

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
module GeoNeighborhoods

using Meshes: MetricBall, metric, BoundedNeighborSearchMethod, Point, PointSet, Domain
using Meshes: centroid, nelements, crs, lentype, embeddim, coords
using GeoTables: domain, georef
import Meshes: maxneighbors, searchdists!

using CoordRefSystems: CRS
import CoordRefSystems

using NearestNeighbors: BallTree, inrange
import Tables
import GeoStatsModels

using StaticArrays: SVector, SMatrix
using LinearAlgebra: I, UniformScaling, Diagonal, norm
using Unitful: Unitful, @u_str, ustrip, unit

include("sectors.jl")
include("constraints.jl")
include("spec.jl")
include("searcher.jl")
include("multipass.jl")
include("estimate.jl")

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

  # neighborhood specification
  SearchNeighborhood,
  metricball,

  # search
  NeighborhoodSearch,
  MultiPass,
  MultiPassSearch,
  searchreport,
  searchreport!,

  # estimation
  interpolate

end
