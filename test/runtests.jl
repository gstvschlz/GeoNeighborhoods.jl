using Neighborhoods
using StaticArrays
using LinearAlgebra: norm
using Unitful: @u_str, ustrip
using Rotations: RotZ
using Distances: evaluate
using Meshes: metric, maxneighbors, Point, PointSet, KBallSearch, search, searchdists, centroid, coords
using GeoTables: georef, domain
using CoordRefSystems: CoordRefSystems
import Tables
using Test

# internals under test that are deliberately not exported
using Neighborhoods: sectorgroups, consecutiveempty, side, cap, satisfied, localprojection

include("fixtures.jl")

@testset "Neighborhoods.jl" begin
  include("sectors.jl")
  include("constraints.jl")
  include("spec.jl")
  include("searcher.jl")
  include("multipass.jl")
end
