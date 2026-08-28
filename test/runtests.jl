using Neighborhoods
using StaticArrays
using LinearAlgebra: norm
using Unitful: @u_str
using Rotations: RotZ
using Distances: evaluate
using Meshes: metric
using Test

# internals under test that are deliberately not exported
using Neighborhoods: sectorgroups, consecutiveempty, side, cap, satisfied, localprojection

@testset "Neighborhoods.jl" begin
  include("sectors.jl")
  include("constraints.jl")
  include("spec.jl")
end
