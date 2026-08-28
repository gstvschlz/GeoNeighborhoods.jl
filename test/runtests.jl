using Neighborhoods
using StaticArrays
using Test

# internals under test that are deliberately not exported
using Neighborhoods: sectorgroups, consecutiveempty

@testset "Neighborhoods.jl" begin
  include("sectors.jl")
end
