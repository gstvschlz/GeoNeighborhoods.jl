@testset "sectors" begin
  @testset "NoSectors" begin
    @test nsectors(NoSectors(), 2) == 1
    @test nsectors(NoSectors(), 3) == 1
    @test sectorid(NoSectors(), SVector(1.0, 2.0, 3.0)) == 1
    @test sectorid(NoSectors(), SVector(-1.0, -2.0, -3.0)) == 1
  end

  @testset "Octants" begin
    @test nsectors(Octants(), 2) == 4
    @test nsectors(Octants(), 3) == 8

    # every sign combination lands in its own sector, and they cover 1:8
    signs = [SVector(sx, sy, sz) for sx in (1.0, -1.0), sy in (1.0, -1.0), sz in (1.0, -1.0)]
    ids = [sectorid(Octants(), u) for u in signs]
    @test sort(vec(ids)) == 1:8

    @test sectorid(Octants(), SVector(1.0, 1.0, 1.0)) == 1
    @test sectorid(Octants(), SVector(-1.0, 1.0, 1.0)) == 2
    @test sectorid(Octants(), SVector(1.0, -1.0, 1.0)) == 3
    @test sectorid(Octants(), SVector(-1.0, -1.0, -1.0)) == 8

    # the centre has no direction
    @test sectorid(Octants(), SVector(0.0, 0.0, 0.0)) == 1

    # 2D quadrants
    @test sort([sectorid(Octants(), SVector(sx, sy)) for sx in (1.0, -1.0), sy in (1.0, -1.0)] |> vec) == 1:4
  end

  @testset "Azimuthal" begin
    @test_throws ArgumentError Azimuthal(1)
    @test nsectors(Azimuthal(12), 2) == 12
    @test nsectors(Azimuthal(12), 3) == 12
    @test nsectors(Azimuthal(6, halves=true), 3) == 12
    # no third axis to split on in 2D
    @test nsectors(Azimuthal(6, halves=true), 2) == 6

    # four wedges: one per axis direction, counter-clockwise from +x
    a = Azimuthal(4)
    @test sectorid(a, SVector(1.0, 0.0)) == 1
    @test sectorid(a, SVector(0.0, 1.0)) == 2
    @test sectorid(a, SVector(-1.0, 0.0)) == 3
    @test sectorid(a, SVector(0.0, -1.0)) == 4

    # ids stay in range all the way around, including the wrap point
    for φ in range(0, 2π, length=257)
      u = SVector(cos(φ), sin(φ))
      @test sectorid(a, u) in 1:4
    end

    # halves put the lower hemisphere in the second block of ids
    h = Azimuthal(4, halves=true)
    @test sectorid(h, SVector(1.0, 0.0, 1.0)) == 1
    @test sectorid(h, SVector(1.0, 0.0, -1.0)) == 5
    @test sectorid(h, SVector(0.0, 1.0, -1.0)) == 6

    # without halves the third coordinate is irrelevant
    @test sectorid(a, SVector(1.0, 0.0, 5.0)) == sectorid(a, SVector(1.0, 0.0, -5.0))
  end

  @testset "consecutiveempty" begin
    # cyclic: a run may wrap around the end of the range
    a = Azimuthal(8)
    @test consecutiveempty(a, 3, [1, 1, 0, 0, 0, 1, 1, 1]) == 3
    @test consecutiveempty(a, 3, [0, 1, 1, 1, 1, 1, 0, 0]) == 3   # wraps 7,8,1
    @test consecutiveempty(a, 3, [1, 1, 1, 1, 1, 1, 1, 1]) == 0
    @test consecutiveempty(a, 3, zeros(Int, 8)) == 8

    # each half wraps independently, never across the equator
    h = Azimuthal(4, halves=true)
    @test consecutiveempty(h, 3, [1, 1, 0, 0, 0, 0, 1, 1]) == 2

    # non-cyclic schemes degrade to a total count
    @test consecutiveempty(Octants(), 3, [1, 0, 1, 0, 1, 0, 1, 0]) == 4
    @test consecutiveempty(Octants(), 3, [1, 1, 1, 1, 1, 1, 1, 1]) == 0
    @test consecutiveempty(NoSectors(), 3, [0]) == 1
  end

  @testset "sector schemes print readably" begin
    @test sprint(show, NoSectors()) == "NoSectors()"
    @test sprint(show, Octants()) == "Octants()"
    @test sprint(show, Azimuthal(8)) == "Azimuthal(8)"
    @test sprint(show, Azimuthal(8, halves=true)) == "Azimuthal(8, halves=true)"
  end
end
