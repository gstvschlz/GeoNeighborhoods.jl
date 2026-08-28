@testset "multipass" begin
  data = drillholes()

  tight = SearchNeighborhood((20, 20, 20), minsamples=8)
  medium = SearchNeighborhood((40, 40, 40), minsamples=5)
  wide = SearchNeighborhood((90, 90, 90), minsamples=2)

  @testset "construction" begin
    m = MultiPass(tight, medium, wide)
    @test length(m) == 3
    @test m[1] === tight
    @test ndims(m) == 3
    @test collect(m) == [tight, medium, wide]
    @test MultiPass([tight, medium]) isa MultiPass

    @test_throws ArgumentError MultiPass(())
    @test_throws ArgumentError MultiPass((tight, "not a neighbourhood"))
    @test_throws ArgumentError MultiPass(tight, SearchNeighborhood((10, 10)))
  end

  @testset "expanding a base neighbourhood" begin
    m = MultiPass(tight, (1, 2, 4))
    @test length(m) == 3
    @test m[1].radii == tight.radii
    @test m[2].radii == 2 .* tight.radii
    @test m[3].radii == 4 .* tight.radii
    @test all(p -> p.minsamples == tight.minsamples, m)

    r = MultiPass(tight, (1, 2, 3), minsamples=(8, 5, 2))
    @test [p.minsamples for p in r] == [8, 5, 2]
    # everything else rides along unchanged
    @test all(p -> p.maxsamples == tight.maxsamples, r)

    @test_throws ArgumentError MultiPass(tight, ())
    @test_throws ArgumentError MultiPass(tight, (1, 2), minsamples=(8,))
  end

  @testset "the first satisfied pass wins" begin
    m = NeighborhoodSearch(data, MultiPass(tight, medium, wide))
    @test maxneighbors(m) == 24

    # a location on a hole is served by the tightest pass
    r = searchreport(Point(3.0, 2.0, 0.0), m)
    @test r.reason == GeoNeighborhoods.Accepted
    @test r.pass == 1
    @test length(r.indices) ≥ 8

    # midway between four collars, out of reach of pass 1
    r = searchreport(Point(25.5, 24.5, 0.0), m)
    @test r.reason == GeoNeighborhoods.Accepted
    @test r.pass == 2
    @test length(r.indices) ≥ 5

    # beyond every pass: reported against the widest one
    r = searchreport(Point(5000.0, 5000.0, 0.0), m)
    @test isempty(r.indices)
    @test r.reason == GeoNeighborhoods.NoCandidates
    @test r.pass == 3
  end

  @testset "a single pass behaves like the neighbourhood itself" begin
    one = NeighborhoodSearch(data, MultiPass(medium))
    bare = NeighborhoodSearch(data, medium)
    p = Point(3.0, 2.0, 0.0)
    @test searchreport(p, one).indices == searchreport(p, bare).indices
    @test searchreport(p, one).pass == 1
  end

  @testset "passes keep their own rules" begin
    # pass 1 demands eight octants and cannot get them; pass 2 has no such rule
    strict = SearchNeighborhood((30, 30, 30), sectors=Octants(), minsectors=8)
    m = NeighborhoodSearch(data, MultiPass(strict, medium))
    r = searchreport(Point(3.0, 2.0, 0.0), m)
    @test r.reason == GeoNeighborhoods.Accepted
    @test r.pass == 2
  end
end
