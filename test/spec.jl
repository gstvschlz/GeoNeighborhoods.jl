@testset "spec" begin
  @testset "construction" begin
    s = SearchNeighborhood((100, 50, 20))
    @test ndims(s) == 3
    @test s.radii == (100.0u"m", 50.0u"m", 20.0u"m")   # plain numbers mean metres
    @test s.sectors === NoSectors()
    @test nsectors(s) == 1
    @test s.minsamples == 1
    @test s.maxsamples == 24
    @test s.maxpersector == GeoNeighborhoods.UNBOUNDED
    @test s.maxemptyconsecutive == GeoNeighborhoods.UNBOUNDED
    @test isnothing(s.split)
    @test isempty(s.category)

    # units are honoured, mixed units promote to a common one, and the
    # keyword form agrees with the positional one
    @test all(SearchNeighborhood((1u"km", 500u"m", 0.2u"km")).radii .≈ (1000u"m", 500u"m", 200u"m"))
    @test allequal(typeof.(SearchNeighborhood((1u"km", 500u"m", 0.2u"km")).radii))
    @test SearchNeighborhood(radii=(10, 5)).radii == SearchNeighborhood((10, 5)).radii

    @test ndims(SearchNeighborhood((10, 5))) == 2

    # a single rule need not be wrapped in a collection
    one = SearchNeighborhood((10, 5), category=CategoryRule(:BHID, maxper=3))
    @test length(one.category) == 1
    many = SearchNeighborhood((10, 5), category=[CategoryRule(:BHID, maxper=3), CategoryRule(:ROCK, match=:block)])
    @test length(many.category) == 2
  end

  @testset "validation" begin
    @test_throws ArgumentError SearchNeighborhood((100,))
    @test_throws ArgumentError SearchNeighborhood((100, 50, 20, 10))
    @test_throws ArgumentError SearchNeighborhood((100, 0, 20))
    @test_throws ArgumentError SearchNeighborhood((100, -50, 20))
    @test_throws ArgumentError SearchNeighborhood((10, 5), minsamples=0)
    @test_throws ArgumentError SearchNeighborhood((10, 5), minsamples=10, maxsamples=4)
    @test_throws ArgumentError SearchNeighborhood((10, 5), minpersector=0)
    @test_throws ArgumentError SearchNeighborhood((10, 5), maxpersector=0)

    # more sectors demanded than the scheme has
    @test_throws ArgumentError SearchNeighborhood((10, 5, 2), sectors=Octants(), minsectors=9)
    @test SearchNeighborhood((10, 5, 2), sectors=Octants(), minsectors=8) isa SearchNeighborhood

    # a rule that can never fire is a mistake, not a no-op
    @test_throws ArgumentError SearchNeighborhood((10, 5, 2), sectors=Octants(), maxemptyconsecutive=8)
    @test SearchNeighborhood((10, 5, 2), sectors=Octants(), maxemptyconsecutive=7) isa SearchNeighborhood

    # mutually unsatisfiable quota combinations are caught up front
    @test_throws ArgumentError SearchNeighborhood(
      (10, 5, 2), sectors=Octants(), minsectors=8, minpersector=2, maxsamples=10
    )
    @test_throws ArgumentError SearchNeighborhood(
      (10, 5, 2), sectors=Octants(), maxpersector=1, minsamples=9
    )
  end

  @testset "local projection" begin
    # axis-aligned: the projection just divides by the radii
    s = SearchNeighborhood((100, 50, 20))
    P = localprojection(s, u"m")
    @test P * SVector(100.0, 0.0, 0.0) ≈ SVector(1.0, 0.0, 0.0)
    @test P * SVector(0.0, 25.0, 0.0) ≈ SVector(0.0, 0.5, 0.0)
    @test norm(P * SVector(0.0, 0.0, 20.0)) ≈ 1

    # rotated: points on the rotated semi-axes still land on the unit sphere
    R = RotZ(π / 6)
    r = SearchNeighborhood((100, 50, 20), rotation=R)
    Pr = localprojection(r, u"m")
    for (i, len) in enumerate((100.0, 50.0, 20.0))
      axis = R * SVector(ntuple(j -> j == i ? len : 0.0, 3))
      @test norm(Pr * axis) ≈ 1
    end

    # and the projection agrees with the MetricBall's Mahalanobis metric,
    # which is what the tree will actually index on
    m = metric(metricball(r))
    origin = SVector(0.0, 0.0, 0.0)
    for d in (SVector(30.0, -10.0, 5.0), SVector(-80.0, 40.0, -12.0), SVector(1.0, 1.0, 1.0))
      @test evaluate(m, origin, d) ≈ norm(Pr * d)
    end

    # unit conversion: radii in km, offsets asked for in m
    k = SearchNeighborhood((1u"km", 1u"km", 1u"km"))
    @test norm(localprojection(k, u"m") * SVector(1000.0, 0.0, 0.0)) ≈ 1
    @test norm(localprojection(k, u"km") * SVector(1.0, 0.0, 0.0)) ≈ 1
  end

  @testset "scale" begin
    s = SearchNeighborhood((100, 50, 20), sectors=Octants(), maxpersector=3)
    t = GeoNeighborhoods.scale(s, 2)
    @test t.radii == (200.0u"m", 100.0u"m", 40.0u"m")
    @test t.sectors === s.sectors
    @test t.maxpersector == s.maxpersector
    @test_throws ArgumentError GeoNeighborhoods.scale(s, 0)
    @test_throws ArgumentError GeoNeighborhoods.scale(s, -1)
  end

  @testset "show" begin
    str = sprint(show, SearchNeighborhood((100, 50, 20), sectors=Octants(), maxpersector=3, minsectors=2))
    @test occursin("Octants", str)
    @test occursin("≤3/sector", str)
    @test occursin("≥2 sectors", str)
  end
end
