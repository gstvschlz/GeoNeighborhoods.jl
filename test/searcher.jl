@testset "searcher" begin
  data = drillholes()
  origin = Point(0.0, 0.0, 0.0)

  @testset "agrees with Meshes on the unconstrained case" begin
    # with no sectors and no rules the search is "k nearest inside the
    # ellipsoid", which is exactly what KBallSearch computes
    s = SearchNeighborhood((95, 95, 62), maxsamples=12)
    ours = NeighborhoodSearch(data, s)
    theirs = KBallSearch(domain(data), 12, metricball(s))

    r = searchreport(origin, ours)
    want, wantd = searchdists(origin, theirs)
    @test r.reason == Neighborhoods.Accepted
    @test r.pass == 1
    @test length(r.indices) == length(want)

    # composites sit symmetrically about z = 0, so the twelfth place is a tie
    # between the samples above and below. Compare the distance multiset --
    # which sample wins a tie is not part of either contract.
    @test sort(ustrip.(r.distances)) ≈ sort(ustrip.(wantd))
    @test length(symdiff(Set(r.indices), Set(want))) ≤ 2
  end

  @testset "sector and hole quotas spread the selection" begin
    ball = (95, 95, 62)

    plain = NeighborhoodSearch(data, SearchNeighborhood(ball, maxsamples=12))
    holes = tally(data, searchreport(origin, plain).indices, :BHID)
    # the pathology: everything comes from the hole through the block
    @test length(holes) == 1
    @test first(values(holes)) == 12

    spread = NeighborhoodSearch(
      data,
      SearchNeighborhood(
        ball,
        sectors=Octants(),
        maxsamples=12,
        maxpersector=2,
        category=CategoryRule(:BHID, maxper=3)
      )
    )
    r = searchreport(origin, spread)
    @test r.reason == Neighborhoods.Accepted
    @test length(r.indices) == 12

    holes = tally(data, r.indices, :BHID)
    @test length(holes) ≥ 4                 # informed from several holes
    @test maximum(values(holes)) ≤ 3        # and none of them dominates
  end

  @testset "Azimuthal cannot separate above from below" begin
    # A tall, narrow ellipsoid reaching only the hole through the block, so the
    # sector scheme is the only thing deciding how many samples are usable.
    # Azimuthal measures azimuth about the vertical: every composite in a
    # vertical hole shares one wedge, so a per-sector cap throttles the hole as
    # a whole. halves=true and Octants both split above from below.
    withsectors(x) =
      NeighborhoodSearch(data, SearchNeighborhood((30, 30, 60), sectors=x, maxsamples=24, maxpersector=2))

    flat = searchreport(origin, withsectors(Azimuthal(8))).indices
    @test length(flat) == 2
    @test length(tally(data, flat, :BHID)) == 1     # one hole is all there is

    halved = searchreport(origin, withsectors(Azimuthal(8, halves=true))).indices
    @test length(halved) == 4

    oct = searchreport(origin, withsectors(Octants())).indices
    @test length(oct) == 4

    # the split is genuinely above/below, not an arbitrary relabelling
    zs = [_xyz(data, i)[3] for i in halved]
    @test count(≥(0), zs) == 2
    @test count(<(0), zs) == 2
  end

  @testset "quotas are never exceeded" begin
    s = SearchNeighborhood(
      (95, 95, 62),
      sectors=Octants(),
      maxsamples=20,
      maxpersector=2,
      category=CategoryRule(:BHID, maxper=3)
    )
    got = searchreport(origin, NeighborhoodSearch(data, s)).indices

    @test length(got) ≤ 20
    @test all(≤(3), values(tally(data, got, :BHID)))

    proj = localprojection(s, u"m")
    secs = [sectorid(Octants(), proj * SVector(_xyz(data, i)...)) for i in got]
    @test all(≤(2), values(_counts(secs)))
  end

  @testset "distances are sorted and normalised" begin
    m = NeighborhoodSearch(data, SearchNeighborhood((95, 95, 62), maxsamples=10))
    dists = searchreport(origin, m).distances
    @test issorted(dists)
    # normalised: 1.0 is the ellipsoid surface, so everything accepted is inside
    @test all(d -> ustrip(d) ≤ 1, dists)
  end

  @testset "rejection rules" begin
    # nothing within reach at all
    far = NeighborhoodSearch(data, SearchNeighborhood((5, 5, 5)))
    rep = searchreport!(
      Vector{Int}(undef, 24),
      Vector{typeof(1.0u"m")}(undef, 24),
      Point(10_000.0, 10_000.0, 10_000.0),
      far
    )
    @test rep.n == 0
    @test rep.reason == Neighborhoods.NoCandidates

    # samples exist, but not enough of them
    few = NeighborhoodSearch(data, SearchNeighborhood((10, 10, 10), minsamples=50, maxsamples=60))
    r = searchreport(origin, few)
    @test isempty(r.indices)
    @test r.reason == Neighborhoods.TooFewSamples

    # too few distinct sectors filled: a ball reaching only the nearest hole
    # can populate just the two octants above and below it
    sect = NeighborhoodSearch(
      data, SearchNeighborhood((30, 30, 30), sectors=Octants(), minsectors=8, maxsamples=24)
    )
    r = searchreport(origin, sect)
    @test isempty(r.indices)
    @test r.reason == Neighborhoods.TooFewSectors

    # a run of empty sectors around the block
    run = NeighborhoodSearch(
      data, SearchNeighborhood((30, 30, 30), sectors=Azimuthal(16), maxemptyconsecutive=1, maxsamples=24)
    )
    r = searchreport(origin, run)
    @test isempty(r.indices)
    @test r.reason == Neighborhoods.EmptySectorRun

    # one side of the plane empty: search above the highest composite so that
    # every sample in reach lies below it
    split = NeighborhoodSearch(data, SearchNeighborhood((30, 30, 30), split=HalfSpace((0, 0, 1))))
    r = searchreport(Point(3.0, 2.0, 92.0), split)
    @test isempty(r.indices)
    @test r.reason == Neighborhoods.HalfSpaceEmpty

    # a category minimum the data cannot meet
    lonely = NeighborhoodSearch(
      data, SearchNeighborhood((12, 12, 12), category=CategoryRule(:BHID, mindistinct=3))
    )
    r = searchreport(origin, lonely)
    @test isempty(r.indices)
    @test r.reason == Neighborhoods.CategoryUnmet
  end

  @testset "hard domain matching" begin
    s = SearchNeighborhood((95, 95, 62), maxsamples=12, category=CategoryRule(:ROCK, match=:block))
    m = NeighborhoodSearch(data, s)

    # the block's own value has to come from somewhere
    @test_throws ArgumentError searchreport(origin, m)

    r = searchreport(origin, m; blockvals=(; ROCK="WEST"))
    @test r.reason == Neighborhoods.Accepted
    @test all(==("WEST"), keys(tally(data, r.indices, :ROCK)))

    r = searchreport(origin, m; blockvals=(; ROCK="EAST"))
    @test all(==("EAST"), keys(tally(data, r.indices, :ROCK)))
  end

  @testset "reports describe what was assembled" begin
    s = SearchNeighborhood(
      (95, 95, 62), sectors=Octants(), maxsamples=12, maxpersector=2,
      category=CategoryRule(:BHID, maxper=3)
    )
    r = searchreport(origin, NeighborhoodSearch(data, s))
    @test r.nsectors == length(_counts([
      sectorid(Octants(), localprojection(s, u"m") * SVector(_xyz(data, i)...)) for i in r.indices
    ]))
    @test r.ndistinct == length(tally(data, r.indices, :BHID))

    # a rejected search still explains what it managed to gather
    strict = SearchNeighborhood((30, 30, 30), sectors=Octants(), minsectors=8)
    r = searchreport(origin, NeighborhoodSearch(data, strict))
    @test r.reason == Neighborhoods.TooFewSectors
    @test isempty(r.indices)
    @test 0 < r.nsectors < 8

    # no category rules means nothing to count
    plain = searchreport(origin, NeighborhoodSearch(data, SearchNeighborhood((30, 30, 30))))
    @test plain.ndistinct == 0
  end

  @testset "construction errors" begin
    # a 2D neighbourhood cannot search 3D data
    @test_throws ArgumentError NeighborhoodSearch(data, SearchNeighborhood((10, 10)))
    # category rules need a table
    @test_throws ArgumentError NeighborhoodSearch(
      domain(data), SearchNeighborhood((10, 10, 10), category=CategoryRule(:BHID, maxper=2))
    )
    # and the column has to exist
    @test_throws ArgumentError NeighborhoodSearch(
      data, SearchNeighborhood((10, 10, 10), category=CategoryRule(:NOPE, maxper=2))
    )
  end
end
