@testset "estimate" begin
  data = drillholes()
  grid = CartesianGrid((-75.0, -75.0, -40.0), (75.0, 75.0, 40.0), dims=(6, 6, 4))
  model = Kriging(SphericalVariogram(range=120.0u"m"))

  @testset "reproduces GeoStatsModels.fitpredict when unconstrained" begin
    # A ball small enough that everything inside it fits within maxsamples, so
    # neither search has to break a tie -- both take exactly the same samples.
    # The estimates must then agree, which is the whole claim of this package:
    # it changes which samples reach kriging and nothing else.
    ball = (30, 30, 30)
    s = SearchNeighborhood(ball, maxsamples=60)

    ours = interpolate(data, grid, model; search=s, vars=(:Au,))
    theirs = GeoStatsModels.fitpredict(
      model,
      georef((Au=getproperty(values(data), :Au),), domain(data)),
      grid;
      neighbors=true,
      minneighbors=1,
      maxneighbors=60,
      neighborhood=metricball(s)
    )

    a = getproperty(values(ours), :Au)
    b = getproperty(values(theirs), :Au)
    @test length(a) == length(b) == nelements(grid)
    @test ismissing.(a) == ismissing.(b)

    both = .!ismissing.(a)
    @test any(both)                       # the comparison is not vacuous
    @test all(isapprox.(a[both], b[both], rtol=1e-6))
  end

  @testset "the neighborhood changes the answer" begin
    ball = (95, 95, 62)
    plain = interpolate(data, grid, model; search=SearchNeighborhood(ball, maxsamples=12), vars=(:Au,))
    spread = interpolate(
      data, grid, model;
      search=SearchNeighborhood(
        ball, sectors=Octants(), maxsamples=12, maxpersector=2,
        category=CategoryRule(:BHID, maxper=3)
      ),
      vars=(:Au,)
    )

    a = getproperty(values(plain), :Au)
    b = getproperty(values(spread), :Au)
    both = .!ismissing.(a) .& .!ismissing.(b)
    @test any(both)
    @test !all(isapprox.(a[both], b[both], rtol=1e-8))
  end

  @testset "rejected locations are missing, not guessed" begin
    strict = SearchNeighborhood((20, 20, 20), sectors=Octants(), minsectors=8, minsamples=8)
    est = interpolate(data, grid, model; search=strict, vars=(:Au,), diagnostics=true)
    tab = values(est)

    au = getproperty(tab, :Au)
    reject = getproperty(tab, :reject)
    nsamples = getproperty(tab, :nsamples)

    @test all(ismissing.(au) .== (reject .≠ GeoNeighborhoods.Accepted))
    @test all(nsamples[reject .≠ GeoNeighborhoods.Accepted] .== 0)
    @test any(ismissing, au)              # the rule really does bite
  end

  @testset "diagnostics describe each location" begin
    s = SearchNeighborhood(
      (95, 95, 62), sectors=Octants(), maxsamples=12, maxpersector=2,
      category=CategoryRule(:BHID, maxper=3)
    )
    est = interpolate(data, grid, model; search=s, vars=(:Au,), diagnostics=true)
    tab = values(est)

    for col in (:pass, :nsamples, :nsectors, :ndistinct, :reject)
      @test hasproperty(tab, col)
    end
    @test all(getproperty(tab, :pass) .== 1)
    @test all(getproperty(tab, :nsamples) .≤ 12)
    @test all(getproperty(tab, :nsectors) .≤ 8)

    # and they agree with searching that location directly
    m = NeighborhoodSearch(data, s)
    r = searchreport(centroid(grid, 1), m)
    @test getproperty(tab, :nsamples)[1] == length(r.indices)
    @test getproperty(tab, :nsectors)[1] == r.nsectors
    @test getproperty(tab, :ndistinct)[1] == r.ndistinct

    # without the flag the output carries only the variables
    bare = values(interpolate(data, grid, model; search=s, vars=(:Au,)))
    @test Tables.columnnames(Tables.columns(bare)) == (:Au,)
  end

  @testset "multipass records which pass fired" begin
    m = MultiPass(
      SearchNeighborhood((20, 20, 20), minsamples=8),
      SearchNeighborhood((45, 45, 45), minsamples=5),
      SearchNeighborhood((120, 120, 90), minsamples=2)
    )
    est = interpolate(data, grid, model; search=m, vars=(:Au,), diagnostics=true)
    tab = values(est)

    pass = getproperty(tab, :pass)
    @test all(1 .≤ pass .≤ 3)
    @test length(unique(pass)) > 1        # the passes are actually being used

    # a tighter first pass can only shift work later, never earlier
    wide = interpolate(data, grid, model; search=MultiPass(m[2], m[3]), vars=(:Au,), diagnostics=true)
    @test count(ismissing, getproperty(values(wide), :Au)) ≤ count(ismissing, getproperty(tab, :Au))
  end

  @testset "hard domain matching uses the target's own code" begin
    # tag each block with the domain it sits in, the way a block model does
    codes = [ustrip(coords(centroid(grid, i)).x) < 0 ? "WEST" : "EAST" for i in 1:nelements(grid)]
    blocks = georef((ROCK=codes,), grid)

    s = SearchNeighborhood(
      (95, 95, 62), maxsamples=12, category=CategoryRule(:ROCK, match=:block)
    )
    est = interpolate(data, blocks, model; search=s, vars=(:Au,))
    @test nelements(domain(est)) == nelements(grid)

    # a bare domain cannot answer "which rock am I in?"
    @test_throws ArgumentError interpolate(data, grid, model; search=s, vars=(:Au,))
  end

  @testset "options" begin
    s = SearchNeighborhood((95, 95, 62), maxsamples=12)

    # unknown variables are caught rather than silently dropped
    @test_throws ArgumentError interpolate(data, grid, model; search=s, vars=(:Ag,))

    # a diagnostics column that would shadow a variable is refused
    shadow = georef((pass=fill(1.0, nelements(domain(data))),), domain(data))
    @test_throws ArgumentError interpolate(shadow, grid, model; search=s, diagnostics=true)

    # probabilistic prediction returns distributions
    est = interpolate(data, grid, model; search=s, vars=(:Au,), prob=true)
    au = getproperty(values(est), :Au)
    @test all(x -> ismissing(x) || x isa Normal, au)
  end
end
