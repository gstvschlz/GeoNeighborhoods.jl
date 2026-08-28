@testset "constraints" begin
  @testset "HalfSpace" begin
    h = HalfSpace((0, 0, 1))
    @test side(h, [0.0, 0.0, 5.0]) == 1
    @test side(h, [0.0, 0.0, -5.0]) == 2
    # a sample on the plane counts as positive
    @test side(h, [3.0, 4.0, 0.0]) == 1

    # the normal need not be axis-aligned or unit length
    d = HalfSpace([2.0, 2.0, 0.0])
    @test side(d, [1.0, 0.0, 0.0]) == 1
    @test side(d, [-1.0, 0.5, 0.0]) == 2

    @test_throws ArgumentError HalfSpace((0, 0, 0))
  end

  @testset "CategoryRule construction" begin
    r = CategoryRule(:BHID, maxper=3, mindistinct=2)
    @test r.column == :BHID
    @test r.maxper == 3
    @test r.mindistinct == 2
    @test r.matchblock == false
    @test isnothing(r.quotas)

    @test CategoryRule(:ROCK, match=:block).matchblock

    @test_throws ArgumentError CategoryRule(:X, match=:cell)
    @test_throws ArgumentError CategoryRule(:X, maxper=0)
    @test_throws ArgumentError CategoryRule(:X, mindistinct=0)
    @test_throws ArgumentError CategoryRule(:X, quotas=Dict("A" => 5:2))
    @test_throws ArgumentError CategoryRule(:X, quotas=Dict("A" => -1:2))
    @test_throws ArgumentError CategoryRule(:X, quotas=Dict{String,UnitRange{Int}}())
    # a rule that constrains nothing is a mistake, not a no-op
    @test_throws ArgumentError CategoryRule(:X)
  end

  @testset "cap" begin
    @test cap(CategoryRule(:X, maxper=3), "A") == 3
    @test cap(CategoryRule(:X, mindistinct=2), "A") == typemax(Int)

    # the tighter of maxper and the quota maximum wins, in both directions
    r = CategoryRule(:X, maxper=3, quotas=Dict("A" => 0:1, "B" => 0:9))
    @test cap(r, "A") == 1
    @test cap(r, "B") == 3
    @test cap(r, "unlisted") == 3
  end

  @testset "satisfied" begin
    r = CategoryRule(:BHID, mindistinct=2)
    @test satisfied(r, Dict("h1" => 3, "h2" => 1))
    @test !satisfied(r, Dict("h1" => 12))
    # a value present with a zero count does not count as distinct
    @test !satisfied(r, Dict("h1" => 3, "h2" => 0))

    q = CategoryRule(:ZONE, quotas=Dict("A" => 2:6, "B" => 0:4))
    @test satisfied(q, Dict("A" => 2, "B" => 0))
    @test satisfied(q, Dict("A" => 6))          # B has no floor
    @test !satisfied(q, Dict("A" => 1, "B" => 4))
    @test !satisfied(q, Dict("B" => 4))         # A missing entirely

    # no minimum requirements ⇒ always satisfied
    @test satisfied(CategoryRule(:X, maxper=2), Dict("A" => 1))
  end
end
