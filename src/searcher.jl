"""
    NeighborhoodSearch(data, spec)

A [`SearchNeighborhood`](@ref) compiled against `data` — a geotable, or any
`Domain` when the specification carries no [`CategoryRule`](@ref).

This is a `Meshes.BoundedNeighborSearchMethod`, so it can be used anywhere the
ecosystem accepts a searcher, and it is immutable and safe to share across
threads.

Searching returns the accepted sample indices. When a rejection rule fires it
returns **no** neighbours, which is what makes the caller leave the location
uninterpolated.

!!! note "Hard domain matching"
    A `CategoryRule` with `match=:block` needs the value of the location being
    estimated, which the plain `Meshes` search interface has no way to pass.
    Use [`searchreport!`](@ref) with `blockvals=(; COLUMN=value)`, or
    [`interpolate`](@ref), which supplies them from the target geotable.
"""
struct NeighborhoodSearch{Dim,T,C<:CRS,S<:SearchNeighborhood,TR,P,V,U} <: BoundedNeighborSearchMethod
  spec::S
  tree::TR
  proj::P
  coords::Vector{SVector{Dim,T}}
  catvals::V
  unit::U
end

function NeighborhoodSearch(data, spec::SearchNeighborhood{Dim}) where {Dim}
  dom = _domainof(data)
  embeddim(dom) == Dim ||
    throw(ArgumentError("neighbourhood is $(Dim)D but the data is $(embeddim(dom))D"))

  C = crs(dom)
  u = unit(lentype(C))
  n = nelements(dom)

  X = [_svec(centroid(dom, i)) for i in 1:n]
  T = eltype(eltype(X))

  tree = BallTree(X, metric(metricball(spec)))
  proj = localprojection(spec, u)
  catvals = _categorycolumns(data, spec)

  NeighborhoodSearch{Dim,T,C,typeof(spec),typeof(tree),typeof(proj),typeof(catvals),typeof(u)}(
    spec, tree, proj, X, catvals, u
  )
end

_domainof(data::Domain) = data
_domainof(data) = domain(data)

function _categorycolumns(data, spec::SearchNeighborhood)
  isempty(spec.category) && return ()
  data isa Domain &&
    throw(ArgumentError("category rules need a geotable with the referenced columns, got a bare Domain"))
  cols = Tables.columns(values(data))
  names = Tables.columnnames(cols)
  map(spec.category) do rule
    rule.column ∈ names ||
      throw(ArgumentError("column :$(rule.column) is not in the data; available: $(collect(names))"))
    collect(Tables.getcolumn(cols, rule.column))
  end
end

maxneighbors(m::NeighborhoodSearch) = m.spec.maxsamples

"""
    spec(searcher)

The [`SearchNeighborhood`](@ref) the `searcher` was compiled from.
"""
spec(m::NeighborhoodSearch) = m.spec

function searchdists!(neighbors, distances, pₒ::Point, m::NeighborhoodSearch; mask=nothing)
  searchreport!(neighbors, distances, pₒ, m; mask).n
end

"""a search outcome: how many neighbours, why it stopped, and which pass produced it"""
_report(n, reason, pass=1, nsectors=0, ndistinct=0) = (; n, reason, pass, nsectors, ndistinct)

"""
    searchreport!(neighbors, distances, pₒ, method; mask=nothing, blockvals=nothing)

Like `Meshes.searchdists!`, but returns a named tuple `(; n, reason, pass)`
describing how the search ended:

- `n` — how many neighbours were accepted; always `0` unless `reason` is `Accepted`.
- `reason` — see [`RejectReason`](@ref).
- `pass` — which [`MultiPass`](@ref) pass produced this, always `1` for a single
  neighbourhood.
- `nsectors` — how many sectors reached `minpersector`.
- `ndistinct` — how many distinct values of the **first** [`CategoryRule`](@ref)
  column were drawn, or `0` when there are no category rules.

`nsectors` and `ndistinct` describe what was assembled even when a rule then
rejected it, which is what makes a refusal explainable.

`blockvals` supplies the estimated location's own category values, keyed by
column name — for example `(; ROCK="WEST")`. It is required when any rule uses
`match=:block`. Keying by name rather than position means the passes of a
[`MultiPass`](@ref) may carry different rules.
"""
function searchreport!(
  neighbors,
  distances,
  pₒ::Point,
  m::NeighborhoodSearch{Dim,T,C};
  mask=nothing,
  blockvals=nothing
) where {Dim,T,C}
  s = m.spec
  x = _svec(convert(C, coords(pₒ)))

  # the Mahalanobis metric makes the ellipsoid interior the unit ball, so a
  # radius of one is the whole neighbourhood
  cand = inrange(m.tree, x, one(T))
  isempty(cand) && return _report(0, NoCandidates)

  keep = Int[]
  dists = T[]
  secs = Int[]
  sides = Int[]
  @inbounds for i in cand
    isnothing(mask) || mask[i] || continue
    _matchesblock(m, s, i, blockvals) || continue
    d = m.coords[i] - x
    u = m.proj * d
    push!(keep, i)
    push!(dists, norm(u))
    push!(secs, sectorid(s.sectors, u))
    push!(sides, isnothing(s.split) ? 1 : side(s.split, d))
  end
  isempty(keep) && return _report(0, NoCandidates)

  nsec = nsectors(s)
  seccount = zeros(Int, nsec)
  sidecount = zeros(Int, 2)
  tallies = _emptytallies(m.catvals)
  taken = 0

  @inbounds for k in sortperm(dists)
    taken == s.maxsamples && break
    seccount[secs[k]] ≥ s.maxpersector && continue
    _withincaps(m, s, tallies, keep[k]) || continue

    taken += 1
    neighbors[taken] = keep[k]
    # the Mahalanobis distance is dimensionless; it carries the length unit
    # only so that the value matches what Meshes' own ball searches return
    distances[taken] = dists[k] * m.unit
    seccount[secs[k]] += 1
    sidecount[sides[k]] += 1
    _tally!(m, s, tallies, keep[k])
  end

  # reported whether or not a rule then rejects the selection
  nfilled = count(≥(s.minpersector), seccount)
  ndistinct = isempty(tallies) ? 0 : count(>(0), values(first(tallies)))
  reject(reason) = _report(0, reason, 1, nfilled, ndistinct)

  taken < s.minsamples && return reject(TooFewSamples)
  s.minsectors > 0 && nfilled < s.minsectors && return reject(TooFewSectors)
  if s.maxemptyconsecutive ≠ UNBOUNDED
    consecutiveempty(s.sectors, Dim, seccount) > s.maxemptyconsecutive &&
      return reject(EmptySectorRun)
  end
  if !isnothing(s.split)
    (sidecount[1] == 0 || sidecount[2] == 0) && return reject(HalfSpaceEmpty)
  end
  for (r, rule) in enumerate(s.category)
    satisfied(rule, tallies[r]) || return reject(CategoryUnmet)
  end

  _report(taken, Accepted, 1, nfilled, ndistinct)
end

"""
    searchreport(pₒ, method; mask=nothing, blockvals=nothing)

Allocating form of [`searchreport!`](@ref): returns a named tuple
`(; indices, distances, reason, pass, nsectors, ndistinct)`.
"""
function searchreport(pₒ::Point, m::BoundedNeighborSearchMethod; kwargs...)
  k = maxneighbors(m)
  neighbors = Vector{Int}(undef, k)
  distances = Vector{typeof(one(Float64) * _lenunit(m))}(undef, k)
  r = searchreport!(neighbors, distances, pₒ, m; kwargs...)
  (; indices=view(neighbors, 1:r.n), distances=view(distances, 1:r.n), r.reason, r.pass, r.nsectors, r.ndistinct)
end

_lenunit(m::NeighborhoodSearch) = m.unit

# same convention as Meshes' internal svec: plain numbers in the CRS units
_svec(p::Point) = _svec(coords(p))
_svec(c::CRS) = SVector(CoordRefSystems.raw(c))

_emptytallies(catvals) = map(v -> Dict{eltype(v),Int}(), catvals)

function _matchesblock(m, s, i, blockvals)
  @inbounds for (r, rule) in enumerate(s.category)
    rule.matchblock || continue
    (isnothing(blockvals) || !haskey(blockvals, rule.column)) && throw(
      ArgumentError(
        "CategoryRule on :$(rule.column) uses match=:block, which needs the " *
        "estimated location's own value; pass blockvals=(; $(rule.column)=…), " *
        "or use interpolate with a geotable target"
      )
    )
    m.catvals[r][i] == blockvals[rule.column] || return false
  end
  true
end

function _withincaps(m, s, tallies, i)
  @inbounds for (r, rule) in enumerate(s.category)
    v = m.catvals[r][i]
    get(tallies[r], v, 0) ≥ cap(rule, v) && return false
  end
  true
end

function _tally!(m, s, tallies, i)
  @inbounds for r in eachindex(s.category)
    v = m.catvals[r][i]
    tallies[r][v] = get(tallies[r], v, 0) + 1
  end
end

function Base.show(io::IO, m::NeighborhoodSearch)
  print(io, "NeighborhoodSearch(", length(m.coords), " samples, ", m.spec, ")")
end
