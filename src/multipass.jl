"""
    MultiPass(pass₁, pass₂, ...)
    MultiPass(passes)
    MultiPass(base, factors; minsamples=nothing)

An ordered sequence of [`SearchNeighborhood`](@ref)s. Each pass is a complete
neighbourhood in its own right, not merely a wider radius: it carries its own
sectors, quotas and rejection rules.

The first pass whose constraints are satisfied wins, and the report says which
one that was. A location no pass can inform is left uninterpolated.

Passes are tried in the order given — put the tightest first.

## Expanding a base neighbourhood

`MultiPass(base, factors)` scales `base` by each factor in turn, which covers
the common case of relaxing the same neighbourhood outwards. Pass `minsamples`
to relax the sample floor alongside the radii:

```julia
MultiPass(base, (1, 2, 3), minsamples=(8, 5, 2))
```

## Example

```julia
MultiPass(
  SearchNeighborhood(( 60, 60, 26), sectors=Octants(), minsamples=8),
  SearchNeighborhood((110, 110, 48), sectors=Octants(), minsamples=5),
  SearchNeighborhood((180, 180, 78), minsamples=2),
)
```
"""
struct MultiPass{P<:Tuple}
  passes::P

  function MultiPass(passes::P) where {P<:Tuple}
    isempty(passes) && throw(ArgumentError("MultiPass needs at least one pass"))
    all(p -> p isa SearchNeighborhood, passes) ||
      throw(ArgumentError("every pass must be a SearchNeighborhood"))
    allequal(ndims(p) for p in passes) ||
      throw(ArgumentError("all passes must have the same dimension, got $(ndims.(passes))"))
    new{P}(passes)
  end
end

MultiPass(passes::AbstractVector) = MultiPass(Tuple(passes))
MultiPass(pass::SearchNeighborhood, rest::SearchNeighborhood...) = MultiPass((pass, rest...))

function MultiPass(base::SearchNeighborhood, factors; minsamples=nothing)
  fs = Tuple(factors)
  isempty(fs) && throw(ArgumentError("MultiPass needs at least one scale factor"))
  if !isnothing(minsamples)
    length(minsamples) == length(fs) ||
      throw(ArgumentError("got $(length(minsamples)) minsamples for $(length(fs)) factors"))
  end

  passes = map(enumerate(fs)) do (i, f)
    p = scale(base, f)
    isnothing(minsamples) ? p : _withminsamples(p, minsamples[i])
  end
  MultiPass(passes)
end

function _withminsamples(s::SearchNeighborhood, n::Int)
  SearchNeighborhood(
    s.radii, s.rotation, s.sectors, n, s.maxsamples,
    s.minpersector, s.maxpersector, s.minsectors, s.maxemptyconsecutive,
    s.split, s.category
  )
end

Base.length(m::MultiPass) = length(m.passes)
Base.getindex(m::MultiPass, i) = m.passes[i]
Base.iterate(m::MultiPass, s...) = iterate(m.passes, s...)
Base.ndims(m::MultiPass) = ndims(first(m.passes))

function Base.show(io::IO, m::MultiPass)
  print(io, "MultiPass(")
  for (i, p) in enumerate(m.passes)
    i > 1 && print(io, ", ")
    print(io, "pass ", i, ": ", p)
  end
  print(io, ")")
end

"""
    MultiPassSearch

A [`MultiPass`](@ref) compiled against data. Build one with
`NeighborhoodSearch(data, multipass)`.

Like [`NeighborhoodSearch`](@ref), it is a `Meshes.BoundedNeighborSearchMethod`
and is safe to share across threads.
"""
struct MultiPassSearch{S<:Tuple} <: BoundedNeighborSearchMethod
  searches::S
end

NeighborhoodSearch(data, m::MultiPass) =
  MultiPassSearch(map(p -> NeighborhoodSearch(data, p), m.passes))

maxneighbors(m::MultiPassSearch) = maximum(maxneighbors, m.searches)
_lenunit(m::MultiPassSearch) = _lenunit(first(m.searches))

npasses(m::MultiPassSearch) = length(m.searches)

function searchdists!(neighbors, distances, pₒ::Point, m::MultiPassSearch; mask=nothing)
  searchreport!(neighbors, distances, pₒ, m; mask).n
end

function searchreport!(neighbors, distances, pₒ::Point, m::MultiPassSearch; mask=nothing, blockvals=nothing)
  # a rejected pass may have written into the buffers; the next pass overwrites
  # what it needs, and a report of n = 0 makes the rest unreadable by contract
  outcome = _report(0, NoCandidates, length(m.searches))
  for (p, s) in enumerate(m.searches)
    r = searchreport!(neighbors, distances, pₒ, s; mask, blockvals)
    r.reason == Accepted && return _report(r.n, r.reason, p, r.nsectors, r.ndistinct)
    outcome = _report(0, r.reason, p, r.nsectors, r.ndistinct)
  end
  # the widest pass gives the most informative failure
  outcome
end

function Base.show(io::IO, m::MultiPassSearch)
  print(io, "MultiPassSearch(", npasses(m), " passes)")
end
