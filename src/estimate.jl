# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const DIAGNOSTICS = (:pass, :nsamples, :nsectors, :ndistinct, :reject)

"""
    interpolate(data, target, model; search, kwargs...)

Estimate `data` onto `target` with a GeoStatsModels `model`, choosing the
samples for each location with `search` — a [`SearchNeighborhood`](@ref), a
[`MultiPass`](@ref), or an already-compiled searcher.

`target` is a `Domain` (a grid or point set), or a geotable when a
[`CategoryRule`](@ref) uses `match=:block` and the estimated locations carry
their own domain code.

The model is untouched: it is fitted to the samples the search returns and
asked to predict, exactly as `GeoStatsModels.fitpredict` does. **A location
whose search is rejected is left `missing`** rather than estimated from
whatever happened to be nearby.

## Options

- `vars` — which variables to estimate (default: every column of `data`).
  Kriging cannot combine non-numeric columns, so name the grade columns when
  the table also carries hole ids or domain codes.
- `point` — estimate on point support (default `true`).
- `prob` — return distributions instead of means (default `false`).
- `diagnostics` — append `pass`, `nsamples`, `nsectors`, `ndistinct` and
  `reject` columns explaining each location (default `false`).

## Example

```julia
model = Kriging(SphericalVariogram(range=120.0))
search = SearchNeighborhood((100, 100, 40), sectors=Octants(), maxpersector=2)

estimate = interpolate(samples, grid, model; search, vars=(:Au,), diagnostics=true)
```

!!! note "Difference from `GeoStatsModels.fitpredict`"
    `fitpredict` rescales the data, domain and model before solving, for
    numerical conditioning. This function does not, so that it depends only on
    the public model interface. With coordinates far from the origin and a
    near-singular system, prefer rescaling the inputs yourself.
"""
function interpolate(data, target, model; search, vars=nothing, point=true, prob=false, diagnostics=false)
  names = _variables(data, vars)
  # the searcher needs every column the category rules mention, while the model
  # must see only the estimated variables -- same support, same element order
  support = _support(data, point)
  samples = georef(values(data), support)
  dat = _select(data, names, support)
  dom = target isa Domain ? target : domain(target)

  if diagnostics
    clash = intersect(names, DIAGNOSTICS)
    isempty(clash) ||
      throw(ArgumentError("diagnostics would overwrite the variable(s) $(collect(clash)); rename them or pass diagnostics=false"))
  end

  searcher = search isa BoundedNeighborSearchMethod ? search : NeighborhoodSearch(samples, search)
  blockvals = _blockvalues(target, search)

  n = nelements(dom)
  k = maxneighbors(searcher)
  ℒ = typeof(one(Float64) * _lenunit(searcher))

  columns = [Vector{Any}(undef, n) for _ in names]
  pass = Vector{Int}(undef, n)
  nsamples = Vector{Int}(undef, n)
  nsectors = Vector{Int}(undef, n)
  ndistinct = Vector{Int}(undef, n)
  reject = Vector{RejectReason}(undef, n)

  _foreachchunk(1:n) do chunk
    # buffers are per task, so the searcher stays shareable and read-only
    neighbors = Vector{Int}(undef, k)
    distances = Vector{ℒ}(undef, k)

    for ind in chunk
      pₒ = centroid(dom, ind)
      bv = isnothing(blockvals) ? nothing : blockvals(ind)
      r = searchreport!(neighbors, distances, pₒ, searcher; blockvals=bv)

      if r.n > 0
        ndata = view(dat, view(neighbors, 1:r.n))
        fitted = GeoStatsModels.fit(model, ndata)
        geom = point ? pₒ : dom[ind]
        for (j, var) in enumerate(names)
          columns[j][ind] = prob ?
                            GeoStatsModels.predictprob(fitted, var, geom) :
                            GeoStatsModels.predict(fitted, var, geom)
        end
      else
        for j in eachindex(names)
          columns[j][ind] = missing
        end
      end

      pass[ind] = r.pass
      nsamples[ind] = r.n
      nsectors[ind] = r.nsectors
      ndistinct[ind] = r.ndistinct
      reject[ind] = r.reason
    end
  end

  # narrow Vector{Any} to whatever the predictions actually are
  estimated = (; zip(names, map(c -> map(identity, c), columns))...)
  table = diagnostics ? merge(estimated, (; pass, nsamples, nsectors, ndistinct, reject)) : estimated

  georef(table, dom)
end

# ---------------------------------------------------------------- helpers

function _variables(data, vars)
  cols = Tables.columns(values(data))
  available = Tables.columnnames(cols)
  isnothing(vars) && return Tuple(available)
  names = Tuple(Symbol.(vars isa Symbol ? (vars,) : vars))
  for v in names
    v ∈ available ||
      throw(ArgumentError("variable :$v is not in the data; available: $(collect(available))"))
  end
  names
end

"point support collapses each element to its centroid; volume support keeps it"
function _support(data, point)
  dom = domain(data)
  point ? PointSet([centroid(dom, i) for i in 1:nelements(dom)]) : dom
end

"the estimated variables alone, so kriging never meets a hole id or rock code"
function _select(data, names, support)
  cols = Tables.columns(values(data))
  georef((; (v => Tables.getcolumn(cols, v) for v in names)...), support)
end

"""
Returns a function mapping a target index to the block's own category values,
keyed by column, or `nothing` when no rule needs them.
"""
function _blockvalues(target, search)
  rules = _matchrules(search)
  isempty(rules) && return nothing

  target isa Domain && throw(
    ArgumentError(
      "CategoryRule on :$(first(rules).column) uses match=:block, so the target " *
      "must be a geotable carrying that column, not a bare Domain"
    )
  )

  cols = Tables.columns(values(target))
  available = Tables.columnnames(cols)
  for rule in rules
    rule.column ∈ available ||
      throw(ArgumentError("target has no column :$(rule.column); available: $(collect(available))"))
  end

  names = Tuple(unique(rule.column for rule in rules))
  vectors = map(v -> Tables.getcolumn(cols, v), names)
  ind -> (; zip(names, map(c -> c[ind], vectors))...)
end

_matchrules(s::SearchNeighborhood) = filter(r -> r.matchblock, collect(s.category))
_matchrules(m::MultiPass) = reduce(vcat, _matchrules(p) for p in m.passes)
_matchrules(m::NeighborhoodSearch) = _matchrules(m.spec)
_matchrules(m::MultiPassSearch) = reduce(vcat, _matchrules(s) for s in m.searches)

"split `inds` into one chunk per thread and run `f` on each concurrently"
function _foreachchunk(f, inds; nchunks=Threads.nthreads())
  n = length(inds)
  n == 0 && return nothing
  nchunks = clamp(nchunks, 1, n)
  Threads.@sync for c in 1:nchunks
    lo = 1 + ((c - 1) * n) ÷ nchunks
    hi = (c * n) ÷ nchunks
    Threads.@spawn f(view(inds, lo:hi))
  end
  nothing
end
