# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RejectReason

Why a neighbourhood search declined to produce a sample set. `Accepted` means
it did not. Reported per location when diagnostics are requested, so that a
block left `missing` can be explained rather than guessed at.
"""
@enum RejectReason begin
  Accepted
  NoCandidates
  TooFewSamples
  TooFewSectors
  EmptySectorRun
  HalfSpaceEmpty
  CategoryUnmet
end

"""
    HalfSpace(normal)

Split the ellipsoid with the plane through the block centre whose normal is
`normal`, and require samples on both sides. A block informed entirely from one
side is not interpolated.

The normal is given in **world coordinates**, not in the rotated frame, so
`HalfSpace((0, 0, 1))` means "above and below the horizontal plane" regardless
of how the ellipsoid is oriented.

Samples lying exactly on the plane count towards the positive side.
"""
struct HalfSpace{V<:AbstractVector}
  normal::V

  function HalfSpace(normal::V) where {V<:AbstractVector}
    iszero(normal) && throw(ArgumentError("HalfSpace normal must be non-zero"))
    new{V}(normal)
  end
end

HalfSpace(normal::Tuple) = HalfSpace(collect(float.(normal)))

"""
    side(halfspace, d)

Which side of the plane the world-space offset `d` falls on: `1` for the
positive side (including the plane itself), `2` for the negative side.
"""
side(h::HalfSpace, d) = _dot(h.normal, d) ≥ 0 ? 1 : 2

function _dot(a, b)
  s = zero(eltype(a)) * zero(eltype(b))
  @inbounds for i in eachindex(a, b)
    s += a[i] * b[i]
  end
  s
end

"""
    CategoryRule(column; maxper=nothing, mindistinct=nothing, match=nothing, quotas=nothing)

A constraint tying sample selection to a categorical column of the sample
table, such as a drillhole id or a domain code.

- `maxper` — accept at most this many samples sharing any one value. The
  classic "no more than 3 composites from the same hole" rule.
- `mindistinct` — require samples from at least this many different values,
  otherwise do not interpolate.
- `match=:block` — hard domain matching: only accept samples whose value equals
  the value of the location being estimated. Requires the target to carry the
  same column; see [`interpolate`](@ref).
- `quotas` — explicit per-value bounds, as `value => min:max`. The maximum is
  enforced while selecting; the minimum is checked afterwards and rejects the
  location when unmet.

`maxper` and a quota maximum for the same value are both honoured — the tighter
one wins.

## Examples

```julia
CategoryRule(:BHID, maxper=3, mindistinct=2)
CategoryRule(:ROCK, match=:block)
CategoryRule(:ZONE, quotas=Dict("A" => 2:6, "B" => 0:4))
```
"""
struct CategoryRule{Q<:Union{Nothing,AbstractDict}}
  column::Symbol
  maxper::Union{Int,Nothing}
  mindistinct::Union{Int,Nothing}
  matchblock::Bool
  quotas::Q
end

function CategoryRule(column::Symbol; maxper=nothing, mindistinct=nothing, match=nothing, quotas=nothing)
  matchblock = if isnothing(match)
    false
  elseif match === :block
    true
  else
    throw(ArgumentError("match must be :block or nothing, got $(repr(match))"))
  end

  isnothing(maxper) || maxper ≥ 1 || throw(ArgumentError("maxper must be ≥ 1, got $maxper"))
  isnothing(mindistinct) ||
    mindistinct ≥ 1 ||
    throw(ArgumentError("mindistinct must be ≥ 1, got $mindistinct"))

  if !isnothing(quotas)
    isempty(quotas) && throw(ArgumentError("quotas must not be empty"))
    for (value, range) in quotas
      first(range) ≤ last(range) ||
        throw(ArgumentError("quota for $(repr(value)) is empty: $range"))
      first(range) ≥ 0 || throw(ArgumentError("quota for $(repr(value)) must be non-negative: $range"))
    end
  end

  if isnothing(maxper) && isnothing(mindistinct) && isnothing(quotas) && !matchblock
    throw(ArgumentError("CategoryRule on :$column constrains nothing"))
  end

  CategoryRule(column, maxper, mindistinct, matchblock, quotas)
end

"""
    cap(rule, value)

Largest number of samples the `rule` allows for a given category `value`,
combining `maxper` with any per-value quota maximum. `typemax(Int)` when
unbounded.
"""
function cap(rule::CategoryRule, value)
  c = isnothing(rule.maxper) ? typemax(Int) : rule.maxper
  if !isnothing(rule.quotas)
    q = get(rule.quotas, value, nothing)
    isnothing(q) || (c = min(c, last(q)))
  end
  c
end

"""
    satisfied(rule, tally)

Whether the per-value counts in `tally` meet the rule's minimum requirements:
`mindistinct` and any per-value quota floors. Only meaningful once selection
has finished.
"""
function satisfied(rule::CategoryRule, tally::AbstractDict)
  if !isnothing(rule.mindistinct)
    count(>(0), values(tally)) ≥ rule.mindistinct || return false
  end
  if !isnothing(rule.quotas)
    for (value, range) in rule.quotas
      get(tally, value, 0) ≥ first(range) || return false
    end
  end
  true
end
