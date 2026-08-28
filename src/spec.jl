# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const UNBOUNDED = typemax(Int)

"""
    SearchNeighborhood(radii; kwargs...)
    SearchNeighborhood(; radii, kwargs...)

A declarative description of which samples may inform a location.

This is **plain data** — it holds no tree and no state, so it can be built from
an external parameter set, printed, compared and serialised. Compile it against
data with [`NeighborhoodSearch`](@ref) to actually search.

## Geometry

- `radii` — semi-axes of the ellipsoid, as a 2- or 3-tuple. Plain numbers are
  interpreted as metres; `Unitful` lengths are taken as given.
- `rotation` — orientation of the ellipsoid. Accepts `I` (axis-aligned, the
  default), any rotation from Rotations.jl, or a mining convention from
  GeoStatsBase.jl such as `GslibAngles`, `DatamineAngles`, `VulcanAngles` or
  `MinesightAngles`.
- `sectors` — angular subdivision; see [`SectorScheme`](@ref).

## Quotas

- `minsamples` / `maxsamples` — bounds on the total accepted (default `1` / `24`).
- `minpersector` — how many samples a sector needs before it counts as *filled*
  (default `1`).
- `maxpersector` — cap per sector (default unbounded).
- `minsectors` — how many sectors must be filled (default `0`, no requirement).

## Rejection rules

A location that fails any of these is not interpolated:

- fewer than `minsamples` samples found;
- fewer than `minsectors` filled sectors;
- more than `maxemptyconsecutive` consecutive empty sectors — note the degraded
  meaning under [`Octants`](@ref);
- `split` set and one side of the plane empty; see [`HalfSpace`](@ref);
- a [`CategoryRule`](@ref) minimum unmet.

## Example

```julia
SearchNeighborhood(
  (100, 50, 20),
  rotation     = GslibAngles(30, 0, 0),
  sectors      = Octants(),
  minsamples   = 4, maxsamples   = 24,
  minpersector = 1, maxpersector = 3,
  minsectors   = 3,
  maxemptyconsecutive = 2,
  split    = HalfSpace((0, 0, 1)),
  category = CategoryRule(:BHID, maxper=3, mindistinct=2),
)
```
"""
struct SearchNeighborhood{Dim,ℒ<:Unitful.Length,R,S<:SectorScheme,H,C<:Tuple}
  radii::NTuple{Dim,ℒ}
  rotation::R
  sectors::S
  minsamples::Int
  maxsamples::Int
  minpersector::Int
  maxpersector::Int
  minsectors::Int
  maxemptyconsecutive::Int
  split::H
  category::C
end

function SearchNeighborhood(
  radii;
  rotation=I,
  sectors::SectorScheme=NoSectors(),
  minsamples::Int=1,
  maxsamples::Int=24,
  minpersector::Int=1,
  maxpersector=nothing,
  minsectors::Int=0,
  maxemptyconsecutive=nothing,
  split=nothing,
  category=()
)
  # promote so that mixed units (km alongside m) share one representation
  rs = promote(_aslen.(Tuple(radii))...)
  Dim = length(rs)

  Dim ∈ (2, 3) || throw(ArgumentError("radii must have 2 or 3 entries, got $Dim"))
  all(>(zero(first(rs))), rs) || throw(ArgumentError("all radii must be positive, got $radii"))

  maxper = isnothing(maxpersector) ? UNBOUNDED : maxpersector
  maxempty = isnothing(maxemptyconsecutive) ? UNBOUNDED : maxemptyconsecutive
  cats = category isa CategoryRule ? (category,) : Tuple(category)

  minsamples ≥ 1 || throw(ArgumentError("minsamples must be ≥ 1, got $minsamples"))
  maxsamples ≥ minsamples ||
    throw(ArgumentError("maxsamples ($maxsamples) must be ≥ minsamples ($minsamples)"))
  minpersector ≥ 1 || throw(ArgumentError("minpersector must be ≥ 1, got $minpersector"))
  maxper ≥ 1 || throw(ArgumentError("maxpersector must be ≥ 1, got $maxpersector"))
  minsectors ≥ 0 || throw(ArgumentError("minsectors must be ≥ 0, got $minsectors"))
  maxempty ≥ 0 || throw(ArgumentError("maxemptyconsecutive must be ≥ 0, got $maxemptyconsecutive"))

  n = nsectors(sectors, Dim)
  minsectors ≤ n ||
    throw(ArgumentError("minsectors ($minsectors) exceeds the $n sectors of $(sectors)"))
  if maxempty ≠ UNBOUNDED && maxempty ≥ n
    throw(ArgumentError("maxemptyconsecutive ($maxempty) must be < the $n sectors of $(sectors), otherwise it can never trigger"))
  end
  if minsectors > 0 && minsectors * minpersector > maxsamples
    throw(ArgumentError("minsectors × minpersector ($(minsectors * minpersector)) exceeds maxsamples ($maxsamples): no selection can satisfy both"))
  end
  if maxper ≠ UNBOUNDED && n * maxper < minsamples
    throw(ArgumentError("maxpersector ($maxper) across $n sectors allows at most $(n * maxper) samples, below minsamples ($minsamples)"))
  end

  SearchNeighborhood(
    rs, rotation, sectors, minsamples, maxsamples,
    minpersector, maxper, minsectors, maxempty, split, cats
  )
end

SearchNeighborhood(; radii, kwargs...) = SearchNeighborhood(radii; kwargs...)

_aslen(x::Unitful.Length) = float(x)
_aslen(x::Number) = float(x) * u"m"

"""
    ndims(spec)

Spatial dimension the neighbourhood is defined in — 2 or 3.
"""
Base.ndims(::SearchNeighborhood{Dim}) where {Dim} = Dim

"""
    nsectors(spec)

Number of sectors the neighbourhood divides into.
"""
nsectors(s::SearchNeighborhood) = nsectors(s.sectors, ndims(s))

"""
    metricball(spec)

The `Meshes.MetricBall` describing the ellipsoid, whose Mahalanobis metric
makes the ellipsoid a unit ball.
"""
metricball(s::SearchNeighborhood) = MetricBall(s.radii, s.rotation)

"""
    localprojection(spec, unit)

Matrix mapping a world-space offset from the block centre, expressed in `unit`,
to normalised local coordinates: rotated into the ellipsoid frame and divided by
the radii. The ellipsoid interior is then exactly the unit ball, so `norm(u) ≤ 1`
tests membership and `u` feeds sector assignment directly.
"""
function localprojection(s::SearchNeighborhood{Dim}, unit) where {Dim}
  R = _rotmatrix(s.rotation, Val(Dim))
  r = SVector{Dim}(ustrip.(unit, s.radii))
  SMatrix{Dim,Dim}(Diagonal(inv.(r)) * transpose(R))
end

_rotmatrix(::UniformScaling, ::Val{Dim}) where {Dim} = SMatrix{Dim,Dim}(1.0I)
_rotmatrix(R::AbstractMatrix, ::Val{Dim}) where {Dim} = SMatrix{Dim,Dim}(R)

"""
    scale(spec, factor)

Scale the ellipsoid radii by a strictly positive `factor`, leaving every other
parameter untouched. Used to rescale a neighbourhood alongside the data for
numerical stability, and to build the passes of a [`MultiPass`](@ref).
"""
function scale(s::SearchNeighborhood, factor::Real)
  factor > 0 || throw(ArgumentError("scale factor must be positive, got $factor"))
  SearchNeighborhood(
    factor .* s.radii, s.rotation, s.sectors, s.minsamples, s.maxsamples,
    s.minpersector, s.maxpersector, s.minsectors, s.maxemptyconsecutive,
    s.split, s.category
  )
end

function Base.show(io::IO, s::SearchNeighborhood)
  print(io, "SearchNeighborhood(", s.radii, ", ", nameof(typeof(s.sectors)))
  print(io, ", ", s.minsamples, "-", s.maxsamples, " samples")
  s.maxpersector == UNBOUNDED || print(io, ", ≤", s.maxpersector, "/sector")
  s.minsectors == 0 || print(io, ", ≥", s.minsectors, " sectors")
  s.maxemptyconsecutive == UNBOUNDED || print(io, ", ≤", s.maxemptyconsecutive, " empty")
  isnothing(s.split) || print(io, ", split")
  isempty(s.category) || print(io, ", ", length(s.category), " category rule(s)")
  print(io, ")")
end
