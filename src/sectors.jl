"""
    SectorScheme

How the search ellipsoid is divided into angular sectors.

Sectors are always assigned in the *rotated, radius-normalised* frame — the
frame in which the ellipsoid is a unit sphere — so that sector geometry and
anisotropy can never disagree.

See [`NoSectors`](@ref), [`Octants`](@ref) and [`Azimuthal`](@ref).
"""
abstract type SectorScheme end

"""
    NoSectors()

No angular subdivision: every sample falls in a single sector. Per-sector
quotas then have no effect, and the search reduces to a plain bounded search
inside the ellipsoid.
"""
struct NoSectors <: SectorScheme end

"""
    Octants()

Divide by the sign of each rotated coordinate: 4 quadrants in 2D, 8 octants
in 3D. This is what GSLIB `kt3d` does. Assignment is branch-free and needs no
trigonometry.

!!! note "Consecutive empty sectors"
    Sign-based sectors have no natural cyclic order, so there is no meaningful
    notion of *adjacent* octants. Under this scheme `maxemptyconsecutive`
    degrades to a limit on the **total** number of empty sectors. Use
    [`Azimuthal`](@ref) if you need the run-based rule.
"""
struct Octants <: SectorScheme end

"""
    Azimuthal(n; halves=false)

Divide into `n` wedges by azimuth around the principal axis, measured in the
normalised local frame. Any `n ≥ 2` is allowed.

In 3D the wedges are vertical: a sample directly above the block lands in the
wedge its horizontal projection falls into. Pass `halves=true` to additionally
split by the sign of the third axis, giving `2n` sectors — the upper half in
`1:n` and the lower half in `n+1:2n`. In 2D there is no third axis and `halves`
is ignored.
"""
struct Azimuthal <: SectorScheme
  n::Int
  halves::Bool

  function Azimuthal(n::Int, halves::Bool)
    n ≥ 2 || throw(ArgumentError("Azimuthal needs at least 2 sectors, got $n"))
    new(n, halves)
  end
end

Azimuthal(n::Int; halves::Bool=false) = Azimuthal(n, halves)

"""
    nsectors(scheme, dim)

Number of sectors the `scheme` produces in `dim` dimensions.
"""
nsectors(::NoSectors, dim::Int) = 1
nsectors(::Octants, dim::Int) = 1 << dim
nsectors(s::Azimuthal, dim::Int) = (s.halves && dim == 3) ? 2s.n : s.n

"""
    sectorid(scheme, u)

Sector of the normalised local coordinates `u`, as an integer in
`1:nsectors(scheme, length(u))`.

A sample lying exactly on the block centre has no direction; it is assigned to
sector 1.
"""
sectorid(::NoSectors, u::AbstractVector) = 1

function sectorid(::Octants, u::AbstractVector)
  id = 1
  @inbounds for i in eachindex(u)
    u[i] < 0 && (id += 1 << (i - 1))
  end
  id
end

function sectorid(s::Azimuthal, u::AbstractVector)
  w = _wedge(u[1], u[2], s.n)
  (s.halves && length(u) == 3 && u[3] < 0) ? w + s.n : w
end

function _wedge(x, y, n)
  φ = atan(y, x)
  φ < 0 && (φ += 2π)
  min(n, floor(Int, n * φ / 2π) + 1)
end

"""
    sectorgroups(scheme, dim)

Ranges of sector ids that are cyclically ordered, i.e. within which "adjacent"
is meaningful and a run of empty sectors can wrap around. Returns an empty
vector for schemes with no cyclic order.
"""
sectorgroups(::NoSectors, dim::Int) = UnitRange{Int}[]
sectorgroups(::Octants, dim::Int) = UnitRange{Int}[]
function sectorgroups(s::Azimuthal, dim::Int)
  (s.halves && dim == 3) ? [1:(s.n), (s.n + 1):(2s.n)] : [1:(s.n)]
end

"""
    consecutiveempty(scheme, dim, counts)

Longest run of consecutive empty sectors, given per-sector sample `counts`.

For cyclic schemes the run may wrap around, and each half of a `halves=true`
[`Azimuthal`](@ref) wraps independently. For schemes with no cyclic order the
result degrades to the **total** number of empty sectors — see [`Octants`](@ref).
"""
function consecutiveempty(scheme::SectorScheme, dim::Int, counts::AbstractVector{Int})
  groups = sectorgroups(scheme, dim)
  isempty(groups) && return count(iszero, counts)
  maximum(_longestzerorun(view(counts, g)) for g in groups)
end

"longest circular run of zeros in `c`"
function _longestzerorun(c)
  n = length(c)
  all(iszero, c) && return n
  best = 0
  run = 0
  # two laps so a run spanning the wrap point is counted once, in full
  @inbounds for i in 1:(2n)
    if iszero(c[mod1(i, n)])
      run += 1
      run > best && (best = run)
    else
      run = 0
    end
  end
  min(best, n)
end
