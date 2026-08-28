# ------------------------------------------------------------------
# Regenerates the README figures. Plain Julia, no dependencies.
#   julia assets/make_figures.jl
# The geometry mirrors the documented search rules so the pictures
# stay honest as the design evolves.
# ------------------------------------------------------------------

const BG = "#fbfaf8"
const EDGE = "#e6e2db"
const INK = "#2f2c28"
const MUTED = "#8a837a"
const ELLIPSE = "#3f6fb5"
const GRID = "#d8d2c8"
const TAKEN = "#1f7a4d"
const SKIP = "#bdb6ac"
const BLOCK = "#d1442f"

# ---------------------------------------------------------------- geometry

struct Sample
  x::Float64
  z::Float64
  hole::Int
end

"local (rotated) coordinates of a point relative to the block centre"
function localcoords(x, z, θ)
  c, s = cos(θ), sin(θ)
  (c * x + s * z, -s * x + c * z)
end

"normalised distance: 1.0 is exactly on the ellipsoid"
function normdist(x, z, radii, θ)
  lx, lz = localcoords(x, z, θ)
  hypot(lx / radii[1], lz / radii[2])
end

"sector id in 1:n by azimuth in the normalised local frame"
function azimuthal(x, z, radii, θ, n)
  lx, lz = localcoords(x, z, θ)
  φ = atan(lz / radii[2], lx / radii[1])
  φ < 0 && (φ += 2π)
  min(n, floor(Int, n * φ / 2π) + 1)
end

"quadrant id in 1:4 by sign of the local coordinates — the 2D case of Octants"
function quadrant(x, z, radii, θ)
  lx, lz = localcoords(x, z, θ)
  (lx ≥ 0 ? 1 : 2) + (lz ≥ 0 ? 0 : 2)
end

"""
Greedy selection in ascending normalised distance, skipping a candidate when
its sector is full, its hole is capped, or the global maximum is reached.
"""
function select(cands, dists, sectorof; maxtotal, maxpersector, maxperhole)
  order = sortperm(dists)
  persector = Dict{Int,Int}()
  perhole = Dict{Int,Int}()
  taken = Int[]
  for i in order
    length(taken) == maxtotal && break
    s = sectorof(i)
    h = cands[i].hole
    get(persector, s, 0) ≥ maxpersector && continue
    get(perhole, h, 0) ≥ maxperhole && continue
    push!(taken, i)
    persector[s] = get(persector, s, 0) + 1
    perhole[h] = get(perhole, h, 0) + 1
  end
  taken, persector, perhole
end

"""
Steeply dipping holes on a 45 m section spacing, composited every 4 m — the
ordinary case where down-hole sampling is an order of magnitude denser than
hole spacing, which is exactly what makes a plain nearest-N search degenerate.
The middle hole is placed so that it passes through the block being estimated.
"""
function drillholes()
  collars = [-16.72 + 45.0k for k in -3:3]
  out = Sample[]
  for (h, x0) in enumerate(collars)
    for k in 0:38
      z = 76.0 - 4.0k
      x = x0 + 0.22 * (76.0 - z)
      push!(out, Sample(x, z, h))
    end
  end
  out
end

# ---------------------------------------------------------------- svg

r2(v) = string(round(v, digits=2))

esc(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

function text(x, y, s; size=12, fill=INK, anchor="start", weight="400")
  font = "ui-sans-serif,-apple-system,Segoe UI,Roboto,sans-serif"
  string(
    "<text x=\"", r2(x), "\" y=\"", r2(y), "\" font-family=\"", font,
    "\" font-size=\"", size, "\" font-weight=\"", weight,
    "\" fill=\"", fill, "\" text-anchor=\"", anchor, "\">", esc(s), "</text>"
  )
end

struct View
  ox::Float64
  oy::Float64
  scale::Float64
  xspan::Float64
  zspan::Float64
end
sx(v::View, x) = v.ox + (x + v.xspan) * v.scale
sy(v::View, z) = v.oy + (v.zspan - z) * v.scale

"rotate a point given in the local ellipsoid frame back to world coordinates"
function toworld(lx, lz, θ)
  c, s = cos(θ), sin(θ)
  (c * lx - s * lz, s * lx + c * lz)
end

function ellipsepath(v::View, radii, θ)
  d = IOBuffer()
  for (i, t) in enumerate(range(0, 2π, length=181))
    x, z = toworld(radii[1] * cos(t), radii[2] * sin(t), θ)
    print(d, i == 1 ? "M" : "L", r2(sx(v, x)), " ", r2(sy(v, z)), " ")
  end
  print(d, "Z")
  String(take!(d))
end

"dashed ray from the block centre to the ellipsoid boundary at normalised angle φ"
function sectorray(v::View, radii, θ, φ)
  x, z = toworld(radii[1] * cos(φ), radii[2] * sin(φ), θ)
  string(
    "<line x1=\"", r2(sx(v, 0)), "\" y1=\"", r2(sy(v, 0)),
    "\" x2=\"", r2(sx(v, x)), "\" y2=\"", r2(sy(v, z)),
    "\" stroke=\"", GRID, "\" stroke-width=\"1\" stroke-dasharray=\"3 3\"/>"
  )
end

"""
Draws the panel frame and its heading, then opens a clipped group so that
drillhole traces cannot bleed past the frame. Caller must close it with
`endpanel`.
"""
function panel(io, v::View, w, h, title, subtitle; id="clip")
  println(io, "<rect x=\"", r2(v.ox - 14), "\" y=\"", r2(v.oy - 44), "\" width=\"", w,
    "\" height=\"", h, "\" rx=\"8\" fill=\"", BG, "\" stroke=\"", EDGE, "\"/>")
  println(io, text(v.ox - 2, v.oy - 24, title; size=13, weight="600"))
  println(io, text(v.ox - 2, v.oy - 8, subtitle; size=11, fill=MUTED))
  println(io, "<defs><clipPath id=\"", id, "\"><rect x=\"", r2(v.ox - 14), "\" y=\"", r2(v.oy),
    "\" width=\"", w, "\" height=\"", r2(h - 46), "\" rx=\"6\"/></clipPath></defs>")
  println(io, "<g clip-path=\"url(#", id, ")\">")
end

endpanel(io) = println(io, "</g>")

function drawholes(io, v, samples)
  for h in unique(s.hole for s in samples)
    pts = filter(s -> s.hole == h, samples)
    a, b = first(pts), last(pts)
    println(io, "<line x1=\"", r2(sx(v, a.x)), "\" y1=\"", r2(sy(v, a.z)),
      "\" x2=\"", r2(sx(v, b.x)), "\" y2=\"", r2(sy(v, b.z)),
      "\" stroke=\"", GRID, "\" stroke-width=\"1\"/>")
  end
end

function dot(io, v, s, r, fill; opacity=1.0)
  println(io, "<circle cx=\"", r2(sx(v, s.x)), "\" cy=\"", r2(sy(v, s.z)),
    "\" r=\"", r2(r), "\" fill=\"", fill, "\" opacity=\"", r2(opacity), "\"/>")
end

function drawblock(io, v)
  println(io, "<rect x=\"", r2(sx(v, 0) - 5), "\" y=\"", r2(sy(v, 0) - 5),
    "\" width=\"10\" height=\"10\" fill=\"none\" stroke=\"", BLOCK, "\" stroke-width=\"2\"/>")
end

function svgopen(io, w, h)
  println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", w, "\" height=\"", h,
    "\" viewBox=\"0 0 ", w, " ", h, "\">")
  println(io, "<rect width=\"", w, "\" height=\"", h, "\" fill=\"white\"/>")
end

# ---------------------------------------------------------------- figures

function fig_why(path)
  radii, θ = (95.0, 62.0), 0.0
  samples = drillholes()
  inside = filter(s -> normdist(s.x, s.z, radii, θ) ≤ 1, samples)
  outside = filter(s -> normdist(s.x, s.z, radii, θ) > 1, samples)
  dists = [normdist(s.x, s.z, radii, θ) for s in inside]
  nsec = 8
  secof = i -> azimuthal(inside[i].x, inside[i].z, radii, θ, nsec)

  plain, _, holesA = select(inside, dists, i -> 1; maxtotal=12, maxpersector=99, maxperhole=99)
  sect, secB, holesB = select(inside, dists, secof; maxtotal=12, maxpersector=2, maxperhole=3)

  W, H = 792, 326
  io = IOBuffer()
  svgopen(io, W, H)

  panels = (
    (plain, "Nearest 12 samples", "all twelve come from the one hole that passes through the block", false),
    (sect, "8 sectors · max 2 per sector · max 3 per hole", "same data, same ellipsoid: four holes, spread over six sectors", true)
  )
  for (k, (taken, ttl, sub, rays)) in enumerate(panels)
    v = View(24.0 + (k - 1) * 386, 66.0, 1.35, 130.0, 78.0)
    panel(io, v, 372, 250, ttl, sub; id="why$k")
    rays && for j in 0:(nsec - 1)
      println(io, sectorray(v, radii, θ, 2π * j / nsec))
    end
    println(io, "<path d=\"", ellipsepath(v, radii, θ), "\" fill=\"none\" stroke=\"", ELLIPSE, "\" stroke-width=\"1.6\"/>")
    drawholes(io, v, samples)
    for s in outside
      dot(io, v, s, 2.0, SKIP; opacity=0.45)
    end
    takenset = Set(taken)
    for (i, s) in enumerate(inside)
      i in takenset ? dot(io, v, s, 4.2, TAKEN) : dot(io, v, s, 2.6, SKIP)
    end
    drawblock(io, v)
    endpanel(io)
  end

  y = H - 14
  println(io, "<circle cx=\"34\" cy=\"", y - 4, "\" r=\"4.2\" fill=\"", TAKEN, "\"/>")
  println(io, text(44, y, "selected"; size=11, fill=MUTED))
  println(io, "<circle cx=\"122\" cy=\"", y - 4, "\" r=\"2.6\" fill=\"", SKIP, "\"/>")
  println(io, text(132, y, "available, not used"; size=11, fill=MUTED))
  println(io, "<rect x=\"268\" y=\"", y - 9, "\" width=\"10\" height=\"10\" fill=\"none\" stroke=\"", BLOCK, "\" stroke-width=\"2\"/>")
  println(io, text(284, y, "block being estimated"; size=11, fill=MUTED))
  println(io, text(452, y, "section view · steeply dipping holes · composites every 4 m"; size=11, fill=MUTED))
  println(io, "</svg>")
  write(path, String(take!(io)))

  (plain=holesA, sector=holesB, persector=secB)
end

function fig_sectors(path)
  radii, θ = (92.0, 46.0), deg2rad(-22)
  inside = filter(s -> normdist(s.x, s.z, radii, θ) ≤ 1, drillholes())

  W, H = 792, 300
  io = IOBuffer()
  svgopen(io, W, H)

  specs = (
    ("Octants()", "sign of the rotated coordinates — 4 in 2D, 8 in 3D", 4, true),
    ("Azimuthal(12)", "any number of wedges around the principal axis", 12, false)
  )
  for (k, (ttl, sub, n, isquad)) in enumerate(specs)
    v = View(44.0 + (k - 1) * 386, 66.0, 1.30, 120.0, 78.0)
    panel(io, v, 372, 222, ttl, sub; id="sec$k")
    for j in 0:(n - 1)
      println(io, sectorray(v, radii, θ, 2π * j / n))
    end
    println(io, "<path d=\"", ellipsepath(v, radii, θ), "\" fill=\"none\" stroke=\"", ELLIPSE, "\" stroke-width=\"1.6\"/>")
    for s in inside
      sec = isquad ? quadrant(s.x, s.z, radii, θ) : azimuthal(s.x, s.z, radii, θ, n)
      hue = round(Int, 360 * (sec - 1) / n)
      dot(io, v, s, 3.6, string("hsl(", hue, " 46% 46%)"))
    end
    drawblock(io, v)
    endpanel(io)
  end
  println(io, text(28, H - 14,
    "colour = sector id · ellipsoid rotated 22° · sectors are measured in the rotated, radius-normalised frame, so anisotropy and sectors agree";
    size=11, fill=MUTED))
  println(io, "</svg>")
  write(path, String(take!(io)))
end

function fig_multipass(path)
  θ = deg2rad(-15)
  passes = ((60.0, 26.0), (110.0, 48.0), (180.0, 78.0))
  labels = ("pass 1 · min 8 samples", "pass 2 · min 5 samples", "pass 3 · min 2 samples")
  samples = drillholes()

  W, H = 792, 340
  io = IOBuffer()
  svgopen(io, W, H)
  v = View(33.0, 66.0, 1.20, 304.2, 96.0)
  panel(io, v, 758, 262, "MultiPass",
    "each pass is a complete neighbourhood; the first whose constraints are satisfied wins, and the output records which";
    id="mp")

  drawholes(io, v, samples)
  for s in samples
    dot(io, v, s, 2.2, SKIP)
  end
  passopacity(i) = 1.0 - 0.22 * (i - 1)
  passdash(i) = i == 1 ? "none" : "6 4"
  for (i, radii) in enumerate(passes)
    println(io, "<path d=\"", ellipsepath(v, radii, θ), "\" fill=\"none\" stroke=\"", ELLIPSE,
      "\" stroke-width=\"1.6\" opacity=\"", r2(passopacity(i)), "\" stroke-dasharray=\"", passdash(i), "\"/>")
  end
  drawblock(io, v)
  endpanel(io)

  # legend below the panel — labels on the ellipses themselves collide
  for (i, lbl) in enumerate(labels)
    lx = 28.0 + 232 * (i - 1)
    println(io, "<line x1=\"", r2(lx), "\" y1=\"", H - 44, "\" x2=\"", r2(lx + 30), "\" y2=\"", H - 44,
      "\" stroke=\"", ELLIPSE, "\" stroke-width=\"1.6\" opacity=\"", r2(passopacity(i)),
      "\" stroke-dasharray=\"", passdash(i), "\"/>")
    println(io, text(lx + 38, H - 40, lbl; size=11, fill=INK))
  end
  println(io, text(28, H - 14,
    "a block that pass 1 cannot inform falls through to pass 2, then 3 — and is left missing if no pass qualifies";
    size=11, fill=MUTED))
  println(io, "</svg>")
  write(path, String(take!(io)))
end

# ---------------------------------------------------------------- main

here = @__DIR__
stats = fig_why(joinpath(here, "why-sector-search.svg"))
fig_sectors(joinpath(here, "sectors.svg"))
fig_multipass(joinpath(here, "multipass.svg"))

println("samples per hole, nearest-12 : ", sort(collect(stats.plain)))
println("samples per hole, sectored   : ", sort(collect(stats.sector)))
println("samples per sector, sectored : ", sort(collect(stats.persector)))
println("wrote 3 figures to ", here)
