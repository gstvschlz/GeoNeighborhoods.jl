# Examples

Each example is a Julia script rendered to markdown by
[`generate.jl`](generate.jl). The code in the document is the code that
produced the output beneath it, so the two cannot drift apart.

| | |
|---|---|
| [Why sector search](01-why-sectors.md) | A nearest-12 search takes all twelve samples from one hole. The same ellipsoid with octant and per-hole quotas spreads across five, and moves the kriged estimate. |
| [Sector schemes](02-sector-schemes.md) | `NoSectors`, `Octants` and `Azimuthal` compared on the same data, why `Azimuthal` cannot separate above from below, and how rotation carries the sectors with it. |
| [Category rules](03-category-rules.md) | Caps per drillhole, minimum distinct holes, hard domain matching and explicit per-value quotas — including why `mindistinct` is nearly useless on its own. |
| [Multipass searches](04-multipass.md) | Three neighbourhoods in order, what each covers alone, and why 208 of 676 blocks are still left unestimated. |
| [Working with the existing kriging](05-kriging.md) | Agreement with `GeoStatsModels.fitpredict` to 9e-16, what changes once the neighbourhood does something, and the per-block audit trail. |

## Running them

```bash
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/generate.jl
```

The shared synthetic dataset is [`drillholes.jl`](drillholes.jl) — a 5 × 5 grid
of vertical holes on 45 m centres, composited every 4 m. The test suite uses the
same definition, so the documented numbers and the tested behaviour describe one
dataset.
