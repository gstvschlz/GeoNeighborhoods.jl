# Examples

Each example is a Julia script rendered to markdown by
[`generate.jl`](generate.jl). The code shown in a document is the code that
produced the output beneath it, so prose, code, and results cannot drift apart.

**Table 1.** The five examples, in the order in which they build on one
another.

| Example | What it demonstrates |
|---|---|
| [Why sector search matters](01-why-sectors.md) | A nearest-12 search draws all twelve samples from one hole. The same ellipsoid under octant and per-hole quotas spreads them over five holes and moves the kriged estimate by 1.0 %. |
| [Sector schemes](02-sector-schemes.md) | `NoSectors`, `Octants`, and `Azimuthal` compared on one dataset, why `Azimuthal` cannot separate samples above the block from those below, and how rotation carries the sectors with it. |
| [Category rules](03-category-rules.md) | Caps per drillhole, minimum distinct holes, hard domain matching, and explicit per-value quotas, including why `mindistinct` achieves little on its own. |
| [Multipass searches](04-multipass.md) | Three neighborhoods in order, the coverage each one reaches alone, and why 208 of 676 blocks remain unestimated. |
| [Working with the existing kriging](05-kriging.md) | Agreement with `GeoStatsModels.fitpredict` to 9 × 10⁻¹⁶, what changes once the neighborhood does something, and the per-block audit trail. |

## Running them

```bash
mise run examples
```

The task resolves the example environment against the working copy of the
package and then rewrites every document from its script.

## The dataset

Every example shares one synthetic dataset, [`drillholes.jl`](drillholes.jl): a
5 × 5 grid of vertical holes on 45 m centers, composited every 4 m. The test
suite is built on the same definition, so the documented numbers and the tested
behavior describe one dataset.
