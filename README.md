# AnalyticLoadFlow.jl

Compact Julia reference implementation of an analytic power-series based AC load-flow approach.

## Project Status

AnalyticLoadFlow.jl is provided by SOPTIM AG as a compact reference implementation of an analytic power-series based AC load-flow approach.

The repository is intentionally kept small and mostly static. It is provided for study, reproducibility, and experimentation with the method. SOPTIM AG does not maintain AnalyticLoadFlow.jl as an industrial power-flow product, does not provide commercial support for it, and does not offer it as part of a commercial product or service.

The code is made available under the Apache-2.0 license. This does not imply any warranty, maintenance obligation, support commitment, or product roadmap by SOPTIM AG.

## Feature Matrix

| Feature | Status | Notes |
| --- | --- | --- |
| PQ buses / PQ load flow | Supported | APSLF PQ core |
| PV buses | Supported | outer-loop PV handling and direct PV kernel |
| Q-limit handling | Supported, limited | PV→PQ switching when `Qmin`/`Qmax` are violated |
| NR polish | Optional | rectangular Newton post-processing |
| Padé evaluation | Supported | used for series evaluation / analytical continuation |
| 400 V LV PQ-only example | Included | synthetic educational radial street-feeder case |
| Parametric tiled-grid scaling example | Included | synthetic one-voltage-level sparse Y-bus, configurable with `--buses=N`, includes timing output |
| External case import | Not included | no MATPOWER/CGMES import workflow |
| Transformer tap control / OLTC | Not included | no tap-control or voltage-regulator logic |
| Industrial support/productization | Not included | reference implementation only |

## Changelog

See CHANGELOG.md for notable user-visible changes.

## What is included

- APSLF solver core
- Padé evaluation and stability indicator helpers
- Internal Newton polish where used by the solver core
- Self-contained Y-bus examples
- Theory article in `docs/src/theorie-eng.md`

## What is not included

- Third-party network-framework integration
- External case-file import workflow
- Benchmark or rescue suites
- GUI, Web UI, or API service layer
- Industrial support promise

## Installation

From the Julia package registry:

```julia
using Pkg
Pkg.add("AnalyticLoadFlow")
using AnalyticLoadFlow
```

From a local checkout, run commands with `julia --project=.` from the repository root.

## Run the example

The easiest way to run the synthetic IEEE-118-sized integration case directly is:

```bash
julia --project=. examples/synthetic_118_ybus_demo.jl
```

The direct 400 V low-voltage street-feeder example is:

```bash
julia --project=. examples/lv_400v_streets_ybus_demo.jl
```

The LV 400 V example is a synthetic radial street-feeder case with one transformer/slack bus and PQ loads only. It is intended as a compact educational integration example, not as a real distribution-grid model.

The parametric tiled-grid scaling example can be run with a requested maximum bus count and timing samples:

```bash
julia --project=. examples/tiled_grid_scaling_demo.jl --buses=100 --samples=3
```

The tiled-grid scaling example is synthetic educational data. It is useful for observing scaling behavior and timing on generated sparse Y-bus networks. It is not a benchmark suite and not a real grid model.

The reusable minimal integration template remains available for multiple cases:

```bash
julia --project=. examples/minimal_ybus_demo.jl
julia --project=. examples/minimal_ybus_demo.jl --case=9 --inner=pq
julia --project=. examples/minimal_ybus_demo.jl --case=lv400 --inner=pq
julia --project=. examples/minimal_ybus_demo.jl --case=118 --inner=pq
julia --project=. examples/minimal_ybus_demo.jl --case=all --inner=pq
```

- `--case=9 --inner=pq` runs the small teaching case with outer-loop PV handling and shows limited Q-limit handling for PV buses via PV→PQ switching.
- `--case=lv400` runs the synthetic 400 V low-voltage radial street-feeder case; it is PQ-only and does not demonstrate Q-limit switching.
- `--case=118` runs the synthetic 118-bus-sized integration case; it is generated data, not the official IEEE 118 benchmark.
- `--case=all` runs the 9-bus, LV 400 V, and synthetic 118-bus cases.
- `--inner=pq` uses the outer-loop PV logic with PQ inner solves.

The console examples are thin entry points: reusable demo data and solver helpers live in `src/demo_cases.jl`, while presentation and argument parsing stay in `examples/`. The examples print summaries for users and intentionally return `nothing` from `main`. They provide the same NamedTuple data contract: `Y`, `bustype`, `Pspec`, `Qspec`, `Vm`, `Qmin`, `Qmax`, and `slack`.

The `tiled_grid_scaling_demo.jl` script is the direct entry point for synthetic tiled-grid timing runs. The `lv_400v_streets_ybus_demo.jl` and `synthetic_118_ybus_demo.jl` scripts are the shortest direct entry points for their generated examples. The `minimal_ybus_demo.jl` script remains the console integration template and supports the 9-bus, LV 400 V, synthetic 118-bus, and combined runs. The LV 400 V and 118-bus cases are synthetic generated data; the LV case is not a real grid model, and the 118-bus case is not the official IEEE 118 benchmark data. Tests validate the reusable helpers and solver path directly via `using AnalyticLoadFlow`; they do not include example scripts.

## Architecture

See `docs/src/architecture_overview.md` for the public interface, Y-bus data model, directory structure, solver flow, and feature scope.

## Basic usage

```julia
using AnalyticLoadFlow

spec = (
    Y = Y,
    bustype = bustype,
    Pspec = Pspec,
    Qspec = Qspec,
    Vm = Vm,
    Qmin = Qmin,
    Qmax = Qmax,
    slack = slack_bus,
)

res = solve_pf_apslf(
    spec;
    mode = :direct,
    order = 40,
    use_pade = true,
    nr_polish = true,
)
```

APSLF uses the canonical analytic germ `V(s=0)=1∠0`. This is not a user-selectable Newton-style start value. If Newton polish is enabled, it starts from the APSLF solution.

## Theory

See `docs/src/theorie-eng.md` for the theory article.

## License 

Source code license: Apache-2.0. Refer to `LICENSE` for the full license text.
