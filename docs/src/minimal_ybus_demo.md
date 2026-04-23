# Minimal Demo (Y-Bus Input)

## Overview

The standalone minimal demo shows how to call APSLF with a directly provided Y-bus matrix and plain Julia arrays. It is intended as a compact integration template rather than a full case-import workflow.

The example covers three cases:

* a small 9-bus teaching case, and
* a synthetic 400 V low-voltage radial street-feeder case with PQ loads only, and
* a synthetic 118-bus-sized integration case.

All three cases are built in Julia, use the same NamedTuple data contract, and run without external case files, downloads, MATPOWER import, or additional network packages.

## Data Contract

The example solver wrapper expects a case NamedTuple with these fields:

| Field | Meaning |
| --- | --- |
| `Y` | bus admittance matrix |
| `bustype` | bus type per bus: `:slack`, `:pv`, or `:pq` |
| `Pspec` | active power injection in pu |
| `Qspec` | reactive power injection in pu |
| `Vm` | voltage magnitude setpoints |
| `Qmin` | PV reactive lower limits |
| `Qmax` | PV reactive upper limits |
| `slack` | slack-bus index |

Optional reporting fields are `labels` and `baseMVA`.

## Run the Example

From the project root:

Direct 118-bus integration example:

```bash
julia --project=. examples/synthetic_118_ybus_demo.jl
```

Direct 400 V low-voltage street-feeder example:

```bash
julia --project=. examples/lv_400v_streets_ybus_demo.jl
```

Tiled-grid scaling example:

```bash
julia --project=. examples/tiled_grid_scaling_demo.jl --buses=100 --samples=3
```

Combined minimal demo:

```bash
julia --project=. examples/minimal_ybus_demo.jl --case=9 --inner=pq
julia --project=. examples/minimal_ybus_demo.jl --case=lv400 --inner=pq
julia --project=. examples/minimal_ybus_demo.jl --case=118 --inner=pq
julia --project=. examples/minimal_ybus_demo.jl --case=all --inner=pq
```

`lv_400v_streets_ybus_demo.jl` is the easiest way to run the LV 400 V example directly. `synthetic_118_ybus_demo.jl` is the easiest way to run the synthetic 118-bus integration example directly. `minimal_ybus_demo.jl` remains the reusable integration template and supports multiple cases and PV handling options.

The script runs automatically unless `ENV["APSLF_SUITE_NO_AUTORUN"] = "1"` is set before including it from another Julia file.

## Cases

### 9-Bus Teaching Case

`demo_case_9bus()` builds a small illustrative case directly from branch parameters. Bus 1 is the slack bus, buses 2 and 3 are PV buses with Q-limits, and buses 4 through 9 are PQ buses.

This case prints the full bus-voltage table and is intentionally small enough to show PV/Q-limit behavior in the console output.

The 9-bus data intentionally uses tight Q-limits on PV buses so the outer-loop `--inner=pq` path demonstrates limited Q-limit handling: PV buses are switched to PQ when computed reactive power violates `Qmin`/`Qmax`. This is a compact teaching mechanism, not full industrial voltage-control, generator capability-curve, transformer tap-control, or OLTC logic.

### Synthetic 400 V LV Street-Feeder Case

`demo_case_lv_400v_streets()` builds a synthetic radial low-voltage network with one 400 V transformer/slack bus and three short street feeders. All non-slack buses are PQ loads, so Q-limit PV→PQ switching is not exercised there; there are no PV buses, Q-limit switching demonstrations, transformer tap controls, external case files, or downloads.

The LV 400 V example is intended as a compact educational integration example, not as a real distribution-grid model. Use `examples/lv_400v_streets_ybus_demo.jl` for the direct console path and `examples/minimal_ybus_demo.jl --case=lv400 --inner=pq` when you want to exercise the combined template.

### Parametric Tiled-Grid Scaling Example

`examples/tiled_grid_scaling_demo.jl` uses the existing tiled-grid helpers, including `APSLF.build_tiled_grid_spec`, to build a synthetic one-voltage-level sparse Y-bus network. The user can request the maximum number of buses with `--buses=N`; the actual number may be lower because the helper chooses a rectangular grid with `rows * cols <= requested_buses`.

Run it from the project root with:

```bash
julia --project=. examples/tiled_grid_scaling_demo.jl --buses=100 --samples=3
```

The example prints the requested and actual bus counts, grid dimensions, branch count, case build timing, solve timing min/median/max across measured samples, convergence status, voltage range, and mismatch summary. It is synthetic educational data for observing scaling behavior and timing on generated sparse Y-bus networks. It is not a benchmark suite and not a real grid model.

### Synthetic 118-Bus Integration Case

`demo_case_118bus_synthetic()` builds a synthetic 118-bus-sized case algorithmically. It creates a connected backbone, adds tie-lines for a lightly meshed network, assigns slack/PV/PQ bus types, and fills power specifications and Q-limits programmatically.

The case follows the same data contract as the 9-bus case. For compact output, the demo prints only the first voltage rows plus summary information for the larger case.

The 118-bus case is synthetic and IEEE-118-sized. It is not the official IEEE 118 benchmark data and should not be interpreted as benchmark-equivalent. Use `examples/synthetic_118_ybus_demo.jl` for the direct console path and `examples/minimal_ybus_demo.jl --case=118 --inner=pq` when you want to exercise the combined template.

## Solver Call

The integration wrapper reads the case fields and calls APSLF:

```julia
res = solve_demo_case(case; inner = :pq, order = 40)
```

Internally, this forwards `case.Y`, `case.bustype`, `case.Pspec`, `case.Qspec`, `case.Vm`, `case.Qmin`, `case.Qmax`, and `case.slack` to `APSLF.solve_pf_apslf_with_pv_q_limits(...)`.

## Reusing the Integration Pattern

Applications can construct the same minimal case object and reuse the wrapper and mismatch check:

```julia
case = (
    Y = Y,
    bustype = bustype,
    Pspec = Pspec,
    Qspec = Qspec,
    Vm = Vm,
    Qmin = Qmin,
    Qmax = Qmax,
    slack = slack,
)

res = solve_demo_case(case; inner = :pq)
maxP, maxQ = compute_demo_mismatch(case, res)
```

## Notes on Scope

The example does not load MATPOWER files. The example does not use Sparlectra. The example does not download external datasets.

The synthetic 118-bus case is intended to demonstrate integration-size structure, not benchmark equivalence. Limited Q-limit handling remains PV→PQ switching only; the demos do not add industrial voltage-regulator, tap-control, or OLTC logic.
