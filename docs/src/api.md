# API

## Public Entrypoints

```@docs
APSLF.solve_pf_apslf
APSLF.solve_pf_apslf_with_pv_q_limits
APSLF.apslf_pf_pv_direct
APSLF.apslf_pq
```

## Series Evaluation and Stability Helpers

```@docs
APSLF.APSLFEvaluationOptions
APSLF.evaluate_series
APSLF.is_taylor_good
APSLF.pade_eval
APSLF.poly_roots
APSLF.pade_poles
APSLF.stability_from_Vcoeff
APSLF.st_level
```

## Solver Internals Exposed for Testing and Experimentation

```@docs
APSLF.APSLFPQWorkspace
APSLF.build_apslf_pq_workspace
APSLF.NRRectCache
APSLF.build_nr_rect_cache
APSLF.maybe_sparse_Y
APSLF.calc_injections
APSLF.calc_injections!
APSLF.mismatch_rectangular!
APSLF.build_rect_jac_sparse
APSLF.build_rect_jac_dense
APSLF.nr_refine_step_rect!
```

## Validation and Formatting Utilities

```@docs
APSLF.max_mismatch
APSLF.max_mismatch_on_specY
APSLF.any_pv_at_qlimit
APSLF.print_bus_voltages
APSLF.print_line_flows
APSLF.with_silent
APSLF.safe_get
```

## Post-processing

```@docs
APSLF.total_line_losses
APSLF.line_flows_pi
```
