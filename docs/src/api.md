# API

## Public Entrypoints

```@docs
AnalyticLoadFlow.solve_pf_apslf
AnalyticLoadFlow.solve_pf_apslf_with_pv_q_limits
AnalyticLoadFlow.apslf_pf_pv_direct
AnalyticLoadFlow.apslf_pq
```

## Series Evaluation and Stability Helpers

```@docs
AnalyticLoadFlow.APSLFEvaluationOptions
AnalyticLoadFlow.evaluate_series
AnalyticLoadFlow.is_taylor_good
AnalyticLoadFlow.pade_eval
AnalyticLoadFlow.poly_roots
AnalyticLoadFlow.pade_poles
AnalyticLoadFlow.stability_from_Vcoeff
AnalyticLoadFlow.st_level
```

## Solver Internals Exposed for Testing and Experimentation

```@docs
AnalyticLoadFlow.APSLFPQWorkspace
AnalyticLoadFlow.build_apslf_pq_workspace
AnalyticLoadFlow.NRRectCache
AnalyticLoadFlow.build_nr_rect_cache
AnalyticLoadFlow.maybe_sparse_Y
AnalyticLoadFlow.calc_injections
AnalyticLoadFlow.calc_injections!
AnalyticLoadFlow.mismatch_rectangular!
AnalyticLoadFlow.build_rect_jac_sparse
AnalyticLoadFlow.build_rect_jac_dense
AnalyticLoadFlow.nr_refine_step_rect!
```

## Validation and Formatting Utilities

```@docs
AnalyticLoadFlow.max_mismatch
AnalyticLoadFlow.max_mismatch_on_specY
AnalyticLoadFlow.any_pv_at_qlimit
AnalyticLoadFlow.print_bus_voltages
AnalyticLoadFlow.print_line_flows
AnalyticLoadFlow.with_silent
AnalyticLoadFlow.safe_get
```

## Germ and Embedding

```@docs
AnalyticLoadFlow.apslf_germ
AnalyticLoadFlow.apslf_noload_germ
AnalyticLoadFlow.apslf_row_sums
```

## Branches, Transformers and Phase Shifters

```@docs
AnalyticLoadFlow.pi_branch
AnalyticLoadFlow.transformer_branch
AnalyticLoadFlow.branch_admittances
AnalyticLoadFlow.build_ybus
AnalyticLoadFlow.branch_flows
AnalyticLoadFlow.branch_active_power
AnalyticLoadFlow.print_branch_flows
AnalyticLoadFlow.demo_9bus_branches
AnalyticLoadFlow.demo_case_9bus_pst
AnalyticLoadFlow.solve_pf_pst_regulated
```

## MATPOWER Import

```@docs
AnalyticLoadFlow.parse_matpower_m
AnalyticLoadFlow.matpower_case
```

## Demo Cases

```@docs
AnalyticLoadFlow.build_ybus_from_branches
AnalyticLoadFlow.demo_case_9bus
AnalyticLoadFlow.solve_demo_case
```

## Post-processing

```@docs
AnalyticLoadFlow.total_line_losses
AnalyticLoadFlow.line_flows_pi
```
