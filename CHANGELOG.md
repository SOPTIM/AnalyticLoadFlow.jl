# Changelog

## Version 0.10.0

### Added
- Transformer branches with ratio and phase shift (PST): `pi_branch`, `transformer_branch`, `build_ybus` (dense or sparse), `branch_flows`, `branch_active_power`, `print_branch_flows`, the 9-bus PST case `demo_case_9bus_pst` and `solve_pf_pst_regulated`, an outer secant loop on the angle of a regulated phase shifter (theory Section 6.5). A fixed shift is exact with the default `germ = :deviation`.
- MATPOWER case reader (`parse_matpower_m`, `matpower_case`) that detects the angle unit, the angle sign and the ratio convention from the stored solution; PEGASE 2869 example.
- Notebooks "Transformers and PST" and "Large network (PEGASE)".

### Changed
- The solver is sparse only. `Y` may still be passed dense or sparse; it is converted to `SparseMatrixCSC` once on entry. The dense direct PV kernel, the dense Jacobian and the dense-direct retry in `solve_pf_apslf` are removed. Existing calls keep working: `apslf_pf_pv_direct_sparse` and `inner = :direct_pv_sparse` are the same kernel as `apslf_pf_pv_direct` and `:direct_pv`; `use_sparse`, `sparse_nbus_min` and `dense_fallback_nbus_max` are accepted without effect.
- Precompile and first call are much shorter: the solver compiles only the path a solve takes (inner kernel, germ variant, debug and report output) instead of every branch, and the Padé and Newton steps call `lu` instead of the `\` polyalgorithm. Results are unchanged on the test and benchmark cases.

### Fixed
- After a non-finite PQ inner solve, a successful direct-PV fallback is now used; before, its result was discarded and the original error raised.

## Version 0.9.16

### Changed
- The precompile workload is selectable via `ANALYTICLOADFLOW_PRECOMPILE_WORKLOAD` (`off`, `core` = default, `full`); the default compiles only the default `solve_pf_apslf` call, which cuts the install-time precompile to a fraction.

## Version 0.9.15

### Added
- Precompile workload (PrecompileTools): the solver paths used by the documentation and the workshop notebook are compiled at install time, so first calls no longer pause for compilation.

### Changed
- Requires Julia ≥ 1.12 and is tested on 1.12 and 1.13 (Google Colab currently provides 1.12); the test runner handles `Test.TESTSET_PRINT_ENABLE` as a `ScopedValue` on 1.13.

### Fixed
- Fixed the PQ recursion to use the reflected reciprocal `conj(W^(n-1))` on the right-hand side (theory 1.7, Section 2.4).
- Fixed the sign of the reactive-power unknown in the direct PV kernels; the PV active power is now met without NR polish.
- Fixed the order-0 state: with the new default `germ = :deviation` (or `:noload`) line shunts and `Vslack ≠ 1` are exact, so pure APSLF is a load-flow solution without NR polish.
- Added `pv_secant_damping` (default 1.0) for the outer PV loop.

## Version 0.9.14

### New Features

- Added self-contained demo helpers and console examples for direct Y-bus input.
- Added a parametric synthetic tiled-grid scaling example with configurable requested bus count and compact timing output.

### Fixed
- Fixed misleading CLI wording for synthetic integration cases.
- Fixed Documenter warnings for missing API docstrings and removed the public docs link to the repository-level patent note.
