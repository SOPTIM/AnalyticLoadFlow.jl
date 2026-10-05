# Changelog

## Version 0.10.0

### Added
- Transformers with ratio and phase shift (`transformer_branch`, `build_ybus`, `branch_flows`) and a regulated PST (`solve_pf_pst_regulated`).
- MATPOWER reader (`matpower_case`); detects angle unit, angle sign and ratio convention.
- Notebooks "Transformers and PST" and "Large network (PEGASE)".
- Notebook test `docs/test_notebooks.jl`, runs in CI.

### Changed
- Supported: direct mode, Padé, `germ = :deviation`. Outer mode, Taylor evaluation, `:noload` and `:flat` are experimental.
- `solve_pf_apslf` evaluates with Padé only; `use_pade` has no effect.
- Default `inner = :direct_pv` in `solve_pf_apslf_with_pv_q_limits` and `solve_demo_case`.
- Sparse only; a dense `Y` is converted on entry. `use_sparse`, `sparse_nbus_min` and `:direct_pv_sparse` are still accepted.
- `stability_from_Vcoeff`: no spurious poles, new fields `radius` and `level`.
- Precompile with the default workload: 3.7 s instead of 6.3 s (Julia 1.13).

### Fixed
- Outer mode: `converged` now checks the P/Q mismatch.
- Direct-PV fallback after a non-finite PQ solve: result is used, not discarded.
- Workshop notebooks: outputs match the text again.
- Notebook generator on a fresh CI depot.

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
