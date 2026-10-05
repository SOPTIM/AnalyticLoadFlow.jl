# Changelog

## Version 0.10.0

### Added
- Transformers and phase shifters, regulated PST.
- MATPOWER reader.
- Notebooks "Transformers and PST" and "Large network (PEGASE)".
- Notebook test in CI.

### Changed
- Sparse only.
- Padé only.
- Direct mode is the default.
- Outer mode, Taylor evaluation, `:noload` and `:flat` are experimental.
- `stability_from_Vcoeff`: new fields `radius` and `level`.
- Faster precompile.

### Fixed
- Outer mode: `converged` checks the mismatch.
- Direct-PV fallback result is used.

## Version 0.9.16

### Changed
- Selectable precompile workload (`ANALYTICLOADFLOW_PRECOMPILE_WORKLOAD`).

## Version 0.9.15

### Added
- Precompile workload.
- `pv_secant_damping` for the outer PV loop.

### Changed
- Julia 1.12 and 1.13.

### Fixed
- PQ recursion uses `conj(W^(n-1))`.
- Sign of the reactive-power unknown in the direct PV kernels.
- Exact order-0 state with `germ = :deviation`.

## Version 0.9.14

### Added
- Demo helpers and console examples for Y-bus input.
- Tiled-grid scaling example.

### Fixed
- CLI wording for synthetic cases.
- Documenter warnings.
