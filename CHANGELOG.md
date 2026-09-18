# Changelog

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
