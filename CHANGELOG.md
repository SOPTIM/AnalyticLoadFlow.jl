# Changelog

## Unreleased
- Non-converged result carries `reason` and `diagnosis` (nearest Padé pole, order comparison)
- Warnings for a `:noload` germ near zero and a `:flat` germ with an `s = 0` residual

## Version 0.10.1
- Fixed: NR-polish autodamping test depended on round-off (FMA contraction on LLVM 23), see #20

## Version 0.10.0
- Transformers, phase shifters, regulated PST
- MATPOWER reader
- PST and PEGASE notebooks, notebook test
- Sparse only, Padé only, direct mode as default
- Experimental: outer mode, Taylor, `:noload`, `:flat`
- Stability indicator with radius
- Faster precompile
- Fixes: outer-mode convergence, direct-PV fallback

## Version 0.9.16
- Selectable precompile workload

## Version 0.9.15
- Precompile workload
- Julia 1.12 and 1.13
- `pv_secant_damping`
- Fixes: PQ recursion, PV sign, exact germ

## Version 0.9.14
- Y-bus demo examples
- Tiled-grid scaling example
- Fixes: CLI wording, Documenter warnings
