# Changelog

## Version 0.10.3
- Phase shifters: the default `shift_embedding = :auto` ramps the shift angle with `s`. On case6495rte the Padé pole pair moves from |s| = 0.30 to 0.73; without NR polish the run still reports `converged = false` (voltage error 8e-6 pu at order 24, 1e-8 pu at order 32, #23)
- Shifts above 90 degrees (vector groups) keep the physical sign of the coupling, also next to negative resistances
- The result reports the embedding that ran (`shift_embedding`, `shift_skipped`)

## Version 0.10.2
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
