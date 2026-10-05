```@meta
EditURL = "../../lit/workshop_large_network.jl"
```

# A large network: PEGASE 2869 from a MATPOWER file

> **Level: Advanced.**

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_large_network.ipynb)

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

The PEGASE cases are fictitious but realistic European transmission
networks (380/220 kV) distributed with MATPOWER under CC BY 4.0
([Josz et al., 2016](https://arxiv.org/abs/1603.01533)). `case2869pegase`
has 2869 buses, 509 generators, 4582 branches, hundreds of transformers and
a dozen phase shifters. This notebook reads the MATPOWER file, solves it
with the direct PV kernel and looks at the series behind the result.

## Warm-up

````@example workshop_large_network
using AnalyticLoadFlow
using LinearAlgebra
using SparseArrays
using Printf
using Downloads
const A = AnalyticLoadFlow
mismatch(case, res) = maximum(A.compute_demo_mismatch(case, res))
solve_pf_apslf(A.demo_case_9bus(); order = 20)   # warm-up
println("ready")
````

## 1. Getting the case file

The file is fetched from the MATPOWER repository. Inside a checkout of
AnalyticLoadFlow.jl it goes to `data/_downloaded/` (git-ignored), elsewhere
to a temporary directory.

````@example workshop_large_network
function fetch_case(name)
   dir = isdir("data") ? joinpath("data", "_downloaded") : mktempdir()
   mkpath(dir)
   path = joinpath(dir, name * ".m")
   isfile(path) || Downloads.download("https://raw.githubusercontent.com/MATPOWER/matpower/master/data/$(name).m", path)
   return path
end
path = fetch_case("case2869pegase")
mp = parse_matpower_m(path)
@printf("%s: baseMVA = %.0f, %d bus rows, %d generator rows, %d branch rows\n", mp.name, mp.baseMVA, size(mp.bus, 1), size(mp.gen, 1), size(mp.branch, 1))
````

## 2. Import

[`matpower_case`](@ref) reads the file with the MATPOWER conventions:
branch `angle` in degrees, `ratio` as the tap on the from side. The
network at a glance:

````@example workshop_large_network
case = matpower_case(path)
nbus = size(case.Y, 1)
shifts = [b.shift_deg for b in case.branches if A.is_phase_shifter(b)]
@printf("%d buses, %d branches, %d PV buses\n", nbus, length(case.branches), count(==(:pv), case.bustype))
@printf("%d transformers, %d of them phase shifters (angles from %.2f° to %.2f°)\n",
   count(A.is_transformer, case.branches), length(shifts), minimum(shifts), maximum(shifts))
````

The file also stores a state `(Vm, Va)` in its bus block. It is not a
solution of this network model; its power mismatch on the imported Y-bus:

````@example workshop_large_network
@printf("stored state: max power mismatch %.1f pu\n", case.conventions.ref_mismatch_pu)
````

The solver does not use the stored state.

## 3. Solving

`Y` is sparse, and so is the whole solver: no conversion is needed. The
default embedding (`germ = :deviation`) keeps the flat germ and ramps bus
shunts and transformer deviations up with $s$. The series is evaluated at
$s = 1$ with Padé approximants; no Newton polish is needed.

````@example workshop_large_network
res = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = false, return_coeffs = true)
t = @elapsed res = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = false, return_coeffs = true)
@printf("converged = %s, mode = %s, outer iterations = %d, %.3f s\n", res.converged, res.effective_mode, res.outer_iters, t)
@printf("max mismatch on the physical Y-bus = %.1e pu,  |V| in [%.4f, %.4f] pu\n", mismatch(case, res), minimum(abs.(res.V)), maximum(abs.(res.V)))
````

## 4. Why the embedding matters on a large network

Both exact embeddings of theory Section 6.5 give the same solution when the
series converges at $s = 1$. On a large meshed network they behave very
differently. The no-load state of PEGASE, with all loads switched off, is
far from the operating point (Ferranti rise on long lightly loaded lines),
so the `:noload` path has a singularity close to $s = 0$ and the series
diverges at $s = 1$. The `:deviation` path starts at 1 pu everywhere and
converges. This is why `:deviation` is the default and the supported germ;
`:noload` is experimental.

The coefficients show it directly. `stability_from_Vcoeff` estimates the
radius of convergence by the root test $|V^{(k)}|^{-1/k}$ over the last
orders (above 1 the series converges at $s = 1$, below 1 it does not) and
rates the case RED when it is below 1.

````@example workshop_large_network
for germ in (:deviation, :noload)
   V, _, Vc, _, _ = A.apslf_pf_pv_direct(case.Y, case.bustype, case.Pspec, case.Qspec, case.Vm; slack = case.slack, Vslack = ComplexF64(case.Vm[case.slack], 0), order = 40, self_check = false, germ = germ)
   st = A.stability_from_Vcoeff(Vc; slack = case.slack, order = 40)
   @printf("germ = %-10s |V^(0)| max = %.2f   |V^(10)| max = %.1e   |V^(40)| max = %.1e   radius ≈ %.3f   %s\n",
      germ, maximum(abs.(Vc[:, 1])), maximum(abs.(Vc[:, 11])), maximum(abs.(Vc[:, 41])), st.radius, A.st_level(st))
end
````

## 5. Reactive limits

With `enforce_q_limits = true` the outer loop switches generators that
leave their band to PQ and re-solves; each outer iteration is one full
series evaluation.

````@example workshop_large_network
res_q = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = true)
sw = get(res_q, :switch_log, ())
@printf("converged = %s, outer iterations = %d, PV→PQ switches = %d, max mismatch = %.1e pu\n", res_q.converged, res_q.outer_iters, length(sw), mismatch(case, res_q))
nmax = count(e -> e.side == :max, sw)
@printf("switched at Qmax: %d, at Qmin: %d\n", nmax, length(sw) - nmax)
````

## 6. Branch flows

[`branch_flows`](@ref) works on the imported branch list. The phase
shifters and the most loaded branches:

````@example workshop_large_network
flows = branch_flows(res.V, case.branches; baseMVA = case.baseMVA)
println("phase shifters:")
print_branch_flows(filter(f -> f.shift_deg != 0.0, flows); labels = case.labels)
top = sort(flows; by = f -> -abs(f.Pij_MW))[1:5]
println("\nfive most loaded branches:")
print_branch_flows(top; labels = case.labels)
@printf("\ntotal losses: %.1f MW\n", sum(f.Ploss_MW for f in flows))
````

