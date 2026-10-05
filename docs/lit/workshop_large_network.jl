# Copyright 2026 SOPTIM AG                                                    #src
#                                                                             #src
# Licensed under the Apache License, Version 2.0 (the "License");             #src
# you may not use this file except in compliance with the License.            #src
# You may obtain a copy of the License at                                     #src
#                                                                             #src
#     https://www.apache.org/licenses/LICENSE-2.0                             #src
#                                                                             #src
# Unless required by applicable law or agreed to in writing, software         #src
# distributed under the License is distributed on an "AS IS" BASIS,           #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #src
# See the License for the specific language governing permissions and         #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_large_network.jl                                    #src
# purpose: Literate.jl source of the large-network (PEGASE) notebook.         #src
#          Regenerate with `julia --project=docs docs/generate_notebooks.jl`. #src

# # A large network: PEGASE 2869 from a MATPOWER file
#
# > **Level: Advanced.**
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_large_network.ipynb)
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# The PEGASE cases are fictitious but realistic European transmission
# networks (380/220 kV) distributed with MATPOWER under CC BY 4.0
# ([Josz et al., 2016](https://arxiv.org/abs/1603.01533)). `case2869pegase`
# has 2869 buses, 509 generators, 4582 branches, hundreds of transformers and
# a dozen phase shifters. This notebook reads the MATPOWER file, sorts out its
# transformer conventions, solves it with the sparse direct PV kernel and
# compares the result with the solved state stored in the file.

#nb # ## Setup (Colab)
#nb # This cell installs AnalyticLoadFlow from GitHub (branch `main`) into a
#nb # fresh temporary environment. Run it first, once per session.
#nb using Pkg
#nb Pkg.activate(temp = true)
#nb Pkg.add(url = "https://github.com/SOPTIM/AnalyticLoadFlow.jl", rev = "main")
#nb # For the latest registered release use: Pkg.add("AnalyticLoadFlow")

# ## Warm-up

using AnalyticLoadFlow
using LinearAlgebra
using SparseArrays
using Printf
using Downloads
const A = AnalyticLoadFlow
mismatch(case, res) = maximum(A.compute_demo_mismatch(case, res))
solve_pf_apslf(A.demo_case_9bus(); order = 20)   # warm-up
println("ready")

# ## 1. Getting the case file
#
# The file is fetched from the MATPOWER repository. Inside a checkout of
# AnalyticLoadFlow.jl it goes to `data/_downloaded/` (git-ignored), elsewhere
# to a temporary directory.

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

# ## 2. Transformer conventions
#
# MATPOWER defines the branch `angle` column in degrees and `ratio` as the tap
# on the from side. Converted cases do not always follow this: the PEGASE
# files carry the phase-shift angle in radians with the opposite sign. A
# wrong guess is not a small error, a 0.43 rad shift read as 0.43° is a
# different network. [`matpower_case`](@ref) resolves this from the solved
# state stored in the file: every combination of unit, sign and ratio
# convention is stamped into a Y-bus, and the one that reproduces the stored
# `(Vm, Va)` with the smallest total power mismatch wins.

case = matpower_case(path)
conv = case.conventions
println("angle | sign | ratio     | max mismatch | L1 mismatch of the stored state")
for t in conv.trials
   @printf("%-5s | %+d   | %-9s | %10.2e   | %10.2e %s\n", t.angle_unit, t.angle_sign, t.ratio_convention, t.ref_mismatch_pu, t.ref_mismatch_l1_pu,
      (t.angle_unit, t.angle_sign, t.ratio_convention) == (conv.angle_unit, conv.angle_sign, conv.ratio_convention) ? "<- chosen" : "")
end
shifts = [b.shift_deg for b in case.branches if A.is_phase_shifter(b)]
@printf("\n%d buses, %d branches, %d PV buses, %d transformers, %d phase shifters with angles in [%.1f°, %.1f°]\n",
   size(case.Y, 1), length(case.branches), count(==(:pv), case.bustype), count(A.is_transformer, case.branches), length(shifts), minimum(shifts), maximum(shifts))

# The remaining mismatch of the stored state (about 0.5 pu at the PST buses)
# is a property of the file, not of the import: the stored voltages were not
# produced with exactly this branch model.

# ## 3. Solving
#
# `Y` is sparse, so `solve_pf_apslf` uses the sparse direct PV kernel. The
# default embedding (`germ = :deviation`) keeps the flat germ and ramps bus
# shunts and transformer deviations up with $s$. No Newton polish is needed.

res = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = false, return_coeffs = true)
t = @elapsed res = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = false, return_coeffs = true)
@printf("converged = %s, mode = %s, outer iterations = %d, %.3f s\n", res.converged, res.effective_mode, res.outer_iters, t)
@printf("max mismatch on the physical Y-bus = %.1e pu,  |V| in [%.4f, %.4f] pu\n", mismatch(case, res), minimum(abs.(res.V)), maximum(abs.(res.V)))
dV = abs.(res.V .- case.V_ref)
@printf("distance to the stored state: max %.1e pu (bus %s), mean %.1e pu\n", maximum(dV), case.labels[argmax(dV)], sum(dV) / length(dV))

# ## 4. Why the embedding matters on a large network
#
# Both exact embeddings of theory Section 6.5 give the same solution when the
# series converges at $s = 1$. On a large meshed network they behave very
# differently. The no-load state of PEGASE, with all loads switched off, is
# far from the operating point (Ferranti rise on long lightly loaded lines),
# so the `:noload` path has a Padé pole inside the unit circle and the series
# diverges. The `:deviation` path starts at 1 pu everywhere and converges.

for germ in (:deviation, :noload)
   V, _, Vc, _, _ = A.apslf_pf_pv_direct_sparse(case.Y, case.bustype, case.Pspec, case.Qspec, case.Vm; slack = case.slack, Vslack = ComplexF64(case.Vm[case.slack], 0), order = 40, self_check = false, germ = germ)
   st = A.stability_from_Vcoeff(Vc; slack = case.slack, order = 40)
   @printf("germ = %-10s |V^(0)| max = %.2f   |V^(10)| max = %.1e   |V^(40)| max = %.1e   pole distance %.3f (%s)\n",
      germ, maximum(abs.(Vc[:, 1])), maximum(abs.(Vc[:, 11])), maximum(abs.(Vc[:, 41])), st.dmin, A.st_level(st.dmin))
end

# ## 5. Reactive limits
#
# With `enforce_q_limits = true` the outer loop switches generators that
# leave their band to PQ and re-solves; each outer iteration is one full
# series evaluation.

res_q = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = true)
sw = get(res_q, :switch_log, ())
@printf("converged = %s, outer iterations = %d, PV→PQ switches = %d, max mismatch = %.1e pu\n", res_q.converged, res_q.outer_iters, length(sw), mismatch(case, res_q))
nmax = count(e -> e.side == :max, sw)
@printf("switched at Qmax: %d, at Qmin: %d\n", nmax, length(sw) - nmax)

# ## 6. Branch flows
#
# [`branch_flows`](@ref) works on the imported branch list. The phase
# shifters and the most loaded branches:

flows = branch_flows(res.V, case.branches; baseMVA = case.baseMVA)
println("phase shifters:")
print_branch_flows(filter(f -> f.shift_deg != 0.0, flows); labels = case.labels)
top = sort(flows; by = f -> -abs(f.Pij_MW))[1:5]
println("\nfive most loaded branches:")
print_branch_flows(top; labels = case.labels)
@printf("\ntotal losses: %.1f MW\n", sum(f.Ploss_MW for f in flows))
