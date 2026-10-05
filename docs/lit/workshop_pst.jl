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
# file: docs/lit/workshop_pst.jl                                              #src
# purpose: Literate.jl source of the transformer / phase-shifter notebook     #src
#          (theory Section 6.5). Regenerate with                              #src
#          `julia --project=docs docs/generate_notebooks.jl`.                 #src

# # Transformers and phase shifters in APSLF
#
# > **Level: Advanced**, companion of theory Section 6.5.
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_pst.ipynb)
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# A transformer enters the Y-bus through a complex tap $t = a\,e^{j\varphi}$.
# The ratio $a$ moves voltage, the phase shift $\varphi$ moves active power.
# For the series recursion of APSLF a phase-shifting transformer (PST) is
# special: it drives a circulating current even when all bus voltages are
# equal, so the flat germ of the basic method is no longer an exact
# order-0 solution, and the Y-bus is no longer symmetric. This notebook
# builds the branch model, shows the problem on the row sums, works the two
# consistent embeddings of theory Section 6.5 by hand on a 4-bus network,
# and then uses the built-in kernels for an angle sweep and a regulated PST.

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
const A = AnalyticLoadFlow

mismatch(case, res) = maximum(A.compute_demo_mismatch(case, res))
pol(z) = @sprintf("%.4f ∠ %7.3f°", abs(z), rad2deg(angle(z)))
solve_pf_apslf(demo_case_9bus_pst(shift_deg = 3.0); order = 20)   # warm-up
println("ready")

# ## 1. The branch model
#
# With series admittance $y = 1/(r + jx)$ and tap $t$ on side $i$ the branch
# contributes (theory Section 6.5)
#
# ```math
# I_i = \frac{y}{|t|^2} V_i - \frac{y}{\bar t} V_k, \qquad
# I_k = -\frac{y}{t} V_i + y V_k .
# ```
#
# [`transformer_branch`](@ref) stores the parameters, [`branch_admittances`](@ref AnalyticLoadFlow.branch_admittances)
# returns the four entries. For a pure phase shifter $|t| = 1$: the diagonal
# entries are unchanged, the off-diagonal entries pick up $e^{\pm j\varphi}$.

y = inv(complex(0.01, 0.10))
for (ratio, φ) in ((1.0, 0.0), (0.95, 0.0), (1.0, 10.0))
   br = transformer_branch(1, 2; r = 0.01, x = 0.10, ratio = ratio, shift_deg = φ)
   yii, yij, yji, yjj = A.branch_admittances(br)
   @printf("ratio %.2f  shift %5.1f°:  Y_ii = %s   Y_ij = %s   Y_ji = %s\n", ratio, φ, pol(yii), pol(yij), pol(yji))
end

# ## 2. What the PST does to the germ
#
# The flat germ rests on zero row sums: with all voltages equal, $Y\mathbf{1}$
# is the current vector, and it must vanish. A transformer at nominal tap
# has zero row sums, a PST does not. `build_ybus(...; nominal = true)`
# returns the matrix with every tap forced to $1\angle 0°$ and shunts
# removed, which is the $Y_0$ of Section 6.5.

pst_case = demo_case_9bus_pst(shift_deg = 10.0)     # PST lumped into branch 4-5, tap on bus 4
Y = pst_case.Y
Y0 = pst_case.Y0
I0 = Y * ones(9)                                    # current at the flat state
@printf("flat state, full Y : |I| at bus 4 = %.4f, bus 5 = %.4f, bus 6 = %.4f pu\n", abs(I0[4]), abs(I0[5]), abs(I0[6]))
@printf("flat state, Y0     : max |I| = %.1e pu   (zero row sums)\n", maximum(abs.(Y0 * ones(9))))
@assert maximum(abs.(Y0 * ones(9))) < 1e-12 && abs(I0[4]) > 1   #src
@printf("symmetry: Y[4,5] = %s,  Y[5,4] = %s\n", pol(Y[4, 5]), pol(Y[5, 4]))

# ## 3. The two embeddings by hand
#
# A small PQ-only network makes the recursion transparent: four buses, a PST
# on branch 1–3, one line with charging. Bus 1 is the slack at 1 pu.

br4 = [
   pi_branch(1, 2; r = 0.02, x = 0.10, b = 0.04),
   transformer_branch(1, 3; r = 0.01, x = 0.08, shift_deg = 8.0),
   pi_branch(2, 3; r = 0.03, x = 0.15),
   pi_branch(3, 4; r = 0.02, x = 0.12),
]
Y4 = build_ybus(4, br4)
Y4_0 = build_ybus(4, br4; nominal = true)
S4 = ComplexF64[0, -0.4 - 0.15im, -0.3 - 0.10im, -0.25 - 0.05im]
red = 2:4
Vslack = 1.0 + 0.0im

## shared pieces: inverse series W = 1/V and Padé evaluation at s = 1
function inverse_series!(W, V, n)
   for i in axes(V, 1)
      acc = zero(eltype(V))
      for m = 1:n
         acc += V[i, m+1] * W[i, n-m+1]
      end
      W[i, n+1] = -acc / V[i, 1]
   end
end
evaluate(Vc) = [A.evaluate_series(Vc[i, :], A.APSLFEvaluationOptions(mode = :pade)).voltage for i in axes(Vc, 1)]
order = 30;

# **Variant 1, deviation embedding.** Constant matrix $Y_0$ (zero row sums),
# flat germ, and the deviation $\Delta Y = Y - Y_0$ moves to the right-hand
# side, scaled with $s$ exactly like the injections:
#
# ```math
# Y_{0,\mathrm{red}} V^{(n)} = S^* \odot \overline{W^{(n-1)}} - \bigl(\Delta Y\, V^{(n-1)}\bigr)_{\mathrm{red}}, \qquad V^{(0)} = \mathbf{1}.
# ```
#
# $V^{(n-1)}$ is the full vector including the slack entry ($V_1^{(0)} = V_{\mathrm{slack}}$, zero for $n \ge 1$).

ΔY = Y4 - Y4_0
F0 = lu(Y4_0[red, red])
V1 = zeros(ComplexF64, 4, order + 1); W1 = zeros(ComplexF64, 3, order + 1)
V1[:, 1] .= 1.0; V1[1, 1] = Vslack; W1[:, 1] .= 1.0
for n = 1:order
   rhs = conj.(S4[red]) .* conj.(W1[:, n]) .- (ΔY * V1[:, n])[red]
   V1[red, n+1] = F0 \ rhs
   inverse_series!(W1, V1[red, :], n)
end
Vdev = [Vslack; evaluate(V1[red, :])]
Sdev = Vdev .* conj.(Y4 * Vdev)
@printf("variant 1: max |S - S_spec| at buses 2-4 = %.2e pu\n", maximum(abs.(Sdev[red] .- S4[red])))
@assert maximum(abs.(Sdev[red] .- S4[red])) < 1e-12   #src

# **Variant 2, no-load germ.** Constant matrix is the full $Y$; the germ is the
# solution of the linear no-load problem, and the recursion of Section 4 runs
# unchanged with a non-uniform $W^{(0)} = 1/V^{(0)}$:
#
# ```math
# Y_{\mathrm{red}} V^{(0)} = -Y_{\mathrm{red},1} V_1, \qquad
# Y_{\mathrm{red}} V^{(n)} = S^* \odot \overline{W^{(n-1)}} .
# ```

F = lu(Y4[red, red])
V2 = zeros(ComplexF64, 3, order + 1); W2 = zeros(ComplexF64, 3, order + 1)
V2[:, 1] = F \ (-Y4[red, 1] .* Vslack)
W2[:, 1] .= 1.0 ./ V2[:, 1]
for n = 1:order
   V2[:, n+1] = F \ (conj.(S4[red]) .* conj.(W2[:, n]))
   inverse_series!(W2, V2, n)
end
Vnl = [Vslack; evaluate(V2)]
Snl = Vnl .* conj.(Y4 * Vnl)
@printf("variant 2: max |S - S_spec| at buses 2-4 = %.2e pu\n", maximum(abs.(Snl[red] .- S4[red])))
@assert maximum(abs.(Snl[red] .- S4[red])) < 1e-12   #src
println("no-load germ: ", join([pol(v) for v in V2[:, 1]], ",  "))
@printf("both variants, same function at s = 1: max |V_dev - V_noload| = %.1e\n", maximum(abs.(Vdev .- Vnl)))
@assert maximum(abs.(Vdev .- Vnl)) < 1e-12   #src

# The coefficient paths differ although the limits agree. Theory Section 6.5:
# the two embeddings have different convergence radii.

for n in (1, 2, 5, 10, 20, 30)
   @printf("n = %2d   max |V^(n)|  variant 1: %.2e   variant 2: %.2e\n", n, maximum(abs.(V1[red, n+1])), maximum(abs.(V2[:, n+1])))
end

# The package does the same through the `germ` keyword. Its `:deviation`
# variant (the default and the supported one) uses
# $Y_0 = Y - \mathrm{diag}(Y\mathbf{1})$, which has zero row sums for any $Y$
# and makes the deviation a diagonal matrix; the result is the same
# solution. `:noload` (variant 2) is available as an experimental option.

case4 = (Y = Y4, bustype = [:slack, :pq, :pq, :pq], Pspec = real.(S4), Qspec = imag.(S4), Vm = ones(4), Qmin = fill(-1e9, 4), Qmax = fill(1e9, 4), slack = 1)
for germ in (:deviation, :noload)
   r = solve_pf_apslf(case4; order = order, nr_polish = false, germ = germ)
   @printf("solve_pf_apslf germ = %-10s max |V - V_hand| = %.1e\n", germ, maximum(abs.(r.V .- Vdev)))
   @assert r.converged && maximum(abs.(r.V .- Vdev)) < 1e-12   #src
end

# ## 4. Angle sweep on the 9-bus ring
#
# The 9-bus ring 4-5-6-7-8-9-4 offers two paths from each generator to each
# load. A PST on branch 4–5 shifts active power from one path to the other;
# the infeed at bus 4 stays what generator 1 delivers.

println("shift°     P_45 MW     P_49 MW    P_45+P_49   losses MW   |V| min   |V| max")
for φ in -20.0:5.0:20.0
   c = demo_case_9bus_pst(shift_deg = φ, enforce_q_limits = false)
   r = solve_pf_apslf(c; order = 40, nr_polish = false)
   fl = branch_flows(r.V, c.branches; baseMVA = c.baseMVA)
   p45 = branch_active_power(r.V, c.branches, 4, 5) * c.baseMVA
   p49 = branch_active_power(r.V, c.branches, 4, 9) * c.baseMVA
   @printf("%6.1f  %10.3f  %10.3f  %10.3f  %10.4f  %8.4f  %8.4f\n", φ, p45, p49, p45 + p49, sum(f.Ploss_MW for f in fl), minimum(abs.(r.V)), maximum(abs.(r.V)))
   @assert r.converged && mismatch(c, r) < 1e-8   #src
end

# The ratio tap does the other job. Lowering the tap at bus 4 raises the
# voltage on the far side; the active split barely moves:

println("ratio    |V5|      |V9|     P_45 MW")
for a in (0.95, 1.0, 1.05)
   c = demo_case_9bus_pst(ratio = a, enforce_q_limits = false)
   r = solve_pf_apslf(c; order = 40, nr_polish = false)
   @printf("%.2f   %.4f   %.4f   %8.3f\n", a, abs(r.V[5]), abs(r.V[9]), branch_active_power(r.V, c.branches, 4, 5) * c.baseMVA)
   @assert r.converged && mismatch(c, r) < 1e-8   #src
end

# ## 5. A regulated phase shifter
#
# A regulated PST adjusts $\varphi$ until the active power on its branch meets
# a setpoint. The angle enters $Y$ through $e^{j\varphi}$, i.e. not
# polynomially in $s$, so it cannot be expanded inside the recursion. It is
# handled like PV buses and reactive limits: an outer loop, here the secant
# method, changes the angle and re-runs APSLF (theory Section 6.5).
# [`solve_pf_pst_regulated`](@ref) takes a function that builds the case for a
# given angle.

reg = solve_pf_pst_regulated(φ -> demo_case_9bus_pst(shift_deg = φ, enforce_q_limits = false), 4, 5, 0.6; order = 40, nr_polish = false, verbose = 1)
@printf("\nshift = %.4f°  P_45 = %.6f pu  converged = %s  APSLF solves = %d\n", reg.shift_deg, reg.P_pu, reg.converged, length(reg.history))
@assert reg.converged && abs(reg.P_pu - 0.6) < 1e-5   #src

# Angle limits pin the result and report `converged = false`:

pinned = solve_pf_pst_regulated(φ -> demo_case_9bus_pst(shift_deg = φ, enforce_q_limits = false), 4, 5, 0.6; shift_min = -3.0, shift_max = 3.0, order = 40, nr_polish = false)
@printf("with limits ±3°: shift = %.2f°, P_45 = %.4f pu, converged = %s\n", pinned.shift_deg, pinned.P_pu, pinned.converged)
@assert !pinned.converged && pinned.shift_deg == -3.0   #src

# ## 6. Sparse Y with transformers
#
# `build_ybus(...; sparse_output = true)` returns a `SparseMatrixCSC`, the
# form the solver works on (a dense `Y` is converted on entry). The
# non-symmetric pattern of a PST needs no special treatment.

c = demo_case_9bus_pst(shift_deg = 12.0, sparse_output = true)
r = solve_pf_apslf(c; order = 40, nr_polish = false)
@printf("sparse PST case: Y is %s, converged = %s, max mismatch = %.1e pu\n", typeof(c.Y).name.name, r.converged, mismatch(c, r))
@assert r.converged && mismatch(c, r) < 1e-10   #src
