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
# See the License for the specific language governing permissions and        #src
# limitations under the License.                                              #src
#                                                                             #src
# file: docs/lit/workshop_apslf.jl                                            #src
# purpose: Literate.jl source of the single APSLF workshop notebook: the      #src
#          hand calculations of theory Section 7 against the solver, then the #src
#          9-bus case (series, Padé, germ, PV buses, sparse path).            #src
#          Tables are written twice: an HTML table with a border for the      #src
#          notebook (`#nb`, Colab draws Markdown tables without lines) and a  #src
#          Markdown table for the Documenter page (`#md`). Links into the     #src
#          theory article open in a new tab in the notebook (`#nb`, HTML      #src
#          anchor) and are plain Markdown links on the Documenter page.       #src
#          Regenerate with `julia --project=docs docs/generate_notebooks.jl`. #src

# # APSLF workshop: the theory article against the solver.
#
# > **Level: Newcomer to Advanced.** One notebook, one install; runs in about a minute afterwards.
#
# [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_apslf.ipynb)
#
# > **Note:** This workshop was created with AI assistance and is reviewed
# > and curated by the maintainer; it is not a fully machine-generated text.
#
# This notebook follows the
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md" target="_blank">theory article</a>
#md # [theory article](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md)
# (version 1.7). Every section reference below links to the article; in the
# notebook the links open in a new tab, so you can read the derivation next
# to the code.
#
# **Part 1** takes the two hand calculations of Section 7, the two-bus network
# with a capacitor bank (7.10) and the four-bus π-model network (7.1 to 7.8),
# and checks the solver coefficients digit by digit against the numbers
# printed in the article. The recursion uses the reflected reciprocal
# $\overline{W^{(n-1)}}$ on the right-hand side
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#24-treatment-of-complex-conjugation" target="_blank">Section 2.4</a>),
#md # ([Section 2.4](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#24-treatment-of-complex-conjugation)),
# in the article and in the solver.
#
# **Part 2** uses the 9-bus teaching case that ships with the package to show
# the series, the Padé evaluation, the choice of the germ, PV buses with
# reactive limits and the sparse path.
#
# > **Note:** On Google Colab the install cell takes a few minutes on a
# > fresh session. This notebook targets Julia ≥ 1.12.

#nb # ## Setup (Colab)
#nb # This cell installs AnalyticLoadFlow from GitHub (branch `main`) into a
#nb # fresh temporary environment. Run it first, once per session. To test a
#nb # branch, change `rev`. For a private checkout use a personal access token:
#nb # `Pkg.add(url = "https://USER:TOKEN@github.com/USER/AnalyticLoadFlow.jl", rev = "main")`.
#nb #
#nb # > **If the first line below reports a syntax error**, the notebook was
#nb # > opened with a Python runtime. In Colab choose *Runtime → Change runtime
#nb # > type → Julia* and run the cell again.
#nb println("Julia ", VERSION)
#nb VERSION >= v"1.12" || @warn "AnalyticLoadFlow requires Julia ≥ 1.12; this runtime has $(VERSION). The install below will fail; please report the version shown above."
#nb using Pkg
#nb Pkg.activate(temp = true)
#nb Pkg.add(url = "https://github.com/SOPTIM/AnalyticLoadFlow.jl", rev = "main")
#nb # For the latest registered release use: Pkg.add("AnalyticLoadFlow")

# ## Warm-up and helpers
#
# Three small helpers are used throughout:
#
# - `c(z)` prints a complex number with six decimals, the precision of the
#   article.
# - `show_matrix(name, M)` prints a complex matrix as an aligned table with
#   bus labels.
# - `hand_recursion` is the recursion of the article written down as
#   plainly as possible, **independent of the solver**. It is the reference
#   the solver is checked against.
#
# `hand_recursion` implements
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#73-reduced-system-for-the-non-slack-buses" target="_blank">Section 7.3</a>
#md # [Section 7.3](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#73-reduced-system-for-the-non-slack-buses)
# on the reduced system (slack bus removed). Per order $n \ge 1$ it solves
#
# $$Y^{\mathrm{ser}}_{\mathrm{red}}\, V^{(n)} = S^* \odot \overline{W^{(n-1)}} - y^{\mathrm{sh}}_{\mathrm{red}} \odot V^{(n-1)},
# \qquad V^{(0)} = 1,\; W^{(0)} = 1,$$
#
# and then the coefficient $W^{(n)}$ of the reciprocal series $W(s) = 1/V(s)$
# from the convolution of
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#41-constraint-convolution" target="_blank">Section 4.1</a>.
#md # [Section 4.1](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#41-constraint-convolution).
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>Argument</th><th>Meaning</th></tr>
#nb # <tr><td><code>Yser_red</code></td><td>series matrix of the non-slack buses (zero row sums, constant in s)</td></tr>
#nb # <tr><td><code>ysh_red</code></td><td>shunt vector of the non-slack buses (diagonal of the shunt matrix, ramped with s)</td></tr>
#nb # <tr><td><code>Sred</code></td><td>complex injections P + jQ of the non-slack buses, load negative</td></tr>
#nb # <tr><td><code>order</code></td><td>highest order n to compute</td></tr>
#nb # </table>
#md # | Argument | Meaning |
#md # |:--|:--|
#md # | `Yser_red` | series matrix of the non-slack buses (zero row sums, constant in $s$) |
#md # | `ysh_red` | shunt vector of the non-slack buses (diagonal of the shunt matrix, ramped with $s$) |
#md # | `Sred` | complex injections $P + jQ$ of the non-slack buses, load negative |
#md # | `order` | highest order $n$ to compute |
#
# It returns two matrices `V` and `W` with one row per non-slack bus and
# one column per order, column 1 being order 0.

using AnalyticLoadFlow
using LinearAlgebra
using Printf
const A = AnalyticLoadFlow

## Complex number with six decimals, e.g. "+0.964029 -0.119926im"
c(z) = @sprintf("%+.6f %+.6fim", real(z), imag(z))

## Complex matrix as an aligned table with bus labels
function show_matrix(name, M)
   println(name)
   @printf("%-6s", "")
   for j in axes(M, 2)
      @printf("%-24s", "Bus $j")
   end
   println()
   for i in axes(M, 1)
      @printf("%-6s", "Bus $i")
      for j in axes(M, 2)
         @printf("%-24s", c(M[i, j]))
      end
      println()
   end
end

## The recursion of theory Section 7.3 on the reduced system:
##   Yser_red V^(n) = S* ⊙ conj(W^(n-1)) - ysh_red ⊙ V^(n-1),   V^(0) = 1, W^(0) = 1
## and W^(n) from the convolution of Section 4.1.
## Returns V, W: one row per non-slack bus, columns are the orders 0..order.
function hand_recursion(Yser_red, ysh_red, Sred, order)
   n = length(Sred)
   V = zeros(ComplexF64, n, order + 1)
   W = zeros(ComplexF64, n, order + 1)
   V[:, 1] .= 1                              # order 0: flat germ
   W[:, 1] .= 1                              # W^(0) = 1 / V^(0)
   F = lu(Yser_red)                          # factorize once, reuse for every order
   for k = 1:order
      rhs = conj.(Sred) .* conj.(W[:, k]) .- ysh_red .* V[:, k]
      V[:, k+1] = F \ rhs                     # network equation of order k
      for i = 1:n                             # convolution: W^(k) from V^(1..k) and W^(0..k-1)
         W[i, k+1] = -sum(V[i, m+1] * W[i, k-m+1] for m = 1:k) / V[i, 1]
      end
   end
   return V, W
end
println("ready")

# ## Part 1: Section 7 of the theory article
#
# ### 1. Two-bus network with a capacitor bank (Section 7.10)
#
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#710-hand-calculation-two-bus-network-with-an-explicit-shunt-element" target="_blank">Section 7.10 of the theory article</a>
#md # [Section 7.10 of the theory article](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#710-hand-calculation-two-bus-network-with-an-explicit-shunt-element)
# works a network with one unknown by hand. Bus 1 is the slack at
# $1\angle 0°$, bus 2 carries the load $S_2 = -0.5 - j0.15$, the line has the
# series admittance $y = 1 - j4$ and no charging, and a capacitor bank
# $y^{\mathrm{sh}} = +j0.2$ sits at bus 2. The article prints these
# coefficients:
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>n</th><th>V₂⁽ⁿ⁾</th><th>W₂⁽ⁿ⁾</th></tr>
#nb # <tr><td>1</td><td>−0.017647 − j0.120588</td><td>+0.017647 + j0.120588</td></tr>
#nb # <tr><td>2</td><td>−0.016514 + j0.000415</td><td>+0.002284 + j0.003841</td></tr>
#nb # <tr><td>3</td><td>−0.001338 + j0.000214</td><td>+0.001257 + j0.002113</td></tr>
#nb # </table>
#md # | n | $V_2^{(n)}$ | $W_2^{(n)}$ |
#md # |--:|:--|:--|
#md # | 1 | $-0.017647 - j0.120588$ | $+0.017647 + j0.120588$ |
#md # | 2 | $-0.016514 + j0.000415$ | $+0.002284 + j0.003841$ |
#md # | 3 | $-0.001338 + j0.000214$ | $+0.001257 + j0.002113$ |
#
# and the solution $V_2 = 0.964029 - j0.119926$.
#
# First the hand recursion. The reduced system is one equation: the series
# "matrix" is the scalar $y$, the shunt vector is the scalar $j0.2$.

y = 1.0 - 4.0im                     # series admittance of the line
S2 = -0.5 - 0.15im                  # load at bus 2 (negative = consumption)
ysh2 = 0.2im                        # capacitor bank at bus 2

V2, W2 = hand_recursion([y;;], [ysh2], [S2], 10)   # 1×1 series matrix, 10 orders
println("hand recursion, coefficients of bus 2:")
for n = 1:3
   @printf("  n = %d   V2 = %s   W2 = %s\n", n, c(V2[1, n+1]), c(W2[1, n+1]))
end
println("partial sums at s = 1 (the voltage after N orders):")
for N in (1, 2, 3, 4, 10)
   s = sum(V2[1, 1:N+1])
   @printf("  N = %2d   %s   |V2| = %.4f   angle = %.2f°\n", N, c(s), abs(s), rad2deg(angle(s)))
end

# The coefficients match the article to the printed digits, and the partial
# sums approach $0.964029 - j0.119926$.
#
# #### The same calculation with the solver: `apslf_pq`
#
# `apslf_pq` is the PQ-only solver of the package: one slack bus with fixed
# voltage, every other bus a PQ bus with a specified complex injection.
#
# ```
# V, Vcoeff, Wcoeff = apslf_pq(Y, S; slack = 1, order = 24, use_pade = true, germ = :deviation)
# ```
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>Argument</th><th>Meaning</th></tr>
#nb # <tr><td><code>Y</code></td><td>full bus admittance matrix, slack included</td></tr>
#nb # <tr><td><code>S</code></td><td>complex injections per bus (slack entry ignored), load negative</td></tr>
#nb # <tr><td><code>slack</code></td><td>index of the slack bus</td></tr>
#nb # <tr><td><code>order</code></td><td>number of series coefficients</td></tr>
#nb # <tr><td><code>use_pade</code></td><td>evaluate the series at s=1 with a Padé approximant instead of the plain sum</td></tr>
#nb # <tr><td><code>germ</code></td><td>the order-0 state and embedding, see the next cell</td></tr>
#nb # </table>
#md # | Argument | Meaning |
#md # |:--|:--|
#md # | `Y` | full bus admittance matrix, slack included |
#md # | `S` | complex injections per bus (slack entry ignored), load negative |
#md # | `slack` | index of the slack bus |
#md # | `order` | number of series coefficients |
#md # | `use_pade` | evaluate the series at $s=1$ with a Padé approximant instead of the plain sum |
#md # | `germ` | the order-0 state and embedding, see the next cell |
#
# It returns the voltage vector `V` at $s = 1$ (all buses) and the
# coefficient matrices `Vcoeff` and `Wcoeff` (non-slack buses × orders), the
# same layout as `hand_recursion`.
#
# The solver takes the **full** Y-bus, so the capacitor sits in the diagonal
# entry of bus 2. The default `germ = :deviation` moves the row sums of `Y`
# (here exactly the capacitor) to the right-hand side and ramps them with
# $s$. That is the split of Section 7.10, so the solver coefficients must
# agree with the hand recursion.

Y2 = [y -y; -y y+ysh2]              # full 2×2 Y-bus, capacitor in the diagonal of bus 2
## Injections, one entry per bus, same order as the rows of Y2.
## Bus 1 is the slack: its voltage is fixed and its power follows from the network,
## so the entry is never read (0 is a placeholder). Bus 2 carries the load S2 = -0.5 - j0.15.
S = [0.0im, S2]                     # = [0.0 + 0.0im, -0.5 - 0.15im]

## Y2, S : the 2×2 Y-bus and the injection vector from the cell above
## germ = :deviation is the default; written out here so the embedding is visible
## Vs  : voltage at s = 1 for both buses
## Vc  : coefficients V^(n) of bus 2, columns are the orders 0..10
## Wc  : coefficients W^(n) of the reciprocal series
Vs, Vc, Wc = A.apslf_pq(Y2, S; slack = 1, order = 10, use_pade = true, germ = :deviation)

println("solver coefficients against the hand recursion:")
for n = 1:3
   @printf("  n = %d   V2 = %s   W2 = %s   |Δ to hand| = %.1e\n", n, c(Vc[1, n+1]), c(Wc[1, n+1]), abs(Vc[1, n+1] - V2[1, n+1]))
end
@printf("solver V2 at s = 1: %s   article: +0.964029 -0.119926im   |Δ| = %.1e\n", c(Vs[2]), abs(Vs[2] - (0.964029 - 0.119926im)))

## Residual of the physical load-flow equation at bus 2: conj(V2) * I2 - conj(S2)
## with I2 = y (V2 - V1) + ysh2 V2 the current leaving bus 2. Zero means "is a solution".
residual(V) = abs(conj(V[2]) * (y * (V[2] - 1) + ysh2 * V[2]) - conj(S2))
@printf("residual of the physical equation: %.1e\n", residual(Vs))

# #### What the germ is, and why it matters here
#
# The **germ** is the state at $s = 0$, the point where the power series
# starts: $V^{(0)} = V(s{=}0)$. The recursion only works if the germ solves
# the $s = 0$ equations exactly
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#22-embedding-parameter-s" target="_blank">Section 2.2</a>).
#md # ([Section 2.2](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#22-embedding-parameter-s)).
# With no loads and no shunts the flat state $V = 1$ everywhere is such a
# solution. The capacitor breaks this: at $1$ pu it draws current, so
# $V_2 = 1$ no longer satisfies the $s = 0$ network equation. The `germ`
# keyword offers three ways to deal with it
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#65-transformers-and-phase-shifters" target="_blank">Section 6.5</a>):
#md # ([Section 6.5](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#65-transformers-and-phase-shifters)):
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>germ</th><th>Order-0 state V⁽⁰⁾</th><th>Where the capacitor goes</th><th>Exact at s = 0?</th></tr>
#nb # <tr><td><code>:flat</code></td><td>1 at every bus</td><td>stays in the constant matrix</td><td>no (legacy behaviour before 0.9.15)</td></tr>
#nb # <tr><td><code>:noload</code></td><td>no-load voltage, solution of the linear network without loads</td><td>stays in the constant matrix</td><td>yes (Section 6.5, variant 2)</td></tr>
#nb # <tr><td><code>:deviation</code></td><td>1 at every bus</td><td>moved to the right-hand side, ramped with s</td><td>yes (Section 6.5, variant 1; the split of Section 7.10; default)</td></tr>
#nb # </table>
#md # | germ | Order-0 state $V^{(0)}$ | Where the capacitor goes | Exact at $s=0$? |
#md # |:--|:--|:--|:--|
#md # | `:flat` | 1 at every bus | stays in the constant matrix | no (legacy behaviour before 0.9.15) |
#md # | `:noload` | no-load voltage, solution of the linear network without loads | stays in the constant matrix | yes (Section 6.5, variant 2) |
#md # | `:deviation` | 1 at every bus | moved to the right-hand side, ramped with $s$ | yes (Section 6.5, variant 1; the split of Section 7.10; default) |
#
# The next cell runs all three on the two-bus network and prints the germ,
# the voltage the series arrives at, and the residual of the physical
# equation. Read the residual column: a series that starts from a germ that
# is not exact ends at a state that is **not** a load-flow solution.

println("germ        V2^(0)                   V2 at s = 1              residual")
for germ in (:flat, :noload, :deviation)
   Vg, Vcg, _ = A.apslf_pq(Y2, S; slack = 1, order = 10, germ = germ)
   @printf("%-11s %-24s %-24s %.1e\n", germ, c(Vcg[1, 1]), c(Vg[2]), residual(Vg))
end

# - `:flat` starts at $V_2^{(0)} = 1$, which is not a solution at $s = 0$.
#   The series still converges, but to a state with a residual of about
#   $0.2$: it is not the load flow. This is the remark "What happens without
#   the split" in Section 7.10.
# - `:noload` starts at the no-load voltage $1.0492 - j0.0130$ that the
#   article mentions (the capacitor raises the unloaded bus above 1 pu) and
#   arrives at the correct solution.
# - `:deviation` keeps the flat germ, ramps the capacitor with $s$, and
#   arrives at the same solution. This is the article's split and the
#   solver default.

# ### 2. Four-bus π-model network (Sections 7.1 to 7.8)
#
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#71-network-and-data" target="_blank">Section 7.1 of the theory article</a>
#md # [Section 7.1 of the theory article](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#71-network-and-data)
# defines a four-bus network with five π-model lines. Bus 1 is the slack at
# $1\angle 0°$; the loads are $S_2 = -0.4 - j0.15$, $S_3 = -0.5 - j0.175$
# and $S_4 = -0.3 - j0.1$.
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>Line</th><th>Series admittance yᵢₖ</th><th>Total line charging j bᵢₖ</th><th>Half-shunt per side</th></tr>
#nb # <tr><td>1–2</td><td>2 − j6</td><td>j0.06</td><td>j0.03</td></tr>
#nb # <tr><td>1–3</td><td>1 − j3</td><td>j0.04</td><td>j0.02</td></tr>
#nb # <tr><td>2–3</td><td>1.5 − j4.5</td><td>j0.05</td><td>j0.025</td></tr>
#nb # <tr><td>2–4</td><td>1 − j3</td><td>j0.04</td><td>j0.02</td></tr>
#nb # <tr><td>3–4</td><td>1.2 − j3.6</td><td>j0.06</td><td>j0.03</td></tr>
#nb # </table>
#md # | Line | Series admittance $y_{ik}$ | Total line charging $jb_{ik}$ | Half-shunt per side |
#md # |:--|:--|:--|:--|
#md # | 1–2 | $2 - j6$ | $j0.06$ | $j0.03$ |
#md # | 1–3 | $1 - j3$ | $j0.04$ | $j0.02$ |
#md # | 2–3 | $1.5 - j4.5$ | $j0.05$ | $j0.025$ |
#md # | 2–4 | $1 - j3$ | $j0.04$ | $j0.02$ |
#md # | 3–4 | $1.2 - j3.6$ | $j0.06$ | $j0.03$ |
#
# #### Building the Y-bus: `build_ybus_from_branches`
#
# `build_ybus_from_branches(nbus, branches)` stamps π-model lines into a
# dense `nbus × nbus` Y-bus. Each branch is a tuple `(i, j, r, x, b_total)`:
# the two bus numbers, the series **impedance** $z = r + jx$ and the total
# charging susceptance, which it splits equally on both ends. The article
# gives series **admittances**, so each $y_{ik}$ is inverted first:
# $z_{ik} = 1 / y_{ik}$.

## Lines as in the article: (from, to, series admittance y_ik, total charging b_ik)
lines = [
   (1, 2, 2.0 - 6.0im, 0.06),
   (1, 3, 1.0 - 3.0im, 0.04),
   (2, 3, 1.5 - 4.5im, 0.05),
   (2, 4, 1.0 - 3.0im, 0.04),
   (3, 4, 1.2 - 3.6im, 0.06),
]

## build_ybus_from_branches wants (i, j, r, x, b_total) with the series IMPEDANCE
## z = r + jx = 1 / y_ik. inv(y_ik) does the conversion; real/imag split it into r and x.
branches = NTuple{5,Float64}[(i, j, real(inv(yik)), imag(inv(yik)), b) for (i, j, yik, b) in lines]

Y4 = A.build_ybus_from_branches(4, branches)       # dense 4×4 physical Y-bus
S4 = ComplexF64[0, -0.4 - 0.15im, -0.5 - 0.175im, -0.3 - 0.1im]   # injections per bus; bus 1 is the slack, its entry is a placeholder

show_matrix("physical Y-bus (Section 7.2):", Y4)
println()
println("article, row 2: -2+j6   4.5-j13.425   -1.5+j4.5   -1+j3")

# The off-diagonal entries are $-y_{ik}$; a diagonal entry is the sum of
# the series admittances at the bus plus its half-shunts. Row 2 matches the
# matrix printed in
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#72-physical-y-bus-and-why-a-split-is-useful" target="_blank">Section 7.2</a>.
#md # [Section 7.2](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#72-physical-y-bus-and-why-a-split-is-useful).
#
# #### The split into series and shunt part: `apslf_row_sums`
#
# `apslf_row_sums(Y)` returns the row sums $Y \cdot \mathbf{1}$, one complex
# number per bus. A row sum is the current a bus injects when **every** bus
# sits at $1\angle 0°$: the series terms cancel in that state, so what is
# left is exactly the shunt admittance at the bus (here the half-shunts of
# the lines meeting there). The article calls this the diagonal of the shunt
# matrix (Section 7.2). Subtracting it from the diagonal gives the series
# matrix with zero row sums, the constant part of the embedding; the shunt
# vector is ramped with $s$ on the right-hand side. That is what
# `germ = :deviation` does internally.

Ysh = A.apslf_row_sums(Y4)       # row sums = shunt vector = diagonal of the shunt matrix
for i = 1:4
   @printf("bus %d:  row sum %s\n", i, c(Ysh[i]))
end
println("check bus 2: half-shunts 0.03 + 0.025 + 0.02 = j0.075")

# #### Hand recursion against the article's coefficients
#
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#73-reduced-system-for-the-non-slack-buses" target="_blank">Section 7.3</a>
#md # [Section 7.3](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#73-reduced-system-for-the-non-slack-buses)
# removes the slack bus and works with the $3 \times 3$ series matrix and the
# shunt vector of buses 2 to 4. Sections 7.5 to 7.7 print these
# coefficients:
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>n</th><th>V₂⁽ⁿ⁾</th><th>V₃⁽ⁿ⁾</th><th>V₄⁽ⁿ⁾</th></tr>
#nb # <tr><td>1</td><td>−0.057332 − j0.103422</td><td>−0.072835 − j0.130656</td><td>−0.086243 − j0.156913</td></tr>
#nb # <tr><td>2</td><td>−0.019635 + j0.000848</td><td>−0.025731 + j0.001119</td><td>−0.031963 + j0.001265</td></tr>
#nb # <tr><td>3</td><td>−0.003137 + j0.000151</td><td>−0.004157 + j0.000266</td><td>−0.005249 + j0.000461</td></tr>
#nb # </table>
#md # | n | $V_2^{(n)}$ | $V_3^{(n)}$ | $V_4^{(n)}$ |
#md # |--:|:--|:--|:--|
#md # | 1 | $-0.057332 - j0.103422$ | $-0.072835 - j0.130656$ | $-0.086243 - j0.156913$ |
#md # | 2 | $-0.019635 + j0.000848$ | $-0.025731 + j0.001119$ | $-0.031963 + j0.001265$ |
#md # | 3 | $-0.003137 + j0.000151$ | $-0.004157 + j0.000266$ | $-0.005249 + j0.000461$ |

red = 2:4                                        # the non-slack buses
Yser_red = Y4[red, red] - Diagonal(Ysh[red])     # series matrix: shunts taken out of the diagonal
V4h, W4h = hand_recursion(Yser_red, Ysh[red], S4[red], 40)

## The article's numbers, to compare against
article = Dict(
   1 => [-0.057332 - 0.103422im, -0.072835 - 0.130656im, -0.086243 - 0.156913im],
   2 => [-0.019635 + 0.000848im, -0.025731 + 0.001119im, -0.031963 + 0.001265im],
   3 => [-0.003137 + 0.000151im, -0.004157 + 0.000266im, -0.005249 + 0.000461im],
)
for n = 1:3
   @printf("n = %d   V2 = %s   V3 = %s   V4 = %s   max |Δ to article| = %.1e\n", n, c(V4h[1, n+1]), c(V4h[2, n+1]), c(V4h[3, n+1]), maximum(abs.(V4h[:, n+1] .- article[n])))
end
println("W^(3): ", join(c.(W4h[:, 4]), "  "))
println("article: 0.003912+0.003727im  0.004853+0.006144im  0.005431+0.008989im")

# The differences of a few $10^{-7}$ are the rounding of the article to six
# decimals. Next the partial sum after order 3 and the converged value of
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#78-evaluation-at-s1-after-order-3" target="_blank">Section 7.8</a>:
#md # [Section 7.8](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#78-evaluation-at-s1-after-order-3):

V3sum = [sum(V4h[i, 1:4]) for i = 1:3]           # orders 0..3
println("partial sum N = 3: ", join(c.(V3sum), "  "))
println("article:           0.9199-0.1024im  0.8973-0.1293im  0.8765-0.1552im")
V40 = [sum(V4h[i, :]) for i = 1:3]               # orders 0..40
println("sum of 40 terms:   ", join(c.(V40), "  "))
println("article:           0.918383-0.102390im  0.895249-0.129200im  0.873961-0.155029im")

# #### The solver on the four-bus network
#
# Same call as in the two-bus case, now with the $4 \times 4$ Y-bus. With
# the default `germ = :deviation` the solver builds the same series matrix
# and shunt vector internally, so its coefficients must equal the hand
# recursion to machine precision. `calc_injections(Y, V)` computes the
# complex power $S_i = V_i \cdot \overline{(Y V)_i}$ every bus injects at the
# solved voltages; comparing it with the specified `S4` is the power
# mismatch on the physical network.

## Y4 : the 4×4 physical Y-bus built above with build_ybus_from_branches
## S4 : the injections [0, S2, S3, S4] defined next to it (bus 1 is the slack)
## germ = :deviation is the default; written out here so the embedding is visible
## Vsol : voltages at s = 1 (Padé), Vc4 : coefficient matrix (3 non-slack buses × 41 orders)
Vsol, Vc4, _ = A.apslf_pq(Y4, S4; slack = 1, order = 40, use_pade = true, germ = :deviation)
for n = 1:3
   @printf("solver n = %d   max |Δ to hand| = %.1e   max |Δ to article| = %.1e\n", n, maximum(abs.(Vc4[:, n+1] .- V4h[:, n+1])), maximum(abs.(Vc4[:, n+1] .- article[n])))
end
@printf("solver at s = 1 (Padé): V2 = %s   V3 = %s   V4 = %s\n", c(Vsol[2]), c(Vsol[3]), c(Vsol[4]))

Sinj = A.calc_injections(Y4, Vsol)               # power each bus injects at the solved voltages
@printf("power mismatch on the physical Y-bus: %.1e pu\n", maximum(abs.(Sinj[red] .- S4[red])))

# ### 3. The same check with PV buses
#
# The article's examples are PQ-only. The package also has a **direct PV
# formulation**
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#63-direct-pv-formulation-augmented-real-system" target="_blank">Section 6.3</a>):
#md # ([Section 6.3](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#63-direct-pv-formulation-augmented-real-system)):
# a PV bus fixes $P$ and $|V|$ and leaves $Q$ as an unknown, and the solver
# computes one augmented real linear system per order. A cheap consistency
# check: make bus 3 a PV bus at exactly the magnitude it reached in the PQ
# solution. The PV solver must then reproduce the PQ voltages and return
# the PQ reactive power $Q_3 = -0.175$.
#
# `apslf_pf_pv_direct(Y, bustype, P, Q, Vm; ...)` is the dense kernel,
# `apslf_pf_pv_direct_sparse` the sparse one; both take the bus types, the
# specified $P$ and $Q$ (the $Q$ of PV buses is ignored) and the voltage
# setpoints `Vm` (only used at PV buses). They return the voltages, the
# reactive injections of the PV buses and diagnostics.

Vm3 = abs(Vsol[3])                               # |V3| of the PQ solution becomes the setpoint
bt = [:slack, :pq, :pv, :pq]                     # bus 3 is now a PV bus
P = real.(S4)
Q = imag.(S4)
Vm = [1.0, 1.0, Vm3, 1.0]                        # setpoints; only the PV entry is used
for kern in (A.apslf_pf_pv_direct, A.apslf_pf_pv_direct_sparse)
   Vpv, Qpv, _, _, _ = kern(Y4, bt, P, Q, Vm; slack = 1, order = 40, self_check = false)
   @printf("%-26s max |V - V_pq| = %.1e   Q3 = %.6f (PQ case: %.6f)\n", nameof(kern), maximum(abs.(Vpv .- Vsol)), Qpv[1], Q[3])
end

# All three checks agree with the article to the printed digits, and the
# solver reaches the converged values of Section 7.8 without any Newton
# polish.

# ## Part 2: the 9-bus teaching case
#
# Part 1 used `apslf_pq`, the bare PQ recursion. Real cases have PV buses
# (generators that hold a voltage magnitude), reactive limits and a slack
# that is not at 1 pu. Part 2 uses the high-level entry point
# `solve_pf_apslf` on the 9-bus teaching case that ships with the package
# and looks at what the solver does step by step.
#
# ### 4. Solving a case: `solve_pf_apslf`
#
# #### What goes in: the case
#
# `solve_pf_apslf(case; ...)` works on a small NamedTuple data contract:
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>Field</th><th>Meaning</th></tr>
#nb # <tr><td><code>Y</code></td><td>bus admittance matrix, dense or sparse</td></tr>
#nb # <tr><td><code>bustype</code></td><td><code>:slack</code>, <code>:pv</code> or <code>:pq</code> per bus</td></tr>
#nb # <tr><td><code>Pspec</code>, <code>Qspec</code></td><td>specified injections in pu, generation positive, load negative</td></tr>
#nb # <tr><td><code>Vm</code></td><td>voltage setpoints, used at the slack and the PV buses</td></tr>
#nb # <tr><td><code>Qmin</code>, <code>Qmax</code></td><td>reactive limits of the PV buses</td></tr>
#nb # <tr><td><code>slack</code></td><td>index of the slack bus</td></tr>
#nb # </table>
#md # | Field | Meaning |
#md # |:--|:--|
#md # | `Y` | bus admittance matrix, dense or sparse |
#md # | `bustype` | `:slack`, `:pv` or `:pq` per bus |
#md # | `Pspec`, `Qspec` | specified injections in pu, generation positive, load negative |
#md # | `Vm` | voltage setpoints, used at the slack and the PV buses |
#md # | `Qmin`, `Qmax` | reactive limits of the PV buses |
#md # | `slack` | index of the slack bus |
#
# The three bus types differ in what is **given** and what is **unknown**:
#
# - **slack**: voltage magnitude and angle given, $P$ and $Q$ unknown (it
#   balances the network).
# - **PV**: $P$ and $|V|$ given, $Q$ and the angle unknown. A generator
#   with voltage control.
# - **PQ**: $P$ and $Q$ given, the complex voltage unknown. Loads and
#   passive buses.
#
# `demo_case_9bus()` builds the case. Print it first, so the numbers below
# have a meaning: the slack at 1.04 pu, two PV generators at 1.025 pu with
# narrow reactive bands, three loads and three passive buses.

case = A.demo_case_9bus()

println("bus   type   Pspec    Qspec    Vm      Qmin    Qmax")
for i = 1:9
   ## the ±1e9 limits of slack and PQ buses mean "no limit"; print them as a dash
   band = case.bustype[i] == :pv ? @sprintf("%6.2f  %6.2f", case.Qmin[i], case.Qmax[i]) : "   -       -"
   @printf("%-5s %-6s %6.3f   %6.3f   %.3f  %s\n", case.labels[i], case.bustype[i], case.Pspec[i], case.Qspec[i], case.Vm[i], band)
end
@printf("\nlargest row sum of Y: %.3f pu (line charging, so the flat germ is not exact, see Section 6)\n", maximum(abs.(A.apslf_row_sums(case.Y))))

# #### Solving
#
# The keywords below are the ones you will touch most often:
#
# - `mode`: how PV buses are handled, `:direct` (augmented real system per
#   order) or `:outer` (repeated PQ solves), see Section 7.
# - `order`: number of series coefficients. 40 is generous for a small case.
# - `use_pade`: evaluate the series at $s = 1$ with a Padé approximant.
# - `nr_polish`: a few Newton-Raphson steps **after** the series. Switched
#   off here on purpose: everything you see is the pure series result.
#
# The result is a NamedTuple. The fields used in this notebook:
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>Field</th><th>Meaning</th></tr>
#nb # <tr><td><code>V</code></td><td>complex bus voltages at s = 1</td></tr>
#nb # <tr><td><code>Q</code></td><td>reactive injections; at PV buses the value the solver found</td></tr>
#nb # <tr><td><code>bustype</code></td><td>bus types <em>after</em> limit switching (a limited PV bus shows up as <code>:pq</code>)</td></tr>
#nb # <tr><td><code>converged</code></td><td>the series evaluation succeeded (says nothing about the physics, see Section 6)</td></tr>
#nb # <tr><td><code>effective_mode</code></td><td>the PV mode actually used</td></tr>
#nb # <tr><td><code>outer_iters</code></td><td>number of outer passes (limit switching needs at least two)</td></tr>
#nb # <tr><td><code>switch_log</code></td><td>one entry per PV→PQ switch</td></tr>
#nb # <tr><td><code>Vcoeff</code></td><td>the coefficient matrix, only with <code>return_coeffs = true</code></td></tr>
#nb # </table>
#md # | Field | Meaning |
#md # |:--|:--|
#md # | `V` | complex bus voltages at $s = 1$ |
#md # | `Q` | reactive injections; at PV buses the value the solver found |
#md # | `bustype` | bus types *after* limit switching (a limited PV bus shows up as `:pq`) |
#md # | `converged` | the series evaluation succeeded (says nothing about the physics, see Section 6) |
#md # | `effective_mode` | the PV mode actually used |
#md # | `outer_iters` | number of outer passes (limit switching needs at least two) |
#md # | `switch_log` | one entry per PV→PQ switch |
#md # | `Vcoeff` | the coefficient matrix, only with `return_coeffs = true` |

## mode = :direct   → PV buses via the augmented real system per order (Section 6.3)
## order = 40       → 40 series coefficients per bus
## use_pade = true  → Padé evaluation at s = 1
## nr_polish = false → no Newton step afterwards: pure series result
res = solve_pf_apslf(case; mode = :direct, order = 40, use_pade = true, nr_polish = false)

println("converged = $(res.converged), effective mode = $(res.effective_mode), outer iterations = $(res.outer_iters)")

# #### What comes out, and how to check it
#
# The table shows the solved voltages together with the reactive power of
# the generators. Two things to look at:
#
# - **Bus 1** stays at 1.04 pu and 0°, it is the slack.
# - **Bus 2** holds its setpoint 1.025 pu; **bus 3** does not (1.039 pu).
#   Its reactive limit was hit and it was switched to a PQ bus at
#   $Q_{\min}$. Section 7 looks at this in detail.
#
# The **power mismatch** is the check that matters: at every bus, take the
# specified $P + jQ$ and subtract the power that actually flows into the
# network at the solved voltages, $S_i = V_i \cdot \overline{(Y V)_i}$.
# `compute_demo_mismatch(case, res)` does exactly that on the case's
# physical Y-bus. A mismatch at machine precision means the series alone
# solves the load-flow equations, no Newton step involved.

## largest power mismatch of a result on the physical Y-bus, in pu
mismatch(case, res) = maximum(A.compute_demo_mismatch(case, res))

println("bus   type   |V| (pu)   angle       Q (pu)")
for i = 1:9
   ## res.bustype is the type after limit switching; res.Q the reactive injection found
   note = case.bustype[i] == :pv && res.bustype[i] == :pq ? "  ← PV switched to PQ (Q limit)" : ""
   @printf("%-5s %-6s %.5f   %8.4f°   %7.4f%s\n", case.labels[i], res.bustype[i], abs(res.V[i]), rad2deg(angle(res.V[i])), res.Q[i], note)
end
@printf("\nmax power mismatch on the physical Y-bus: %.2e pu\n", mismatch(case, res))

# ### 5. The series behind the numbers
#
# Every voltage is a power series in the embedding parameter,
# $V_i(s) = \sum_n V_i^{(n)} s^n$, and the load flow is its value at
# $s = 1$. With `return_coeffs = true` the result carries the coefficient
# matrix `Vcoeff` (buses × orders), so the series can be inspected.
#
# The first question about any power series is its **convergence radius**
# $R$: the series converges for $|s| < R$, and it is only useful if $R > 1$.
# For large $n$ the coefficients shrink geometrically, $|V^{(n+1)}| /
# |V^{(n)}| \to 1/R$, so the ratio of consecutive coefficients estimates
# $R$.

res = solve_pf_apslf(case; order = 40, nr_polish = false, return_coeffs = true)
Vcoef = res.Vcoeff                               # 9 buses × 41 orders, column n+1 holds order n

println("decay of the coefficients:")
for n in (0, 1, 2, 5, 10, 20, 40)
   @printf("  n = %2d   max_i |V_i^(n)| = %.2e\n", n, maximum(abs.(Vcoef[:, n+1])))
end

## ratio of consecutive coefficient norms → 1/R
ratio = maximum(abs.(Vcoef[:, 41])) / maximum(abs.(Vcoef[:, 40]))
@printf("\nratio |V^(40)| / |V^(39)| = %.3f   →  convergence radius R ≈ %.2f\n", ratio, 1 / ratio)

# $R \approx 2.5$ means: the operating point $s = 1$ is well inside the
# disc of convergence, and every additional order shrinks the next
# coefficient by a factor of about 2.5. The coefficient of order 40 is
# around $10^{-18}$, far below machine precision, so it contributes nothing
# to the sum at $s = 1$. In other words, `order = 40` is more than this
# case needs.
#
# How many orders **are** needed? The truncation error of the Taylor sum is
# roughly the size of the first omitted coefficient. With a factor 2.5 per
# order, $10^{-10}$ takes about 25 orders and machine precision about 40.
# The Padé evaluation gets there earlier, because the rational function
# extrapolates the tail of the series. The next cell solves the case with
# increasing order and prints the mismatch for both evaluations.

println("order   mismatch Taylor   mismatch Padé")
for o in (5, 10, 15, 20, 30, 40)
   rt = solve_pf_apslf(case; order = o, nr_polish = false, use_pade = false)   # plain sum
   rpd = solve_pf_apslf(case; order = o, nr_polish = false, use_pade = true)   # Padé
   @printf("  %2d      %.1e          %.1e\n", o, mismatch(case, rt), mismatch(case, rpd))
end

# Rule of thumb for this case: with Padé, `order = 20` is accurate to about
# $10^{-13}$ and `order = 30` reaches machine precision. Harder cases (a
# smaller $R$, pole closer to $s = 1$) need more orders, and there the
# advantage of Padé grows.
#
# #### Taylor vs Padé
#
# `evaluate_series(coeffs, options)` evaluates the series of **one** bus at
# $s = 1$. `APSLFEvaluationOptions(mode = :taylor)` sums the polynomial;
# `mode = :pade` turns the polynomial into a rational function, the Padé
# approximant of
#nb # <a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#52-padé-approximation-from-series-to-quotient" target="_blank">Section 5.2</a>.
#md # [Section 5.2](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#52-padé-approximation-from-series-to-quotient).
# A rational function can represent the voltage **beyond** the convergence
# radius of the polynomial (analytic continuation), which is why the
# solver uses Padé by default. On this easy case both must agree.

c5 = Vcoef[5, :]                                 # the 41 coefficients of bus 5 (a load bus)
taylor = A.evaluate_series(c5, A.APSLFEvaluationOptions(mode = :taylor)).voltage   # plain sum
pade = A.evaluate_series(c5, A.APSLFEvaluationOptions(mode = :pade)).voltage       # rational approximant
@printf("bus 5:  Taylor %.8f ∠ %.5f°   Padé %.8f ∠ %.5f°   |Δ| = %.1e\n", abs(taylor), rad2deg(angle(taylor)), abs(pade), rad2deg(angle(pade)), abs(taylor - pade))

# #### The poles as a distance-to-collapse indicator
#
# A Padé approximant is a quotient of two polynomials, and the roots of the
# denominator are its **poles**. They approximate the singularities of the
# true voltage function, and the nearest singularity is where the load
# flow ceases to exist (the "nose" of the PV curve). Its distance to
# $s = 1$ is therefore a heuristic margin to collapse
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#53-properties-and-practical-use" target="_blank">Section 5.3</a>).
#md # ([Section 5.3](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#53-properties-and-practical-use)).
# `stability_from_Vcoeff` computes it from the coefficient matrix and
# returns `dmin` (that distance) and the pole; `st_level` maps `dmin` to a
# traffic light: GRN above 0.3, YEL above 0.1, RED below.
#
# The experiment: scale all injections by a factor and watch the nearest
# pole move towards $s = 1$. Two observations: the distance shrinks as the
# loading grows, and at factor 2 the load flow already has no solution
# while the level is still GRN. The indicator is a heuristic margin, not a
# certificate; read it as a trend, together with `converged` and the
# mismatch.

println("load factor  converged   min |V|   pole distance   level")
for factor in (1.0, 1.5, 2.0, 2.5)
   ## same case, all injections scaled; reactive limits switched off so only the loading changes
   heavy = merge(case, (Pspec = factor .* case.Pspec, Qspec = factor .* case.Qspec, Qmin = fill(-1e9, 9), Qmax = fill(1e9, 9)))
   rh = solve_pf_apslf(heavy; order = 40, nr_polish = false, return_coeffs = true)
   st = A.stability_from_Vcoeff(rh.Vcoeff; slack = 1, order = 40)   # st.dmin = distance of the nearest pole to s = 1
   @printf("   × %.1f      %-9s   %.4f    %.3f           %s\n", factor, rh.converged, minimum(abs.(rh.V)), st.dmin, A.st_level(st.dmin))
end

# ### 6. The germ on the 9-bus case
#
# Section 1 showed on two buses what the germ is: the state at $s = 0$
# from which the series starts, and it has to be an **exact** solution of
# the $s = 0$ equations. On the 9-bus case two things break the flat germ
# $V = 1$: the line charging (non-zero row sums of `Y`, printed in
# Section 4) and the slack at 1.04 pu instead of 1.
#
# `solve_pf_apslf` takes the same `germ` keyword as `apslf_pq`:
#
# - `:deviation` (default): flat germ at the slack voltage, row sums ramped
#   with $s$ (variant 1 of Section 6.5).
# - `:noload`: start from the no-load voltages (variant 2).
# - `:flat`: the plain flat germ on the full `Y`, the behaviour before
#   0.9.15. Not exact here.
#
# What the experiment shows: read `converged` and the mismatch
# **together**. `converged` only reports that the series evaluation
# succeeded, that the coefficients decayed and the Padé approximant could
# be built. It says nothing about whether the result solves the network.
# With the flat germ the series does converge, but to the solution of a
# **different** problem (the one whose $s = 0$ state is $V = 1$), and that
# state is 0.6 pu away from the load flow, with a bus at 0.27 pu. Only the
# mismatch reveals it.
#
# Which germ to use: keep the default `:deviation`. It is exact for any
# `Y`, keeps the germ at nominal voltage and usually has the larger
# convergence radius. `:noload` is the alternative to try when a case with
# strong shunts or transformers does not converge with the default.
# `:flat` exists to reproduce results of versions before 0.9.15; do not use
# it for new work.

println("germ        converged   max mismatch (pu)   min |V| (pu)")
for germ in (:deviation, :noload, :flat)
   rg = solve_pf_apslf(case; order = 40, nr_polish = false, germ = germ)
   @printf("%-11s %-11s %.2e            %.3f\n", germ, rg.converged, mismatch(case, rg), minimum(abs.(rg.V)))
end

# The **Newton polish** (`nr_polish = true`) runs a few Newton-Raphson
# steps from the series result. Before 0.9.15 it was what turned the
# flat-germ result into a solution. With an exact germ the series is
# already at machine precision and the polish has nothing left to do; the
# score is the mismatch before and after.

rp = solve_pf_apslf(case; order = 40, nr_polish = true)              # default germ = :deviation
@printf("with polish: mismatch before %.1e, after %.1e, improved = %s\n", rp.nr_polish_score_before, rp.nr_polish_score_after, rp.nr_polish_improved)

# ### 7. PV buses and reactive limits
#
# A PV bus holds $|V|$ by injecting or absorbing reactive power. A real
# generator can only do that within a band $[Q_{\min}, Q_{\max}]$. When
# the band is exhausted the bus can no longer hold its voltage: it is
# switched to a PQ bus at the violated limit, and its voltage becomes an
# unknown like at any load bus.
#
# First without limit enforcement, to see what bus 3 **would** need:

## enforce_q_limits = false → PV buses hold their setpoint whatever Q that takes
rq = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = false)
@printf("without limits: |V3| = %.4f (setpoint %.3f)   Q3 = %.4f pu   band [%.2f, %.2f]\n", abs(rq.V[3]), case.Vm[3], rq.Q[3], case.Qmin[3], case.Qmax[3])

# Bus 3 would have to absorb more reactive power than its lower limit
# allows. With limits enforced (the default) the solver detects this in
# the first outer pass, fixes $Q_3 = Q_{\min}$, declares bus 3 a PQ bus
# and solves again. That second pass is why `outer_iters = 2` in
# Section 4, and why bus 3 ended above its setpoint: with less reactive
# absorption the voltage rises.
#
# Two ways to handle PV buses
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#6-practical-treatment-of-pv-buses" target="_blank">Section 6</a>):
#md # ([Section 6](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#6-practical-treatment-of-pv-buses)):
#
# - **`mode = :direct`**: $Q$ of the PV buses is an unknown of the series
#   itself; per order one augmented real linear system (Section 6.3).
# - **`mode = :outer`**: every PV bus is treated as a PQ bus, and an outer
#   loop corrects its $Q$ from one series solve to the next until $|V|$
#   matches (Section 6.2).
#
# Both end at the same solution and both switch bus 3. The $Q$ recorded in
# the switch log differs, because it is the value **at the moment of the
# switch**: the direct mode already knows the exact $Q_3$ of the
# unlimited solution, the outer loop only has its current estimate.
#
# Why does the **direct** mode report two outer iterations? The direct
# formulation removes the outer loop for holding $|V|$, that part is
# inside the series. The reactive **limits** are a different matter: which
# buses end up limited is a yes/no decision that can only be made after a
# solution exists. So pass 1 solves with bus 3 as PV, finds $Q_3$ outside
# the band, switches the bus, and pass 2 solves the changed problem. With
# limits enforced, every mode needs at least one pass per switch.

for mode in (:direct, :outer)
   rm_ = solve_pf_apslf(case; mode = mode, order = 40, nr_polish = false)   # limits enforced (default)
   sw = get(rm_, :switch_log, ())                 # one entry per PV→PQ switch
   @printf("mode = %-7s converged = %-5s outer iterations = %d  switches = %d  |V3| = %.4f  Q3 = %.4f  Q2 = %.4f pu\n", mode, rm_.converged, rm_.outer_iters, length(sw), abs(rm_.V[3]), rm_.Q[3], rm_.Q[2])
   for e in sw
      @printf("   outer pass %d: bus %d hit its %s limit (Q at the switch = %.4f pu)\n", e.outer, e.bus, e.side, e.qinj)
   end
end

# Without limit enforcement the difference between the modes becomes
# visible: the direct mode needs **one** pass, the outer loop needs many,
# because it has to iterate the reactive injections of the PV buses until
# their voltages match, and it stops at a tolerance rather than at machine
# precision.

println()
println("enforce_q_limits = false:")
for mode in (:direct, :outer)
   rn = solve_pf_apslf(case; mode = mode, order = 40, nr_polish = false, enforce_q_limits = false)
   @printf("mode = %-7s outer iterations = %2d  max | |V| - Vm | at PV buses = %.1e  mismatch = %.1e pu\n", mode, rn.outer_iters, maximum(abs.(abs.(rn.V[[2, 3]]) .- case.Vm[[2, 3]])), mismatch(case, rn))
end

# On nine buses with two generators that is a matter of milliseconds
# either way. Section 8 repeats the comparison on 118 buses with eleven
# generators, where the outer loop does not get there at all.

# ### 8. Sparse matrices
#
# Nothing in the method depends on the matrix being dense: the recursion
# solves one linear system per order with the **same** matrix, so one
# sparse factorization is reused for every order. Pass a sparse `Y` and
# `solve_pf_apslf` selects the sparse direct PV kernel (the one checked in
# Section 3) automatically. The synthetic 118-bus case that ships with the
# package illustrates it: 11 PV buses, 106 PQ buses. The solve is timed
# twice, because the first call of a Julia function includes compilation.

using SparseArrays
case118 = A.demo_case_118bus_synthetic()
sparse_case = merge(case118, (Y = sparse(case118.Y),))   # same case, Y stored sparse

rs = solve_pf_apslf(sparse_case; order = 40, nr_polish = false)              # warm-up: includes compilation
t = @elapsed rs = solve_pf_apslf(sparse_case; order = 40, nr_polish = false) # timed run
@printf("118 buses, nnz(Y) = %d (of %d entries)\n", nnz(sparse_case.Y), length(sparse_case.Y))
@printf("converged = %s, mode = %s, max mismatch = %.1e pu, solve time %.3f s\n", rs.converged, rs.effective_mode, mismatch(case118, rs), t)

# The mismatch is at machine precision again: same recursion, same germ,
# same Padé evaluation, only the linear algebra changed.
#
# #### Direct vs outer on 118 buses
#
# The same case in both PV modes. The outer loop has to find eleven
# reactive injections at once by successive correction, and the
# corrections interact through the network. Within its 30 passes it does
# not reach the setpoints; the direct mode solves it in one pass, because
# the eleven unknowns are part of the linear system of every order.

pv118 = findall(==(:pv), case118.bustype)        # the eleven PV buses
println("mode     converged   outer passes   max | |V| - Vm | at PV   time")
for mode in (:direct, :outer)
   rm118 = solve_pf_apslf(sparse_case; mode = mode, order = 40, nr_polish = false)
   t118 = @elapsed rm118 = solve_pf_apslf(sparse_case; mode = mode, order = 40, nr_polish = false)
   @printf("%-8s %-11s %2d             %.1e                 %.3f s\n", mode, rm118.converged, rm118.outer_iters, maximum(abs.(abs.(rm118.V[pv118]) .- case118.Vm[pv118])), t118)
end

# This is the practical reason `mode = :direct` is the default: the outer
# loop is the simpler method to explain, the direct formulation is the one
# to use.
