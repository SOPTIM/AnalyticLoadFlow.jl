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
# > fresh session. This notebook targets Julia ≥ 1.13.

#nb # ## Setup (Colab)
#nb # This cell installs AnalyticLoadFlow from GitHub (branch `main`) into a
#nb # fresh temporary environment. Run it first, once per session. To test a
#nb # branch, change `rev`. For a private checkout use a personal access token:
#nb # `Pkg.add(url = "https://USER:TOKEN@github.com/USER/AnalyticLoadFlow.jl", rev = "main")`.
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
S = [0.0im, S2]                     # injections; the slack entry is ignored

## V   : voltage at s = 1 for both buses
## Vc  : coefficients V^(n) of bus 2, columns are the orders 0..10
## Wc  : coefficients W^(n) of the reciprocal series
Vs, Vc, Wc = A.apslf_pq(Y2, S; slack = 1, order = 10, use_pade = true)

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
S4 = ComplexF64[0, -0.4 - 0.15im, -0.5 - 0.175im, -0.3 - 0.1im]   # injections, slack entry unused

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

## Vsol : voltages at s = 1 (Padé), Vc4 : coefficient matrix (3 non-slack buses × 41 orders)
Vsol, Vc4, _ = A.apslf_pq(Y4, S4; slack = 1, order = 40, use_pade = true)
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
# ### 4. Solving a case: `solve_pf_apslf`
#
# `solve_pf_apslf(case; ...)` is the high-level entry point. It works on a
# small NamedTuple data contract:
#
#nb # <table border="1" cellpadding="6" style="border-collapse:collapse">
#nb # <tr><th>Field</th><th>Meaning</th></tr>
#nb # <tr><td><code>Y</code></td><td>bus admittance matrix, dense or sparse</td></tr>
#nb # <tr><td><code>bustype</code></td><td><code>:slack</code>, <code>:pv</code> or <code>:pq</code> per bus</td></tr>
#nb # <tr><td><code>Pspec</code>, <code>Qspec</code></td><td>injections in pu, load negative</td></tr>
#nb # <tr><td><code>Vm</code></td><td>voltage setpoints (used at the slack and PV buses)</td></tr>
#nb # <tr><td><code>Qmin</code>, <code>Qmax</code></td><td>reactive limits of the PV buses</td></tr>
#nb # <tr><td><code>slack</code></td><td>index of the slack bus</td></tr>
#nb # </table>
#md # | Field | Meaning |
#md # |:--|:--|
#md # | `Y` | bus admittance matrix, dense or sparse |
#md # | `bustype` | `:slack`, `:pv` or `:pq` per bus |
#md # | `Pspec`, `Qspec` | injections in pu, load negative |
#md # | `Vm` | voltage setpoints (used at the slack and PV buses) |
#md # | `Qmin`, `Qmax` | reactive limits of the PV buses |
#md # | `slack` | index of the slack bus |
#
# It returns a NamedTuple with the voltages `V`, the reactive injections
# `Q`, `converged`, the mode that was actually used and iteration counts.
# `demo_case_9bus()` builds the 9-bus teaching case that ships with the
# package: two PV generators, six π-lines with charging, slack at 1.04 pu.
# `compute_demo_mismatch(case, res)` evaluates the power mismatch of a
# result on the case's physical Y-bus.

case = A.demo_case_9bus()

## largest power mismatch of a result on the physical Y-bus, in pu
mismatch(case, res) = maximum(A.compute_demo_mismatch(case, res))

## mode = :direct   → PV buses via the augmented real system per order
## nr_polish = false → no Newton step afterwards: what you see is the pure series result
res = solve_pf_apslf(case; mode = :direct, order = 40, use_pade = true, nr_polish = false)
println("converged = $(res.converged), effective mode = $(res.effective_mode), outer iterations = $(res.outer_iters)")
@printf("max mismatch on the physical Y-bus: %.2e pu (no Newton polish involved)\n", mismatch(case, res))
println()
println("bus   type   |V| (pu)   angle")
for i = 1:9
   @printf("%-5s %-6s %.5f   %8.4f°\n", case.labels[i], case.bustype[i], abs(res.V[i]), rad2deg(angle(res.V[i])))
end

# ### 5. The series behind the numbers
#
# With `return_coeffs = true` the result also carries the coefficient matrix
# `Vcoeff` (buses × orders). The coefficients decay geometrically, which is
# what makes the series usable. `evaluate_series(coeffs, options)` sums the
# series of one bus at $s = 1$, either as the plain Taylor sum or as a Padé
# approximant
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#52-padé-approximation-from-series-to-quotient" target="_blank">Section 5.2</a>),
#md # ([Section 5.2](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#52-padé-approximation-from-series-to-quotient)),
# which turns the polynomial into a rational function and extends the
# reach of the series to cases where the Taylor sum does not converge.

res = solve_pf_apslf(case; order = 40, nr_polish = false, return_coeffs = true)
Vcoef = res.Vcoeff                               # 9 buses × 41 orders
println("decay of the coefficients:")
for n in (0, 1, 2, 5, 10, 20, 40)
   @printf("  n = %2d   max_i |V_i^(n)| = %.2e\n", n, maximum(abs.(Vcoef[:, n+1])))
end

c5 = Vcoef[5, :]                                 # the series of bus 5
taylor = A.evaluate_series(c5, A.APSLFEvaluationOptions(mode = :taylor)).voltage
pade = A.evaluate_series(c5, A.APSLFEvaluationOptions(mode = :pade)).voltage
@printf("bus 5:  Taylor %.8f ∠ %.5f°   Padé %.8f ∠ %.5f°   |Δ| = %.1e\n", abs(taylor), rad2deg(angle(taylor)), abs(pade), rad2deg(angle(pade)), abs(taylor - pade))

# On this well-conditioned case both evaluations agree to machine
# precision. The poles of the Padé denominator carry extra information: the
# distance of the nearest pole to $s = 1$ is a heuristic
# distance-to-collapse indicator
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#53-properties-and-practical-use" target="_blank">Section 5.3</a>).
#md # ([Section 5.3](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#53-properties-and-practical-use)).
# `stability_from_Vcoeff` computes it from the coefficient matrix,
# `st_level` maps it to a traffic-light label. Scaling all injections up
# moves the nearest pole towards $s = 1$:

for factor in (1.0, 1.5, 2.0, 2.5)
   ## same case with all injections scaled and the reactive limits switched off
   heavy = merge(case, (Pspec = factor .* case.Pspec, Qspec = factor .* case.Qspec, Qmin = fill(-1e9, 9), Qmax = fill(1e9, 9)))
   rh = solve_pf_apslf(heavy; order = 40, nr_polish = false, return_coeffs = true)
   st = A.stability_from_Vcoeff(rh.Vcoeff; slack = 1, order = 40)
   @printf("load × %.1f: converged = %-5s  min |V| = %.4f  pole distance = %.3f  %s\n", factor, rh.converged, minimum(abs.(rh.V)), st.dmin, A.st_level(st.dmin))
end

# ### 6. The germ on the 9-bus case
#
# The same question as in Section 1, now on a real case. The row sums of
# the 9-bus matrix are non-zero (line charging), and the slack sits at
# 1.04 pu, so the flat germ on the full Y-bus is not an exact order-0 state.
# `solve_pf_apslf` takes the same `germ` keyword as `apslf_pq`: `:deviation`
# (default, variant 1 of Section 6.5), `:noload` (variant 2) or the legacy
# `:flat`. Both exact variants must give the same solution; the legacy germ
# only arrives at one with the Newton polish (`nr_polish = true`).

println("germ        converged   max mismatch (pu)")
for germ in (:deviation, :noload, :flat)
   rg = solve_pf_apslf(case; order = 40, nr_polish = false, germ = germ)
   @printf("%-11s %-11s %.2e\n", germ, rg.converged, mismatch(case, rg))
end

## Newton polish: a few Newton-Raphson steps starting from the series result.
## The score is the mismatch before and after; on an exact germ there is little to gain.
rp = solve_pf_apslf(case; order = 40, nr_polish = true)
@printf("with polish (germ = :deviation): score before %.1e, after %.1e, improved = %s\n", rp.nr_polish_score_before, rp.nr_polish_score_after, rp.nr_polish_improved)

# ### 7. PV buses and reactive limits
#
# Two ways to handle PV buses
#nb # (<a href="https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#6-practical-treatment-of-pv-buses" target="_blank">Section 6</a>):
#md # ([Section 6](https://github.com/SOPTIM/AnalyticLoadFlow.jl/blob/main/docs/src/theorie-eng.md#6-practical-treatment-of-pv-buses)):
# the **outer loop** (`mode = :outer`) treats every PV bus as a PQ bus and
# adjusts its reactive injection from one series solve to the next until
# the voltage magnitude is met; the **direct formulation** (`mode = :direct`)
# solves the augmented real system of Section 6.3 per order. In both modes
# the reactive limits are enforced: a PV bus that would leave its
# $[Q_{\min}, Q_{\max}]$ band is switched to a PQ bus at the violated limit,
# and the result carries a `switch_log`.

for mode in (:direct, :outer)
   rm_ = solve_pf_apslf(case; mode = mode, order = 40, nr_polish = false)
   sw = get(rm_, :switch_log, ())                 # one entry per PV→PQ switch
   @printf("mode = %-7s converged = %-5s outer iterations = %d  PV→PQ switches = %d  Q(bus 2) = %.4f pu\n", mode, rm_.converged, rm_.outer_iters, length(sw), rm_.Q[2])
   for e in sw
      @printf("   outer %d: bus %d hit its %s limit (Q = %.4f pu)\n", e.outer, e.bus, e.side, e.qinj)
   end
end

## Without limit enforcement bus 3 keeps its setpoint and draws whatever Q that needs
rq = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = false)
@printf("enforce_q_limits = false: Q(bus 3) = %.4f pu, band [%.2f, %.2f]\n", rq.Q[3], case.Qmin[3], case.Qmax[3])

# ### 8. Sparse matrices
#
# For larger networks pass a sparse `Y`; `solve_pf_apslf` then selects the
# sparse direct PV kernel automatically. The recursion, the germ and the
# Padé evaluation are identical, only the linear algebra changes. The
# synthetic 118-bus case that ships with the package illustrates it. The
# solve is timed twice: the first call includes compilation.

using SparseArrays
case118 = A.demo_case_118bus_synthetic()
sparse_case = merge(case118, (Y = sparse(case118.Y),))   # same case, sparse Y-bus
rs = solve_pf_apslf(sparse_case; order = 40, nr_polish = false)              # warm-up (compilation)
t = @elapsed rs = solve_pf_apslf(sparse_case; order = 40, nr_polish = false)
@printf("118 buses, nnz(Y) = %d: converged = %s, mode = %s, max mismatch = %.1e pu, %.3f s\n", nnz(sparse_case.Y), rs.converged, rs.effective_mode, mismatch(case118, rs), t)
