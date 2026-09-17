```@meta
EditURL = "../../lit/workshop_apslf.jl"
```

# APSLF workshop: the theory article against the solver, then the 9-bus case

> **Level: Newcomer to Advanced.** One notebook, one install; runs in about a minute afterwards.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/SOPTIM/AnalyticLoadFlow.jl/blob/main/notebooks/workshop_apslf.ipynb)

> **Note:** This workshop was created with AI assistance and is reviewed
> and curated by the maintainer; it is not a fully machine-generated text.

Version 1.7 of the theory article corrected the right-hand side of the
recursion: it needs the reflected reciprocal $\overline{W^{(n-1)}}$, not
$W^{(n-1)}$ (Section 2.4). The same error was in the solver, hidden by the
Newton polish. This notebook takes the two hand calculations of Section 7,
the two-bus network with a capacitor bank (7.10) and the four-bus π-model
network (7.1 to 7.8), and checks the solver coefficients digit by digit
against the numbers printed in the article. It also reproduces the old,
wrong recursion so that the effect of the error is visible.

> **Note:** On Google Colab the install cell takes a few minutes on a
> fresh session. This notebook targets Julia ≥ 1.12.

## Warm-up and helpers

````@example workshop_apslf
using AnalyticLoadFlow
using LinearAlgebra
using Printf
const A = AnalyticLoadFlow

c(z) = @sprintf("%+.6f %+.6fim", real(z), imag(z))

# The recursion of theory Section 7.3 on a reduced system:
#   Yser_red V^(n) = S* ⊙ conj(W^(n-1)) - Ysh_red V^(n-1),   V^(0) = 1, W^(0) = 1
# `conjugate = false` reproduces the error of the earlier article versions
# and of the solver before version 0.9.15.
function hand_recursion(Yser_red, ysh_red, Sred, order; conjugate = true)
   n = length(Sred)
   V = zeros(ComplexF64, n, order + 1); W = zeros(ComplexF64, n, order + 1)
   V[:, 1] .= 1; W[:, 1] .= 1
   F = lu(Yser_red)
   for k = 1:order
      Wprev = conjugate ? conj.(W[:, k]) : W[:, k]
      V[:, k+1] = F \ (conj.(Sred) .* Wprev .- ysh_red .* V[:, k])
      for i = 1:n
         W[i, k+1] = -sum(V[i, m+1] * W[i, k-m+1] for m = 1:k) / V[i, 1]
      end
   end
   return V, W
end
println("ready")
````

## Part 1: Section 7 of the theory article

### 1. Two-bus network with a capacitor bank (Section 7.10)

Bus 1 slack at $1\angle 0°$, bus 2 with $S_2 = -0.5 - j0.15$, line
$y = 1 - j4$, capacitor $y^{\mathrm{sh}} = +j0.2$ at bus 2. The article gives
the coefficients

| n | $V_2^{(n)}$ | $W_2^{(n)}$ |
|--:|:--|:--|
| 1 | $-0.017647 - j0.120588$ | $0.017647 + j0.120588$ |
| 2 | $-0.016514 + j0.000415$ | $0.002284 + j0.003841$ |
| 3 | $-0.001338 + j0.000214$ | $0.001257 + j0.002113$ |

and the solution $V_2 = 0.964029 - j0.119926$.

````@example workshop_apslf
y = 1.0 - 4.0im
S2 = -0.5 - 0.15im
V2, W2 = hand_recursion([y;;], [0.2im], [S2], 10)
println("hand calculation with the reflected reciprocal (article 1.7):")
for n = 1:3
   @printf("  n = %d   V2 = %s   W2 = %s\n", n, c(V2[1, n+1]), c(W2[1, n+1]))
end
println("partial sums at s = 1:")
for N in (1, 2, 3, 4, 10)
   s = sum(V2[1, 1:N+1])
   @printf("  N = %2d   %s   |V2| = %.4f   angle = %.2f°\n", N, c(s), abs(s), rad2deg(angle(s)))
end
````

The same recursion without the conjugation. Order 1 agrees, because
$W^{(0)}$ is real; from order 2 on the numbers differ, and the limit is a
different function: it solves $V \odot (YV) = S^*$ instead of the physical
$\overline{V} \odot (YV) = S^*$.

````@example workshop_apslf
V2w, _ = hand_recursion([y;;], [0.2im], [S2], 10; conjugate = false)
println("without the conjugation (old article / old solver):")
for n = 1:3
   @printf("  n = %d   V2 = %s\n", n, c(V2w[1, n+1]))
end
Vphys(V) = abs(conj(V[2]) * (y * (V[2] - 1) + 0.2im * V[2]) - conj(S2))   # residual of the physical equation
@printf("residual of the physical equation at s = 1:  with conj %.1e   without conj %.1e\n", Vphys([1, sum(V2[1, :])]), Vphys([1, sum(V2w[1, :])]))
````

Now the solver. The full Y-bus keeps the capacitor in the diagonal; the
solver's default embedding `germ = :deviation` moves the row sums
(here exactly the capacitor) to the right-hand side, which is the split
of Section 7.10. The coefficients must therefore agree with the article.

````@example workshop_apslf
Y2 = [y -y; -y y+0.2im]
Vs, Vc, Wc = A.apslf_pq(Y2, [0.0im, S2]; slack = 1, order = 10, use_pade = true)
println("solver coefficients (apslf_pq, germ = :deviation):")
for n = 1:3
   @printf("  n = %d   V2 = %s   W2 = %s   |Δ to hand| = %.1e\n", n, c(Vc[1, n+1]), c(Wc[1, n+1]), abs(Vc[1, n+1] - V2[1, n+1]))
end
@printf("solver V2 at s = 1: %s   article: +0.964029 -0.119926im   |Δ| = %.1e\n", c(Vs[2]), abs(Vs[2] - (0.964029 - 0.119926im)))
@printf("residual of the physical equation: %.1e\n", Vphys(Vs))
````

With the legacy flat germ on the full Y-bus (`germ = :flat`, the behaviour
before 0.9.15) the order-0 state is not a solution, because the capacitor
draws current at 1 pu. The article's remark "What happens without the
split" describes exactly this; the no-load germ (`germ = :noload`) is the
alternative it mentions.

````@example workshop_apslf
for germ in (:flat, :noload, :deviation)
   Vg, Vcg, _ = A.apslf_pq(Y2, [0.0im, S2]; slack = 1, order = 10, germ = germ)
   @printf("germ = %-10s V2^(0) = %s   V2(1) = %s   residual %.1e\n", germ, c(Vcg[1, 1]), c(Vg[2]), Vphys(Vg))
end
````

### 2. Four-bus π-model network (Sections 7.1 to 7.8)

Five lines with series admittance and total charging from the table in
7.1; loads $S_2 = -0.4 - j0.15$, $S_3 = -0.5 - j0.175$, $S_4 = -0.3 - j0.1$.
`build_ybus_from_branches` takes `(i, j, r, x, b_total)`, so the series
admittances are converted to impedances first.

````@example workshop_apslf
lines = [(1, 2, 2.0 - 6.0im, 0.06), (1, 3, 1.0 - 3.0im, 0.04), (2, 3, 1.5 - 4.5im, 0.05), (2, 4, 1.0 - 3.0im, 0.04), (3, 4, 1.2 - 3.6im, 0.06)]
branches = NTuple{5,Float64}[(i, j, real(inv(yik)), imag(inv(yik)), b) for (i, j, yik, b) in lines]
Y4 = A.build_ybus_from_branches(4, branches)
S4 = ComplexF64[0, -0.4 - 0.15im, -0.5 - 0.175im, -0.3 - 0.1im]
println("physical Y-bus (Section 7.2), row 2: ", join(c.(Y4[2, :]), "  "))
Ysh = A.apslf_row_sums(Y4)
println("row sums = shunt matrix diagonal: ", join(c.(Ysh), "  "))
````

Hand recursion of Section 7.3 with the series matrix and the shunt vector,
compared with the article's coefficients:

| n | $V_2^{(n)}$ | $V_3^{(n)}$ | $V_4^{(n)}$ |
|--:|:--|:--|:--|
| 1 | $-0.057332 - j0.103422$ | $-0.072835 - j0.130656$ | $-0.086243 - j0.156913$ |
| 2 | $-0.019635 + j0.000848$ | $-0.025731 + j0.001119$ | $-0.031963 + j0.001265$ |
| 3 | $-0.003137 + j0.000151$ | $-0.004157 + j0.000266$ | $-0.005249 + j0.000461$ |

````@example workshop_apslf
red = 2:4
Yser_red = Y4[red, red] - Diagonal(Ysh[red])
V4h, W4h = hand_recursion(Yser_red, Ysh[red], S4[red], 40)
article = Dict(
   1 => [-0.057332 - 0.103422im, -0.072835 - 0.130656im, -0.086243 - 0.156913im],
   2 => [-0.019635 + 0.000848im, -0.025731 + 0.001119im, -0.031963 + 0.001265im],
   3 => [-0.003137 + 0.000151im, -0.004157 + 0.000266im, -0.005249 + 0.000461im],
)
for n = 1:3
   @printf("n = %d   V2 = %s   V3 = %s   V4 = %s   max |Δ to article| = %.1e\n", n, c(V4h[1, n+1]), c(V4h[2, n+1]), c(V4h[3, n+1]), maximum(abs.(V4h[:, n+1] .- article[n])))
end
println("W^(3): ", join(c.(W4h[:, 4]), "  "), "   (article: 0.003912+0.003727im  0.004853+0.006144im  0.005431+0.008989im)")
````

Partial sum after order 3 (Section 7.8) and the converged value:

````@example workshop_apslf
V3sum = [sum(V4h[i, 1:4]) for i = 1:3]
println("partial sum N = 3: ", join(c.(V3sum), "  "), "   (article: 0.9199-0.1024im  0.8973-0.1293im  0.8765-0.1552im)")
V40 = [sum(V4h[i, :]) for i = 1:3]
println("sum of 40 terms:   ", join(c.(V40), "  "), "   (article: 0.918383-0.102390im  0.895249-0.129200im  0.873961-0.155029im)")
````

And the solver, coefficient by coefficient. `germ = :deviation` is the
article's split; the constant matrix is the series matrix, the shunts are
ramped with $s$ on the right-hand side.

````@example workshop_apslf
Vsol, Vc4, _ = A.apslf_pq(Y4, S4; slack = 1, order = 40, use_pade = true)
for n = 1:3
   @printf("solver n = %d   max |Δ to hand| = %.1e   max |Δ to article| = %.1e\n", n, maximum(abs.(Vc4[:, n+1] .- V4h[:, n+1])), maximum(abs.(Vc4[:, n+1] .- article[n])))
end
Sinj = A.calc_injections(Y4, Vsol)
@printf("solver at s = 1 (Padé): V2 = %s   V3 = %s   V4 = %s\n", c(Vsol[2]), c(Vsol[3]), c(Vsol[4]))
@printf("power mismatch on the physical Y-bus: %.1e pu\n", maximum(abs.(Sinj[red] .- S4[red])))
````

The old recursion without the conjugation produces different coefficients
from order 2 on and ends at a state that does not satisfy the load-flow
equations:

````@example workshop_apslf
V4w, _ = hand_recursion(Yser_red, Ysh[red], S4[red], 40; conjugate = false)
Vw = [1.0 + 0im; [sum(V4w[i, :]) for i = 1:3]]
@printf("old recursion, n = 2: V2 = %s   (article: -0.019635+0.000848im)\n", c(V4w[1, 3]))
@printf("old recursion at s = 1: power mismatch %.2e pu   new recursion: %.1e pu\n", maximum(abs.(A.calc_injections(Y4, Vw)[red] .- S4[red])), maximum(abs.(Sinj[red] .- S4[red])))
````

### 3. The same check with PV buses

The article's examples are PQ-only. The direct PV kernel had a second,
independent error (the sign of the reactive unknown, wrong from order 2
on). Making bus 3 of the four-bus network a PV bus at its converged
magnitude must reproduce the PQ solution exactly, in the dense and in the
sparse kernel.

````@example workshop_apslf
Vm3 = abs(Vsol[3])
bt = [:slack, :pq, :pv, :pq]
P = real.(S4); Q = imag.(S4); Vm = [1.0, 1.0, Vm3, 1.0]
for kern in (A.apslf_pf_pv_direct, A.apslf_pf_pv_direct_sparse)
   Vpv, Qpv, _, _, _ = kern(Y4, bt, P, Q, Vm; slack = 1, order = 40, self_check = false)
   @printf("%-26s max |V - V_pq| = %.1e   Q3 = %.6f (PQ case: %.6f)\n", nameof(kern), maximum(abs.(Vpv .- Vsol)), Qpv[1], Q[3])
end
````

All three checks agree with the article to the printed digits, and the
solver reaches the converged values of Section 7.8 without any Newton
polish.

## Part 2: the 9-bus teaching case

### 4. Solving a case

The solver works on a bus admittance matrix `Y` and a small NamedTuple data
contract: `bustype` (`:slack`, `:pv`, `:pq`), `Pspec`/`Qspec` injections in
pu (load negative), `Vm` voltage setpoints, `Qmin`/`Qmax` for PV buses and
the `slack` index. The 9-bus teaching case ships with the package: two PV
generators, six π-lines with charging, slack at 1.04 pu.

````@example workshop_apslf
case = A.demo_case_9bus()
mismatch(case, res) = maximum(A.compute_demo_mismatch(case, res))
res = solve_pf_apslf(case; mode = :direct, order = 40, use_pade = true, nr_polish = false)
println("converged = $(res.converged), effective mode = $(res.effective_mode), outer iterations = $(res.outer_iters)")
@printf("max mismatch on the physical Y-bus: %.2e pu (no Newton polish involved)\n", mismatch(case, res))
for i = 1:9
   @printf("%-5s %-5s |V| = %.5f pu   angle = %8.4f°\n", case.labels[i], case.bustype[i], abs(res.V[i]), rad2deg(angle(res.V[i])))
end
````

### 5. The series behind the numbers

With `return_coeffs = true` the coefficient matrix is returned. The
coefficients decay geometrically; Padé approximants (Section 5.2) turn the
polynomial into a rational function and extend the reach of the series to
cases where the Taylor sum does not converge at $s=1$.

````@example workshop_apslf
res = solve_pf_apslf(case; order = 40, nr_polish = false, return_coeffs = true)
Vc = res.Vcoeff
for n in (0, 1, 2, 5, 10, 20, 40)
   @printf("n = %2d   max_i |V_i^(n)| = %.2e\n", n, maximum(abs.(Vc[:, n+1])))
end
c5 = Vc[5, :]
taylor = A.evaluate_series(c5, A.APSLFEvaluationOptions(mode = :taylor)).voltage
pade = A.evaluate_series(c5, A.APSLFEvaluationOptions(mode = :pade)).voltage
@printf("bus 5:  Taylor %.8f ∠ %.5f°   Padé %.8f ∠ %.5f°   |Δ| = %.1e\n", abs(taylor), rad2deg(angle(taylor)), abs(pade), rad2deg(angle(pade)), abs(taylor - pade))
````

The poles of the Padé denominator are a heuristic distance-to-collapse
indicator (Section 5.3). Scaling all injections up moves the nearest pole
towards $s = 1$:

````@example workshop_apslf
for factor in (1.0, 1.5, 2.0, 2.5)
   heavy = merge(case, (Pspec = factor .* case.Pspec, Qspec = factor .* case.Qspec, Qmin = fill(-1e9, 9), Qmax = fill(1e9, 9)))
   r = solve_pf_apslf(heavy; order = 40, nr_polish = false, return_coeffs = true)
   st = A.stability_from_Vcoeff(r.Vcoeff; slack = 1, order = 40)
   @printf("load × %.1f: converged = %-5s  min |V| = %.4f  pole distance = %.3f  %s\n", factor, r.converged, minimum(abs.(r.V)), st.dmin, A.st_level(st.dmin))
end
````

### 6. The germ on the 9-bus case

The row sums of the 9-bus matrix are non-zero (line charging), and the
slack sits at 1.04 pu, so the flat germ on the full Y-bus is not an exact
order-0 state (Section 2.2). The `germ` keyword selects the embedding:
`:deviation` (default, Section 6.5 variant 1), `:noload` (variant 2), or
the legacy `:flat`. Both exact variants give the same solution; the legacy
germ only arrives at one with the Newton polish.

````@example workshop_apslf
for germ in (:deviation, :noload, :flat)
   r = solve_pf_apslf(case; order = 40, nr_polish = false, germ = germ)
   @printf("germ = %-10s converged = %-5s  max mismatch = %.2e pu\n", germ, r.converged, mismatch(case, r))
end
r = solve_pf_apslf(case; order = 40, nr_polish = true)
@printf("with polish: score before %.1e, after %.1e, improved = %s\n", r.nr_polish_score_before, r.nr_polish_score_after, r.nr_polish_improved)
````

### 7. PV buses and reactive limits

Two ways to handle PV buses (Section 6): the outer loop that adjusts the
reactive injection until the voltage magnitude is met (`mode = :outer`),
or the direct formulation with one augmented real system per order
(`mode = :direct`). Reactive limits switch PV→PQ at the violated limit.

````@example workshop_apslf
for mode in (:direct, :outer)
   r = solve_pf_apslf(case; mode = mode, order = 40, nr_polish = false)
   sw = get(r, :switch_log, ())
   @printf("mode = %-7s converged = %-5s outer iterations = %d  PV→PQ switches = %d  Q(bus 2) = %.4f pu\n", mode, r.converged, r.outer_iters, length(sw), r.Q[2])
   for e in sw
      @printf("   outer %d: bus %d hit its %s limit (Q = %.4f pu)\n", e.outer, e.bus, e.side, e.qinj)
   end
end
r = solve_pf_apslf(case; order = 40, nr_polish = false, enforce_q_limits = false)
@printf("enforce_q_limits = false: Q(bus 3) = %.4f pu, band [%.2f, %.2f]\n", r.Q[3], case.Qmin[3], case.Qmax[3])
````

### 8. Sparse matrices

For larger networks pass a sparse `Y`; `solve_pf_apslf` then selects the
sparse direct PV kernel. The recursion, the germ and the Padé evaluation
are identical. The synthetic 118-bus case illustrates it.

````@example workshop_apslf
using SparseArrays
case118 = A.demo_case_118bus_synthetic()
sparse_case = merge(case118, (Y = sparse(case118.Y),))
r = solve_pf_apslf(sparse_case; order = 40, nr_polish = false)
t = @elapsed r = solve_pf_apslf(sparse_case; order = 40, nr_polish = false)
@printf("118 buses, nnz(Y) = %d: converged = %s, mode = %s, max mismatch = %.1e pu, %.3f s\n", nnz(sparse_case.Y), r.converged, r.effective_mode, mismatch(case118, r), t)
````

