# =============================================================================
# File: test/test_precision.jl
# Date: 2026-11-11
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Tests for `precision = :double_double` (double-double arithmetic
#          against BigFloat, the direct PV kernel in double-double, the
#          unchanged Float64 path), the NR polish exit below the tolerance,
#          the single pass in direct mode and the precision-limit diagnosis.
#
# Copyright 2026 SOPTIM AG
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =============================================================================

using Test
using Logging
using LinearAlgebra
using AnalyticLoadFlow

const A = AnalyticLoadFlow

@testset "Double-double arithmetic against BigFloat" begin
   setprecision(BigFloat, 400) do
      big_(x::A.DD) = BigFloat(x)
      big_(z::A.CDD) = Complex(BigFloat(z.re), BigFloat(z.im))
      relerr(x, ref) = Float64(abs(big_(x) - ref) / abs(ref))
      # deterministic operands with a full low word
      dd(k) = A.DD(sin(1.0 + k) * 10.0^(k % 7 - 3)) * A.DD(1.0 / 3.0) + A.DD(cos(2.0 * k) * 1e-17)
      tol = 1e-30
      for k = 1:200
         a, b = dd(k), dd(k + 1000)
         ba, bb = big_(a), big_(b)
         @test abs(a.lo) <= eps(a.hi)
         @test relerr(a + b, ba + bb) < tol
         @test relerr(a - b, ba - bb) < tol
         @test relerr(a * b, ba * bb) < tol
         @test relerr(a / b, ba / bb) < tol
         @test relerr(a * b.hi, ba * b.hi) < tol
         za, zb = A.CDD(a, dd(k + 2000)), A.CDD(b, dd(k + 3000))
         bza, bzb = big_(za), big_(zb)
         @test relerr(za + zb, bza + bzb) < tol
         @test relerr(za * zb, bza * bzb) < 1e-28   # the products may cancel
         @test relerr(za / zb, bza / bzb) < 1e-28
         @test big_(conj(za)) == conj(bza)
      end
      # cancellation: the low words carry the result
      x = A.DD(1.0) / A.DD(3.0)
      y = x - A.DD(x.hi)
      @test Float64(y) == x.lo
      @test relerr(x * A.DD(3.0) - A.DD(1.0) + A.DD(1e-20), big_(x) * 3 - 1 + BigFloat(1e-20)) < 1e-10
      @test Float64(A.DD(0.1) + A.DD(0.2)) == 0.1 + 0.2
      @test ComplexF64(A.CDD(1.5 - 2.0im)) == 1.5 - 2.0im
   end
end

@testset "Double-double LU and Padé" begin
   n = 12
   M = [A.CDD(complex(1.0 / (i + j - 1), i == j ? 0.5 : 0.0)) for i = 1:n, j = 1:n]
   xref = [complex(Float64(i), -0.5 * i) for i = 1:n]
   b = [reduce(+, [M[i, j] * xref[j] for j = 1:n]) for i = 1:n]
   @test A._lu_solve_dd!(copy(M), b)
   @test maximum(abs(ComplexF64(b[i]) - xref[i]) for i = 1:n) < 1e-13
   @test !A._lu_solve_dd!(fill(zero(A.CDD), 2, 2), [A.CDD(1.0), A.CDD(1.0)])

   # geometric series 1/(1 - 0.5 s): [1/1] is exact
   c = ComplexF64[0.5^k for k = 0:8]
   cdd = A.CDD.(c)
   @test ComplexF64(A._pade_value_dd(cdd[1:3], 1, 1)) ≈ 2.0 atol = 1e-15
   @test ComplexF64(A._pade_value_dd(cdd, 8, 0)) == sum(c)
   # exp(s) at s = 1: the double-double Padé is not limited by the Float64 Hankel solve
   ce = [A.CDD(1.0)]
   for k = 1:24
      push!(ce, ce[end] / A.CDD(Float64(k)))
   end
   @test abs(ComplexF64(A._pade_value_dd(ce, 12, 12)) - exp(1.0)) < 1e-15
end

@testset "precision option: double-double direct kernel" begin
   for c in (demo_case_9bus(), demo_case_9bus_pst(), demo_case_118bus_synthetic())
      r64 = solve_pf_apslf(c; order = 40, nr_polish = false, enforce_q_limits = false, return_coeffs = true)
      rdd = solve_pf_apslf(c; order = 40, nr_polish = false, enforce_q_limits = false, return_coeffs = true, precision = :double_double)
      @test r64.precision == :float64
      @test rdd.precision == :double_double
      @test rdd.converged
      @test max(A.compute_demo_mismatch(c, rdd)...) < 1e-12
      @test maximum(abs.(rdd.V .- r64.V)) < 1e-10
      @test maximum(abs.(rdd.Vcoeff .- r64.Vcoeff)) < 1e-9 * maximum(abs.(r64.Vcoeff))

      # Float64 path: the explicit default is the default, bit for bit
      r64x = solve_pf_apslf(c; order = 40, nr_polish = false, enforce_q_limits = false, return_coeffs = true, precision = :float64)
      @test r64x.V == r64.V
      @test r64x.Vcoeff == r64.Vcoeff
   end

   # with Q limits and PV/PQ switching
   c = demo_case_9bus()
   r64 = solve_pf_apslf(c; order = 40, nr_polish = false)
   rdd = solve_pf_apslf(c; order = 40, nr_polish = false, precision = :double_double)
   @test rdd.converged == r64.converged
   @test rdd.bustype == r64.bustype
   @test maximum(abs.(rdd.V .- r64.V)) < 1e-9

   @test_throws ArgumentError solve_pf_apslf(c; order = 24, precision = :quad)
   @test_throws ArgumentError solve_pf_apslf(c; mode = :outer, order = 24, precision = :double_double)
end

@testset "Direct mode without a possible PV switch runs one pass" begin
   c = demo_case_9bus()
   kw = (order = 4, nr_polish = false, outer_fallback_nbus_max = 0, return_coeffs = true)
   # not converged at order 4: one pass with limits off or infinite, max_outer with enforced finite limits
   r_off = solve_pf_apslf(c; kw..., enforce_q_limits = false)
   @test !r_off.converged
   @test r_off.outer_iters == 1
   c_inf = merge(c, (Qmin = fill(-Inf, length(c.Qmin)), Qmax = fill(Inf, length(c.Qmax))))
   r_inf = solve_pf_apslf(c_inf; kw...)
   @test !r_inf.converged
   @test r_inf.outer_iters == 1
   @test r_inf.V == r_off.V
   r_lim = solve_pf_apslf(c; kw..., max_outer = 7)
   @test r_lim.converged || r_lim.outer_iters == 7
   # converged runs are unchanged
   @test solve_pf_apslf(c; order = 40, nr_polish = false, enforce_q_limits = false).outer_iters == 1
end

@testset "NR polish stops below the tolerance" begin
   c = demo_case_9bus()
   r0 = solve_pf_apslf(c; order = 40, nr_polish = false, enforce_q_limits = false)
   r1 = @test_logs min_level = Logging.Warn solve_pf_apslf(c; order = 40, nr_polish = true, enforce_q_limits = false, nr_tol = 1e-6)
   @test r1.nr_polish_score_before < 1e-6
   # the candidate is below nr_tol: no Newton step, no rejection
   @test r1.V == r0.V
   @test !r1.nr_polish_rejected
   @test !r1.nr_polish_failed_no_acceptable_step
   @test r1.nr_polish_success
   @test r1.converged
end

@testset "Diagnosis: Float64 precision limit" begin
   # sqrt(1 + s/a): branch point inside the unit circle, s = 1 needs the
   # continuation; relative noise on the coefficients makes the Padé value stagnate
   function coeffs(N, noise)
      C = zeros(ComplexF64, 3, N + 1)
      C[1, 1] = 1.0
      for (i, a) in ((2, 0.7), (3, 0.75))
         c = 1.0
         for n = 0:N
            C[i, n+1] = c * (1 + noise * cos(n^2 + 3i))
            c *= (0.5 - n) / (n + 1) / a
         end
      end
      return C
   end
   d = A.apslf_convergence_diagnosis(coeffs(40, 1e-7); slack = 1)
   @test d.reason == :precision_limit
   @test abs(d.pole_min) < 1.0
   @test occursin("precision = :double_double", d.message)
   @test occursin("nr_polish", d.message)
   @test occursin("48 or 56", d.message)
   @test !occursin("higher order (e.g.", d.message)

   # the same coefficients from a double-double run: the order is the limit
   ddd = A.apslf_convergence_diagnosis(coeffs(40, 1e-7); slack = 1, precision = :double_double)
   @test ddd.reason == :pade_not_converged
   @test occursin("higher order (e.g. 48)", ddd.message)

   # clean coefficients: no precision limit reported
   @test A.apslf_convergence_diagnosis(coeffs(40, 0.0); slack = 1).reason == :mismatch_above_tol
end

@testset "Float64 stagnates, double-double converges (three-bus case, #25)" begin
   # Strong capacitive injection behind a stiff link: the embedding has a Padé
   # pole inside the unit circle (|s| about 0.2), the coefficients grow with
   # the order, and the stiff link turns the Float64 rounding of the two bus
   # voltages into a mismatch far above the tolerance, at every order.
   z = 0.02 + 0.2im
   ystiff = -1e4im
   y = 1 / z
   Y = ComplexF64[y -y 0; -y y+ystiff -ystiff; 0 -ystiff ystiff]
   c = (Y = Y, bustype = [:Slack, :PQ, :PQ], Pspec = [0.0, 0.25, 0.25], Qspec = [0.0, 3.0, 3.0], Vm = ones(3), Qmin = fill(-Inf, 3), Qmax = fill(Inf, 3), slack = 1)
   mismatch(r) = maximum(A.max_mismatch(c.Y, r.bustype, c.Pspec, r.Q, r.V; slack = 1))
   for order in (32, 40, 48)
      r64 = solve_pf_apslf(c; order = order, nr_polish = false, outer_fallback_nbus_max = 0)
      @test !r64.converged
      @test mismatch(r64) > 1e-8
      @test abs(r64.diagnosis.pole_min) < 1.0
      @test r64.outer_iters == 1
      rdd = solve_pf_apslf(c; order = order, nr_polish = false, outer_fallback_nbus_max = 0, precision = :double_double)
      @test rdd.converged
      @test mismatch(rdd) < 1e-9
      @test maximum(abs.(rdd.V .- r64.V)) < 1e-6
   end
end

const _RTE_6495 = normpath(joinpath(@__DIR__, "..", "data", "_downloaded", "case6495rte.m"))
if isfile(_RTE_6495)
   @testset "case6495rte (local download): double-double reaches the tolerance without polish" begin
      c = matpower_case(_RTE_6495)
      r64 = solve_pf_apslf(c; order = 48, nr_polish = false, enforce_q_limits = false)
      @test !r64.converged
      @test r64.reason == :precision_limit
      rdd = solve_pf_apslf(c; order = 48, nr_polish = false, enforce_q_limits = false, precision = :double_double)
      @test rdd.converged
      @test rdd.effective_mode == :direct
      @test !rdd.nr_polish_enabled
      @test max(A.compute_demo_mismatch(c, rdd)...) < 1e-8
   end
else
   @info "case6495rte test skipped (no data/_downloaded/case6495rte.m); run examples/pegase_matpower_demo.jl --case=case6495rte to download it"
end
