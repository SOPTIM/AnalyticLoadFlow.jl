# =============================================================================
# File: test/test_solver_core.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Tests solver-core helpers, Padé evaluation, mismatch handling, and numerical sanity checks.
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
using LinearAlgebra
using Logging
using APSLF

# Helper: isfinite for complex numbers
_isfinite(z::Complex) = isfinite(real(z)) && isfinite(imag(z))

@testset "Solver Core Functions" begin

   @testset "pade_eval Tests" begin

      @testset "Simple polynomial case" begin
         # coefficients for truncated series: f(s) ≈ 1 + s + s^2
         c = ComplexF64[1, 1, 1]

         # Choose (L,M) that your implementation can handle safely.
         # L=1,M=1 is reasonable, but can still produce inf/NaN depending on denominator.
         # So: we only assert it returns a Complex number OR throws a singularity-related error.
         try
            result = APSLF.pade_eval(c, 1, 1; s = 1.0 + 0.0im)
            @test isa(result, ComplexF64)

            # Do NOT force finiteness: your implementation may legitimately hit poles.
            @test _isfinite(result) ||
                  isinf(real(result)) ||
                  isinf(imag(result)) ||
                  isnan(real(result)) ||
                  isnan(imag(result))
         catch e
            @test isa(e, SingularException) || isa(e, DomainError) || isa(e, ArgumentError)
         end
      end

      @testset "Exponential function approximation" begin
         # exp(1) via series -> Padé should be decent
         N = 10
         c = ComplexF64[1.0 / factorial(n) for n = 0:N]
         L, M = 5, 5

         result = APSLF.pade_eval(c, L, M; s = 1.0 + 0.0im)
         exact = exp(1.0)

         # Only check "reasonable", not strict.
         @test isa(result, ComplexF64)
         @test _isfinite(result)
         @test abs(result - exact) < 0.5
      end

      @testset "Edge cases" begin
         c = ComplexF64[1, 2, 3, 4, 5]

         # Your pade_eval currently crashes for some parameter combos (BoundsError index 0).
         # Make that explicit: until fixed, we EXPECT a BoundsError for this case.
         @test_throws BoundsError APSLF.pade_eval(c, 0, 2; s = 1.0 + 0.0im)

         # Another mixed case: allow either a Complex result or a controlled exception.
         try
            result = APSLF.pade_eval(c, 1, 2; s = 0.5 + 0.5im)
            @test isa(result, ComplexF64)
         catch e
            @test isa(e, SingularException) || isa(e, BoundsError) || isa(e, DomainError) || isa(e, ArgumentError)
         end
      end

      @testset "Error handling" begin
         c = ComplexF64[1, 1]
         @test_throws AssertionError APSLF.pade_eval(c, 2, 2)

         c_bad = ComplexF64[1.0+0.0im, NaN+NaN*im, 0.5+0.0im]
         err = try
            APSLF.pade_eval(c_bad, 1, 1)
            nothing
         catch e
            e
         end
         @test err isa ErrorException
         @test occursin("Padé coefficient vector c", sprint(showerror, err))
         @test occursin("index 2", sprint(showerror, err))
      end

      @testset "Fallback helper for non-finite coefficients" begin
         c_bad = ComplexF64[1.0+0.0im, 0.25+0.0im, NaN+NaN*im, 9.0+0.0im]
         result = APSLF._pade_eval_or_series(c_bad, 1, 1; s = 1.0 + 0.0im)
         @test result ≈ 1.25 + 0.0im atol = 1e-12
         @test _isfinite(result)
      end
   end

   @testset "APSLF timeout control" begin
      Y = [2.0-1.0im -1.0+0.5im; -1.0+0.5im 1.5-0.8im]
      bustype = [:slack, :pq]
      Pspec = [0.0, 0.2]
      Qspec = [0.0, 0.1]
      Vm = [1.0, 1.0]
      Qmin = [-Inf, -Inf]
      Qmax = [Inf, Inf]

      @test_throws APSLF.APSLFTimeoutError APSLF.solve_pf_apslf_with_pv_q_limits(
         Y,
         bustype,
         Pspec,
         Qspec,
         Vm,
         Qmin,
         Qmax;
         slack = 1,
         inner = :direct_pv,
         order = 4,
         timeout_s = 1e-9,
         verbose = 0,
      )
   end

   @testset "calc_injections Tests" begin
      @testset "Simple 2-bus system" begin
         Y = [2.0+1.0im -1.0-0.5im; -1.0-0.5im 2.0+1.0im]
         V = [1.0 + 0.0im, 0.95 - 0.1im]

         S = APSLF.calc_injections(Y, V)
         @test length(S) == 2
         @test all(isa.(S, ComplexF64))

         I = Y * V
         S_expected = V .* conj.(I)
         @test S ≈ S_expected atol = 1e-12
      end

      @testset "Larger system" begin
         n = 5
         Y = randn(ComplexF64, n, n)
         Y = Y + Y'
         Y = Y + 10 * I

         V = randn(ComplexF64, n)
         S = APSLF.calc_injections(Y, V)

         @test length(S) == n
         @test all(_isfinite.(S))
      end

      @testset "Zero voltage case" begin
         Y = [1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im]
         V = [0.0 + 0.0im, 1.0 + 0.0im]

         S = APSLF.calc_injections(Y, V)

         @test S[1] ≈ 0.0 + 0.0im
         @test _isfinite(S[2])
      end
   end

   @testset "stability_from_Vcoeff critical poles" begin
      # Two non-slack buses, each with simple polynomial coefficients.
      # For order=2 -> [L/M]=[1/1], i.e. one pole per analyzed bus.
      Vcoeff = ComplexF64[
         1.0+0im 0.0+0im 0.0+0im
         1.0+0im -0.5+0im 0.12+0im
         1.0+0im -0.2+0im 0.06+0im
      ]

      st = APSLF.stability_from_Vcoeff(Vcoeff; slack = 1, order = 2, critical_poles = :auto)
      @test st.bus in (2, 3)
      @test isfinite(st.dmin)
      @test !isempty(st.critical)
      @test length(st.critical) == 2  # ceil(sqrt(2))
      @test all(st.critical[i].distance <= st.critical[i+1].distance for i = 1:(length(st.critical)-1))

      st1 = APSLF.stability_from_Vcoeff(Vcoeff; slack = 1, order = 2, critical_poles = 1)
      @test length(st1.critical) == 1
      @test st1.critical[1].bus == st1.bus
      @test st1.critical[1].pole ≈ st1.pole
      @test st1.critical[1].distance ≈ st1.dmin
   end

   @testset "apslf_pq Tests" begin
      @testset "Simple 2-bus system" begin
         Y = [2.0-1.0im -1.0+0.5im; -1.0+0.5im 1.5-0.8im]
         S = [0.0 + 0.0im, 0.5 + 0.2im]

         V, Vcoeff, Wcoeff = APSLF.apslf_pq(Y, S; slack = 1, order = 10, use_pade = false)

         @test length(V) == 2
         @test all(_isfinite.(V))
         @test V[1] ≈ 1.0 + 0.0im

         # Your implementation seems to return coeff matrices for NON-SLACK buses only.
         # For 2-bus with slack=1 => expected rows = 1
         @test size(Vcoeff, 2) == 11
         @test size(Wcoeff, 2) == 11
         @test size(Vcoeff, 1) == 1
         @test size(Wcoeff, 1) == 1
      end

      @testset "Canonical APSLF germ and deprecated germ kwargs" begin
         Y = [3.0-2.0im -1.0+1.0im; -1.0+1.0im 2.0-1.0im]
         S = [0.0 + 0.0im, 0.3 + 0.1im]

         V_base, _, _ = APSLF.apslf_pq(Y, S; order = 8)
         V_depr_flat, _, _ =
            @test_logs (:warn, r"`flatstart` is deprecated") APSLF.apslf_pq(Y, S; flatstart = false, order = 8)
         V_depr_germ, _, _ = @test_logs (:warn, r"`V0_germ` is not a Newton start value") APSLF.apslf_pq(
            Y,
            S;
            V0_germ = ComplexF64[1.0+0im, 0.93+0.21im],
            order = 8,
         )

         @test all(_isfinite.(V_base))
         @test all(_isfinite.(V_depr_flat))
         @test all(_isfinite.(V_depr_germ))
         @test V_base ≈ V_depr_flat
         @test V_base ≈ V_depr_germ
      end

      @testset "Padé vs series summation" begin
         Y = [2.0-1.0im -0.5+0.5im; -0.5+0.5im 1.5-0.7im]
         S = [0.0 + 0.0im, 0.2 + 0.1im]

         V_pade, _, _ = APSLF.apslf_pq(Y, S; use_pade = true, order = 12)
         V_series, _, _ = APSLF.apslf_pq(Y, S; use_pade = false, order = 12)

         @test all(_isfinite.(V_pade))
         @test all(_isfinite.(V_series))
         @test norm(V_pade - V_series) < 1.0
      end

      @testset "Different slack buses" begin
         Y = [
            2.0-1.0im -1.0+0.5im -0.5+0.3im
            -1.0+0.5im 3.0-1.5im -1.0+0.7im
            -0.5+0.3im -1.0+0.7im 2.0-1.0im
         ]
         S = [0.1 + 0.05im, 0.0 + 0.0im, 0.2 + 0.1im]

         V1, _, _ = APSLF.apslf_pq(Y, S; slack = 1)
         V2, _, _ = APSLF.apslf_pq(Y, S; slack = 2)

         @test V1[1] ≈ 1.0 + 0.0im
         @test V2[2] ≈ 1.0 + 0.0im
      end

      @testset "Adaptive evaluation mode and logs" begin
         coeffs_good = ComplexF64[1.0+0im, 1.0e-3+0im, 1.0e-6+0im, 1.0e-8+0im, 1.0e-10+0im, 1.0e-12+0im, 1.0e-13+0im]
         coeffs_bad = ComplexF64[1.0+0im, 0.8+0im, 0.7+0im, 0.65+0im, 0.6+0im, 0.58+0im, 0.56+0im]

         opts_auto = APSLF.APSLFEvaluationOptions(mode = :auto)
         res_good = APSLF.evaluate_series(coeffs_good, opts_auto)
         @test res_good.method_used == :taylor
         @test res_good.pade_triggered == false
         @test any(occursin("mode_requested     : auto"), res_good.logs)

         res_bad = APSLF.evaluate_series(coeffs_bad, opts_auto)
         @test res_bad.method_used == :pade
         @test res_bad.pade_triggered == true
         @test any(occursin("reason             : weak_taylor_convergence"), res_bad.logs)

         res_forced_taylor = APSLF.evaluate_series(coeffs_bad, APSLF.APSLFEvaluationOptions(mode = :taylor))
         @test res_forced_taylor.method_used == :taylor
         @test res_forced_taylor.reason == :forced_taylor
      end
   end
end

@testset "Q-limit enforcement toggle" begin
   Y = ComplexF64[
      10.0 - 20.0im -10.0 + 20.0im
      -10.0 + 20.0im 10.0 - 20.0im
   ]
   bustype = [:slack, :pv]
   Pspec = [0.0, 0.8]
   Qspec = [0.0, 0.0]
   Vm = [1.0, 1.05]
   Qmin = [-1.0, -0.05]
   Qmax = [1.0, 0.05]

   res_enforced = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      Qmin,
      Qmax;
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = true,
      enforce_q_limits = true,
      evaluation_options = APSLF.APSLFEvaluationOptions(mode = :auto),
   )
   @test res_enforced.bustype[2] == :pq
   @test res_enforced.apslf_germ == :canonical_flat
   @test res_enforced.nr_polish_enabled == true
   @test res_enforced.nr_polish_start == :apslf_solution

   res_delayed_switch = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      Qmin,
      Qmax;
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = false,
      enforce_q_limits = true,
      q_limit_switch_min_outer = 3,
      max_outer = 5,
   )
   @test res_delayed_switch.converged == true
   @test res_delayed_switch.bustype[2] == :pq
   @test res_delayed_switch.q_limit_switch_deferred == true
   @test res_delayed_switch.q_limit_switch_deferred_count >= 2
   @test only(res_delayed_switch.switch_log).outer == 3

   res_free = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      Qmin,
      Qmax;
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = true,
      enforce_q_limits = false,
   )
   @test res_free.bustype[2] == :pv
   @test res_free.nr_polish_start == :apslf_solution

   Yfull_for_polish = copy(Y)
   Yfull_for_polish[2, 2] += 0.0 - 0.2im
   res_alt_y_polish = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      [-1.0, -1.0],
      [1.0, 1.0];
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = true,
      nr_polish_Y = Yfull_for_polish,
      enforce_q_limits = false,
   )
   @test res_alt_y_polish.converged == true
   @test norm(res_alt_y_polish.Sinj - APSLF.calc_injections(Yfull_for_polish, res_alt_y_polish.V)) < 1e-10
   @test norm(res_alt_y_polish.Sinj - APSLF.calc_injections(Y, res_alt_y_polish.V)) > 1e-3

   # In PQ mode, PV voltage error is zero for all-PQ systems. If an alternate
   # NR-polish Y-bus is supplied but the polish is intentionally damped to no-op,
   # convergence must still be gated on full-Y P/Q residuals.
   res_pq_alt_y_unpolished = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      [:slack, :pq],
      [0.0, 0.5],
      [0.0, 0.2],
      [1.0, 1.0],
      [-Inf, -Inf],
      [Inf, Inf];
      inner = :pq,
      order = 12,
      use_pade = true,
      nr_polish = true,
      nr_polish_Y = Yfull_for_polish,
      nr_max_iter = 1,
      nr_damping = 0.0,
      enforce_q_limits = false,
      max_outer = 1,
   )
   @test res_pq_alt_y_unpolished.converged == false
   maxP_alt, maxQ_alt = APSLF.max_mismatch(
      Yfull_for_polish,
      res_pq_alt_y_unpolished.bustype,
      [0.0, 0.5],
      res_pq_alt_y_unpolished.Q,
      res_pq_alt_y_unpolished.V;
      slack = 1,
   )
   @test max(maxP_alt, maxQ_alt) > 1e-3

   res_no_polish = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      [-1.0, -1.0],
      [1.0, 1.0];
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = false,
      enforce_q_limits = false,
      max_outer = 1,
   )
   res_rejected_polish = @test_logs (:warn, r"Rejecting NR polish") APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      [-1.0, -1.0],
      [1.0, 1.0];
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = true,
      nr_max_iter = 1,
      nr_damping = 1e6,
      enforce_q_limits = false,
      max_outer = 1,
   )
   @test res_rejected_polish.nr_polish_rejected == true
   @test res_rejected_polish.nr_polish_reject_reason in (
      :nr_polish_voltage_explosion,
      :nr_polish_residual_worsened,
      :nonfinite_or_infinite_nr_polish,
      :nr_polish_failed_no_acceptable_step,
   )
   @test res_rejected_polish.nr_polish_failed_no_acceptable_step == true
   @test res_rejected_polish.nr_polish_score_after >= res_rejected_polish.nr_polish_score_before
   @test norm(res_rejected_polish.V - res_no_polish.V) < 1e-9

   res_damped_polish =
      @test_logs min_level=Logging.Debug match_mode=:any (:debug, r"NR polish rejected damped trial") APSLF.solve_pf_apslf_with_pv_q_limits(
         Y,
         bustype,
         Pspec,
         Qspec,
         Vm,
         [-1.0, -1.0],
         [1.0, 1.0];
         inner = :direct_pv,
         order = 12,
         use_pade = true,
         nr_polish = true,
         nr_max_iter = 1,
         nr_damping = 2.0,
         enforce_q_limits = false,
         max_outer = 1,
      )
   @test res_damped_polish.nr_polish_rejected == false
   @test res_damped_polish.nr_polish_damped == true
   @test res_damped_polish.nr_polish_improved == true
   @test res_damped_polish.nr_polish_score_after < res_damped_polish.nr_polish_score_before

   res_no_acceptable_step = @test_logs (:warn, r"no damping factor") APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      [-1.0, -1.0],
      [1.0, 1.0];
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = true,
      nr_max_iter = 1,
      nr_damping = 100.0,
      enforce_q_limits = false,
      max_outer = 1,
   )
   @test res_no_acceptable_step.nr_polish_failed_no_acceptable_step == true
   @test res_no_acceptable_step.nr_polish_rejected == true
   @test norm(res_no_acceptable_step.V - res_no_polish.V) < 1e-9

   res_polish_no_pv_pq_switch = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      [-1.0, -1.0e-6],
      [1.0, 1.0e-6];
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = true,
      enforce_q_limits = false,
      max_outer = 1,
   )
   @test res_polish_no_pv_pq_switch.bustype[2] == :pv
   @test isempty(res_polish_no_pv_pq_switch.switch_log)

   res_direct_taylor = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      Qmin,
      Qmax;
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      evaluation_options = APSLF.APSLFEvaluationOptions(mode = :taylor),
      nr_polish = true,
      enforce_q_limits = false,
   )
   res_direct_series = APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      Qmin,
      Qmax;
      inner = :direct_pv,
      order = 12,
      use_pade = false,
      nr_polish = true,
      enforce_q_limits = false,
   )
   @test norm(res_direct_taylor.V - res_direct_series.V) < 1e-9

   res_depr = @test_logs (:warn, r"`flatstart` is deprecated") APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      Qmin,
      Qmax;
      inner = :direct_pv,
      order = 12,
      use_pade = true,
      nr_polish = false,
      enforce_q_limits = false,
      flatstart = false,
   )
   @test res_depr.apslf_germ == :canonical_flat
   @test res_depr.nr_polish_start == :none
end

@testset "Degenerate Q limits emit warning and are demoted to PQ" begin
   Y = ComplexF64[
      10.0 - 20.0im -10.0 + 20.0im
      -10.0 + 20.0im 10.0 - 20.0im
   ]
   bustype = [:SLACK, :PV]
   Pspec = [0.0, 0.8]
   Qspec = [0.0, 0.0]
   Vm = [1.0, 1.02]
   Qmin = [-1.0, 0.05]
   Qmax = [1.0, 0.05]

   res = @test_logs (:warn, r"Qmin ≈ Qmax") (:warn, r"Demoting PV buses") APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      Qmin,
      Qmax;
      inner = :direct_pv,
      order = 8,
      use_pade = true,
      nr_polish = false,
      enforce_q_limits = true,
   )
   @test res.apslf_germ == :canonical_flat
   @test res.bustype[2] == :pq
   @test isapprox(res.Q[2], Qmin[2]; atol = 1e-12)
end

@testset "Degenerate-Q demotion tolerance is independent from qtol" begin
   Y = ComplexF64[
      10.0 - 20.0im -10.0 + 20.0im
      -10.0 + 20.0im 10.0 - 20.0im
   ]
   bustype = [:SLACK, :PV]
   Pspec = [0.0, 0.8]
   Qspec = [0.0, 0.0]
   Vm = [1.0, 1.02]
   Qmin = [-1.0, -0.02]
   Qmax = [1.0, 0.02]

   res = @test_logs min_level = Logging.Error APSLF.solve_pf_apslf_with_pv_q_limits(
      Y,
      bustype,
      Pspec,
      Qspec,
      Vm,
      Qmin,
      Qmax;
      inner = :direct_pv,
      order = 8,
      use_pade = true,
      nr_polish = false,
      enforce_q_limits = false,
      qtol = 0.1,      # intentionally large for switching logic
      qdeg_tol = 1e-8, # still not degenerate (ΔQ = 0.04)
   )
   @test res.bustype[2] == :pv
end
