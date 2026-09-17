# =============================================================================
# File: test/test_germ_fix.jl
# Date: 2026-09-17
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Tests for the reflected-reciprocal PQ recursion, the sign
#          consistency of the direct PV kernels and the exact order-0 state
#          (germ = :deviation / :noload): pure APSLF without NR polish must be
#          a load-flow solution.
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
using SparseArrays
using AnalyticLoadFlow

const A = AnalyticLoadFlow

_series_3bus() = A.build_ybus_from_branches(3, NTuple{5,Float64}[(1, 2, 0.02, 0.10, 0.0), (2, 3, 0.03, 0.15, 0.0), (1, 3, 0.01, 0.08, 0.0)])

@testset "PQ recursion uses the reflected reciprocal (pure APSLF is exact)" begin
   Y = _series_3bus()
   S = ComplexF64[0, -0.8 - 0.3im, -0.6 - 0.2im]
   for germ in (:flat, :noload, :deviation)
      V, _, _ = A.apslf_pq(Y, S; slack = 1, order = 30, use_pade = true, germ = germ)
      @test maximum(abs.(A.calc_injections(Y, V)[2:3] .- S[2:3])) < 1e-12
   end
   for germ in (:noload, :deviation)
      ws = A.build_apslf_pq_workspace(Y; slack = 1, order = 30, germ = germ)
      V, _, _ = A.apslf_pq_solve!(ws, S)
      @test maximum(abs.(A.calc_injections(Y, V)[2:3] .- S[2:3])) < 1e-12
   end
   ws = A.build_apslf_pq_workspace(Y; slack = 1, order = 20, germ = :deviation)
   @test_throws ArgumentError A.apslf_pq_solve!(ws, S; germ = :noload)
   @test_throws ArgumentError A.apslf_germ(Y, [2, 3], 1, 1.0 + 0im, :nonsense)
end

@testset "direct PV kernels are exact without NR polish (dense and sparse)" begin
   Y = _series_3bus()
   bt = [:slack, :pv, :pq]
   P = [0.0, 0.5, -0.6]
   Q = [0.0, 0.0, -0.2]
   for (Vsl, germ) in ((1.0, :flat), (1.0, :noload), (1.05, :noload), (1.05, :deviation)), kern in (A.apslf_pf_pv_direct, A.apslf_pf_pv_direct_sparse)
      Vm = [Vsl, 1.02, 1.0]
      V, Qpv, _, _, _ = kern(Y, bt, P, Q, Vm; slack = 1, Vslack = ComplexF64(Vsl, 0), order = 30, germ = germ, self_check = false)
      Sinj = A.calc_injections(Y, V)
      @test maximum(abs.(real.(Sinj[2:3]) .- P[2:3])) < 1e-12
      @test abs(imag(Sinj[3]) - Q[3]) < 1e-12
      @test abs(abs(V[2]) - 1.02) < 1e-12
      @test abs(Qpv[1] - imag(Sinj[2])) < 1e-10   # internal Q convention equals the physical injection
   end
end

@testset "9-bus case (shunts, Vslack ≠ 1) is exact without NR polish for both embeddings" begin
   case = A.demo_case_9bus()
   for germ in (:deviation, :noload), Ymat in (case.Y, sparse(case.Y)), inner in (:pq, :direct_pv, :direct_pv_sparse)
      res = A.solve_pf_apslf_with_pv_q_limits(
         Ymat, case.bustype, case.Pspec, case.Qspec, case.Vm, case.Qmin, case.Qmax;
         slack = 1, Vslack = ComplexF64(case.Vm[1], 0), inner = inner, order = 40, use_pade = true,
         nr_polish = false, germ = germ, max_outer = 30, use_sparse = issparse(Ymat),
      )
      maxP, maxQ = A.compute_demo_mismatch(case, res)
      @test res.converged
      @test max(maxP, maxQ) < 1e-8
      @test res.apslf_germ == germ
   end
   res_dev = A.solve_demo_case(case; inner = :direct_pv, nr_polish = false, germ = :deviation)
   res_nl = A.solve_demo_case(case; inner = :direct_pv, nr_polish = false, germ = :noload)
   @test maximum(abs.(res_dev.V .- res_nl.V)) < 1e-8   # same function at s = 1
   # the legacy flat germ on the full Y-bus is not exact here
   res_flat = A.solve_pf_apslf_with_pv_q_limits(
      case.Y, case.bustype, case.Pspec, case.Qspec, case.Vm, case.Qmin, case.Qmax;
      slack = 1, Vslack = ComplexF64(case.Vm[1], 0), inner = :direct_pv, order = 40, nr_polish = false, germ = :flat, max_outer = 5,
   )
   @test first(A.compute_demo_mismatch(case, res_flat)) > 1e-3
   @test A.apslf_row_sums(case.Y) ≈ vec(sum(case.Y, dims = 2))
end
