# =============================================================================
# File: test/test_minimal_example.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Tests reusable demo-case helpers without including console example scripts.
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
using AnalyticLoadFlow

@testset "demo helper cases" begin
   case = AnalyticLoadFlow.demo_case_9bus()
   @test size(case.Y, 1) == 9
   @test size(case.Y, 2) == 9
   @test length(case.bustype) == 9
   @test case.slack == 1

   res = redirect_stdout(devnull) do
      AnalyticLoadFlow.solve_demo_case(case; inner = :pq)
   end

   @test all(isfinite, real.(res.V))
   @test all(isfinite, imag.(res.V))

   maxP, maxQ = AnalyticLoadFlow.compute_demo_mismatch(case, res)
   @test isfinite(maxP)
   @test isfinite(maxQ)
   @test max(maxP, maxQ) < 1e-6

   case_lv = AnalyticLoadFlow.demo_case_lv_400v_streets()
   @test size(case_lv.Y) == (16, 16)
   @test length(case_lv.bustype) == 16
   @test case_lv.slack == 1
   @test count(bt -> bt == :slack, case_lv.bustype) == 1
   @test all(case_lv.bustype[i] == :pq for i in eachindex(case_lv.bustype) if i != case_lv.slack)

   res_lv = redirect_stdout(devnull) do
      AnalyticLoadFlow.solve_demo_case(case_lv; inner = :pq, order = 40, nr_polish = true)
   end

   @test all(isfinite, real.(res_lv.V))
   @test all(isfinite, imag.(res_lv.V))
   @test minimum(abs.(res_lv.V)) >= 0.95
   @test maximum(abs.(res_lv.V)) <= 1.01

   maxP_lv, maxQ_lv = AnalyticLoadFlow.compute_demo_mismatch(case_lv, res_lv)
   @test isfinite(maxP_lv)
   @test isfinite(maxQ_lv)
   @test max(maxP_lv, maxQ_lv) < 1e-5

   case118 = AnalyticLoadFlow.demo_case_118bus_synthetic()
   @test size(case118.Y) == (118, 118)
   @test length(case118.bustype) == 118
   @test length(case118.Pspec) == 118
   @test length(case118.Qspec) == 118
   @test length(case118.Vm) == 118
   @test length(case118.Qmin) == 118
   @test length(case118.Qmax) == 118
   @test case118.slack == 1
   @test count(bt -> bt == :slack, case118.bustype) == 1
   @test count(bt -> bt == :pv, case118.bustype) >= 5
   @test count(bt -> bt == :pq, case118.bustype) >= 100

   res118 = redirect_stdout(devnull) do
      AnalyticLoadFlow.solve_demo_case(case118; inner = :pq, order = 40, nr_polish = true)
   end

   @test all(isfinite, real.(res118.V))
   @test all(isfinite, imag.(res118.V))

   maxP118, maxQ118 = AnalyticLoadFlow.compute_demo_mismatch(case118, res118)
   @test isfinite(maxP118)
   @test isfinite(maxQ118)
   @test max(maxP118, maxQ118) < 1e-5
end
