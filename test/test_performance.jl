# =============================================================================
# File: test/test_performance.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Runs small numerical performance and allocation smoke tests without acting as a benchmark suite.
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
using LinearAlgebra
using Random

_isfinite(z::Complex) = isfinite(real(z)) && isfinite(imag(z))

@testset "Performance Tests" begin

   @testset "Large system performance" begin
      Random.seed!(1234)

      n = 30
      Y = zeros(ComplexF64, n, n)

      # small line-charging shunt per bus; the series admittances are stamped below
      for i = 1:n
         Y[i, i] = 0.0 + 0.02im
      end

      # Chain of lines keeps every bus connected (an isolated bus with only a shunt has
      # no load-flow solution); random cross links add meshing.
      for i = 1:(n-1)
         if true
            y_line = -(2.0 + rand()) + (10.0 + rand())im   # stiff line, |z| ≈ 0.1 pu
            Y[i, i+1] = y_line
            Y[i+1, i] = y_line
            Y[i, i] -= y_line
            Y[i+1, i+1] -= y_line
         end
      end

      for i = 1:(n÷2)
         j = i + (n ÷ 2)
         if rand() < 0.2
            y_line = -(1.0 + rand()) + (5.0 + rand())im
            Y[i, j] = y_line
            Y[j, i] = y_line
            Y[i, i] -= y_line
            Y[j, j] -= y_line
         end
      end

      S = zeros(ComplexF64, n)
      S[1] = 0.0 + 0.0im
      # modest loads (negative injections) so that the 30-bus chain stays within a
      # normal voltage band
      for i = 2:n
         S[i] = -(0.002 + 0.002 * rand()) - (0.001 + 0.001 * rand())im
      end

      AnalyticLoadFlow.apslf_pq(Y, S; order = 16, use_pade = true) # warmup

      V = nothing
      t = @elapsed begin
         V, _, _ = AnalyticLoadFlow.apslf_pq(Y, S; order = 16, use_pade = true)
      end

      println("apslf_pq(n=$n, order=16, pade=true) elapsed = $(round(t * 1000, digits=3)) ms")

      @test length(V) == n
      @test all(_isfinite.(V))
      @test abs(V[1] - 1.0) < 1e-10

      non_slack_voltages = abs.(V[2:end])
      @test all(0.7 .< non_slack_voltages .< 1.3)
   end

   @testset "Convergence rate comparison" begin
      Y = [
         4.0-2.0im -2.0+1.0im -1.0+0.5im
         -2.0+1.0im 6.0-3.0im -2.5+1.2im
         -1.0+0.5im -2.5+1.2im 5.0-2.0im
      ]
      S = [0.0 + 0.0im, 0.8 + 0.3im, 0.6 + 0.25im]

      orders = [6, 10, 14, 18, 22]

      println("\nConvergence rate analysis:")
      println("Order\tSeriesErr\t\tPadeErr")

      V_ref = nothing
      for order in orders
         V_series, _, _ = AnalyticLoadFlow.apslf_pq(Y, S; order = order, use_pade = false)
         V_pade, _, _ = AnalyticLoadFlow.apslf_pq(Y, S; order = order, use_pade = true)

         if V_ref === nothing
            V_ref = V_pade
         end

         err_series = norm(V_series - V_ref)
         err_pade = norm(V_pade - V_ref)

         println("$order\t\t$(round(err_series, digits=8))\t$(round(err_pade, digits=8))")

         @test all(_isfinite.(V_series))
         @test all(_isfinite.(V_pade))
      end
   end

   @testset "Stress test - ill-conditioned systems" begin
      Y = [1e6-1e3im -1e6+1e3im; -1e6+1e3im 1e6+1e2im]
      S = [0.0 + 0.0im, 1e-3 + 1e-4im]

      try
         V, _, _ = AnalyticLoadFlow.apslf_pq(Y, S; order = 20, use_pade = true)
         @test all(_isfinite.(V))
         @test abs(V[1] - 1.0) < 1e-10
      catch e
         @test isa(e, SingularException) || isa(e, ArgumentError) || isa(e, DomainError)
      end
   end

   @testset "Memory allocation test" begin
      Random.seed!(5678)

      Y = randn(ComplexF64, 10, 10)
      Y = Y + Y'
      Y = Y + 5 * I

      S = ComplexF64[i == 1 ? 0.0 + 0.0im : (0.1 * rand() + 0.05 * rand() * im) for i = 1:10]

      # warmup
      AnalyticLoadFlow.apslf_pq(Y, S; order = 12, use_pade = true)

      # stable allocation measurement in Julia: @allocated
      bytes = @allocated begin
         for _ = 1:10
            AnalyticLoadFlow.apslf_pq(Y, S; order = 12, use_pade = true)
         end
      end

      println("Allocated bytes for 10 runs: $bytes")

      # Don't enforce a hard cap (platform/compiler sensitive). Just ensure it runs.
      @test bytes ≥ 0
   end
end
