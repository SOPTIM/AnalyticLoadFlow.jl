# =============================================================================
# File: test/test_robustness.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Tests numerical robustness, validation paths, and edge-case solver behavior.
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
using Random
using APSLF

"""
Robustness and edge case tests for APSLF solver
"""

@testset "Edge Cases and Robustness" begin

   @testset "Numerical edge cases" begin

      @testset "Very small system (1 bus)" begin
         Y = reshape([1.0 + 0.5im], 1, 1)
         S = [0.0 + 0.0im]

         V, _, _ = APSLF.apslf_pq(Y, S; slack = 1, order = 4)
         @test V[1] ≈ 1.0 + 0.0im
      end

      @testset "Zero load case" begin
         Y = [2.0-1.0im -1.0+0.5im; -1.0+0.5im 2.0-1.0im]
         S = [0.0 + 0.0im, 0.0 + 0.0im]

         # Either returns a valid solution OR throws a clean SingularException.
         try
            V, _, _ = APSLF.apslf_pq(Y, S; order = 6, use_pade = true)
            @test all(isfinite.(V))
            @test abs(V[1] - 1.0) < 1e-12
         catch e
            @test isa(e, SingularException)
         end
      end

      @testset "Very high load case" begin
         Y = [5.0-2.0im -2.0+1.0im; -2.0+1.0im 5.0-2.0im]
         S = [0.0 + 0.0im, 10.0 + 5.0im]

         try
            V, _, _ = APSLF.apslf_pq(Y, S; order = 20, use_pade = true)
            @test all(isfinite.(V))
         catch e
            @test isa(e, SingularException) || isa(e, DomainError) || isa(e, ArgumentError)
         end
      end

      @testset "Very small admittances" begin
         Y = [1e-6-1e-7im -1e-7+1e-8im; -1e-7+1e-8im 1e-6-1e-7im]
         S = [0.0 + 0.0im, 1e-9 + 1e-10im]

         V, _, _ = APSLF.apslf_pq(Y, S; order = 8)
         @test all(isfinite.(V))
      end

      @testset "Very large admittances" begin
         Y = [1e6-1e5im -1e5+1e4im; -1e5+1e4im 1e6-1e5im]
         S = [0.0 + 0.0im, 1e3 + 1e2im]

         V, _, _ = APSLF.apslf_pq(Y, S; order = 8)
         @test all(isfinite.(V))
      end
   end

   @testset "Padé approximant edge cases" begin

      @testset "Minimum coefficient case" begin
         L, M = 2, 1
         c = [1.0 + 0.0im, 0.5 + 0.0im, 0.25 + 0.0im, 0.125 + 0.0im]
         result = APSLF.pade_eval(c, L, M)
         @test isfinite(result) || isinf(result)
      end

      @testset "High order Padé" begin
         N = 30
         c = [1.0 / factorial(min(n, 20)) + 0.0im for n = 0:N]
         L, M = 15, 15
         try
            result = APSLF.pade_eval(c, L, M)
            @test isfinite(result) || isinf(result)
         catch e
            @test isa(e, SingularException)
         end
      end

      @testset "Complex evaluation points" begin
         c = [1.0 + 0.0im, 1.0 + 0.0im, 0.5 + 0.0im, 0.25 + 0.0im]
         L, M = 2, 1

         for s in [1.0 + 0.0im, 0.5 + 0.5im, 0.0 + 1.0im, -0.5 + 0.3im, 2.0 - 1.0im]
            result = APSLF.pade_eval(c, L, M; s = s)
            @test isfinite(result) || isinf(result)
         end
      end

      @testset "Near-zero denominators" begin
         c = [1.0 + 0.0im, -1.0 + 0.0im, 1.0 + 0.0im, -1.0 + 0.0im, 1.0 + 0.0im]
         L, M = 2, 2

         for s in [1.0 + 1e-10im, 1.0 - 1e-10im, 1.0000001 + 0.0im]
            try
               result = APSLF.pade_eval(c, L, M; s = s)
               @test isfinite(result) || isinf(result)
            catch e
               @test isa(e, SingularException)
            end
         end
      end
   end

   @testset "Memory and performance edge cases" begin

      @testset "Repeated solving" begin
         Y = [2.0-1.0im -1.0+0.5im; -1.0+0.5im 2.0-1.0im]
         S = [0.0 + 0.0im, 0.2 + 0.1im]

         results = Vector{Vector{ComplexF64}}(undef, 20)
         for i = 1:length(results)
            V, _, _ = APSLF.apslf_pq(Y, S; order = 8)
            results[i] = V
         end

         @test all(norm.(results .- Ref(results[1])) .< 1e-12)
      end

      @testset "Large coefficient matrices" begin
         Random.seed!(999)

         n = 20
         Y = randn(ComplexF64, n, n)
         Y = (Y + Y') / 2
         Y = Y + (5 + 2im) * I

         S = zeros(ComplexF64, n)
         S[1] = 0.0 + 0.0im
         for i = 2:n
            S[i] = 0.1 * randn() + 0.05 * randn() * im
         end

         V, Vcoeff, Wcoeff = APSLF.apslf_pq(Y, S; order = 15)

         @test length(V) == n
         @test all(isfinite.(V))

         # Only assert consistency between returned objects.
         @test size(Vcoeff, 2) == 16
         @test size(Wcoeff, 2) == 16
         @test size(Vcoeff, 1) == length(V) || size(Vcoeff, 1) == length(V) - 1
         @test size(Wcoeff, 1) == length(V) || size(Wcoeff, 1) == length(V) - 1
      end
   end
   @testset "Recoverable linear solve errors" begin
      @test APSLF._is_recoverable_linear_solve_error(SingularException(0))
      @test APSLF._is_recoverable_linear_solve_error(DomainError(1.0))
      @test APSLF._is_recoverable_linear_solve_error(ArgumentError("bad"))
      @test !APSLF._is_recoverable_linear_solve_error(BoundsError())
   end
end
