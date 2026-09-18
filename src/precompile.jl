# =============================================================================
# File: src/precompile.jl
# Purpose: Precompile workload. Runs the solver paths the documentation and
#          the workshop notebook use on the small built-in cases, so that the
#          compiled code lands in the package image at install time instead
#          of being compiled on the first call of every function.
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

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
   # Two-bus network with a capacitor bank (theory Section 7.10)
   y = 1.0 - 4.0im
   Y2 = [y -y; -y y+0.2im]
   S2 = [0.0im, -0.5 - 0.15im]
   # Four-bus π-model network (theory Section 7.1)
   branches4 = NTuple{5,Float64}[
      (1, 2, 0.05, 0.15, 0.06), (1, 3, 0.1, 0.3, 0.04), (2, 3, 0.0667, 0.2, 0.05),
      (2, 4, 0.1, 0.3, 0.04), (3, 4, 0.0833, 0.25, 0.06),
   ]
   S4 = ComplexF64[0, -0.4 - 0.15im, -0.5 - 0.175im, -0.3 - 0.1im]
   bt4 = [:slack, :pq, :pv, :pq]
   Vm4 = [1.0, 1.0, 0.9, 1.0]

   @compile_workload begin
      with_logger(NullLogger()) do
         # PQ-only kernel, dense, all three germs, Taylor and Padé
         for germ in (:deviation, :noload, :flat)
            apslf_pq(Y2, S2; slack = 1, order = 8, use_pade = true, germ = germ)
         end
         apslf_pq(Y2, S2; slack = 1, order = 8, use_pade = false)
         Y4 = build_ybus_from_branches(4, branches4)
         apslf_row_sums(Y4)
         V4, Vc4, _ = apslf_pq(Y4, S4; slack = 1, order = 12, use_pade = true)
         calc_injections(Y4, V4)

         # Direct PV kernels, dense and sparse
         apslf_pf_pv_direct(Y4, bt4, real.(S4), imag.(S4), Vm4; slack = 1, order = 12, self_check = false)
         apslf_pf_pv_direct_sparse(Y4, bt4, real.(S4), imag.(S4), Vm4; slack = 1, order = 12, self_check = false)

         # High-level entry point on the 9-bus case: both PV modes, germs,
         # limits on and off, Newton polish, coefficient return
         case = demo_case_9bus()
         res = solve_pf_apslf(case; mode = :direct, order = 12, use_pade = true, nr_polish = false)
         compute_demo_mismatch(case, res)
         solve_pf_apslf(case; mode = :outer, order = 12, nr_polish = false)
         solve_pf_apslf(case; order = 8, nr_polish = false, use_pade = false)
         solve_pf_apslf(case; order = 12, nr_polish = false, enforce_q_limits = false)
         solve_pf_apslf(case; order = 12, nr_polish = false, germ = :noload)
         solve_pf_apslf(case; order = 12, nr_polish = true)
         rc = solve_pf_apslf(case; order = 12, nr_polish = false, return_coeffs = true)
         evaluate_series(rc.Vcoeff[5, :], APSLFEvaluationOptions(mode = :taylor))
         evaluate_series(rc.Vcoeff[5, :], APSLFEvaluationOptions(mode = :pade))
         st = stability_from_Vcoeff(rc.Vcoeff; slack = 1, order = 12)
         st_level(st.dmin)

         # Sparse path with the sparse direct PV kernel
         sparse_case = merge(case, (Y = sparse(case.Y),))
         solve_pf_apslf(sparse_case; order = 12, nr_polish = false)
         solve_pf_apslf(sparse_case; mode = :outer, order = 12, nr_polish = false, max_outer = 3)
      end
   end
end
