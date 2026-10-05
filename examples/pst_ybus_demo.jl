# =============================================================================
# File: examples/pst_ybus_demo.jl
# Date: 2026-09-17
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose:
# Phase-shifting transformer (PST) in APSLF: Y-bus with a complex tap, the two
# embeddings of theory Section 6.5 (deviation embedding vs. no-load germ), an
# angle sweep that shows how the PST redistributes active power in the 9-bus
# ring, and the outer loop for a regulated PST that meets a branch setpoint.
# Run with: julia --project=. examples/pst_ybus_demo.jl [--shift=DEG] [--target=PU] [--order=40] [--sparse]
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

using AnalyticLoadFlow
using LinearAlgebra
using SparseArrays
using Printf

const A = AnalyticLoadFlow

function parse_args(args = ARGS)
   opts = (shift = 10.0, target = 0.6, order = 40, sparse = false)
   for arg in args
      if startswith(arg, "--shift=")
         opts = merge(opts, (shift = parse(Float64, last(split(arg, "="; limit = 2))),))
      elseif startswith(arg, "--target=")
         opts = merge(opts, (target = parse(Float64, last(split(arg, "="; limit = 2))),))
      elseif startswith(arg, "--order=")
         opts = merge(opts, (order = parse(Int, last(split(arg, "="; limit = 2))),))
      elseif arg == "--sparse"
         opts = merge(opts, (sparse = true,))
      else
         error("Unsupported option: $(arg). Use --shift=DEG, --target=PU, --order=N, --sparse.")
      end
   end
   return opts
end

function section(title)
   println("\n", "="^88)
   println(title)
   println("="^88)
end

function solve_case(case; order, germ = :deviation, nr_polish = false)
   return solve_pf_apslf(case; mode = :direct, order = order, use_pade = true, nr_polish = nr_polish, germ = germ)
end

function main(args = ARGS)
   opts = parse_args(args)
   pst = (4, 5)

   # -------------------------------------------------------------------------
   section("1. Branch model: the PST makes the Y-bus non-symmetric")
   # -------------------------------------------------------------------------
   case = demo_case_9bus_pst(shift_deg = opts.shift, pst = pst, enforce_q_limits = false, sparse_output = opts.sparse)
   i, j = pst
   @printf("PST lumped into branch %d-%d: ratio = %.3f, shift = %.2f°, tap on bus %d\n", i, j, case.pst.ratio, case.pst.shift_deg, i)
   pol(z) = @sprintf("%.4f ∠ %.3f°", abs(z), rad2deg(angle(z)))
   @printf("Y[%d,%d] = %s\nY[%d,%d] = %s   (same magnitude, angles differ by 2φ: the matrix is no longer symmetric)\n", i, j, pol(case.Y[i, j]), j, i, pol(case.Y[j, i]))
   rs = A.apslf_row_sums(Matrix{ComplexF64}(case.Y))
   @printf("row sums Y·1 at buses %d, %d: %s, %s  (≠ 0: line charging plus the PST)\n", i, j, pol(rs[i]), pol(rs[j]))
   @printf("row sums of the nominal series matrix Y0: max |Y0·1| = %.1e\n", maximum(abs.(sum(case.Y0, dims = 2))))

   # -------------------------------------------------------------------------
   section("2. Both embeddings of theory Section 6.5 give the same solution without NR polish")
   # -------------------------------------------------------------------------
   results = Dict{Symbol,Any}()
   for germ in (:deviation, :noload)
      res = solve_case(case; order = opts.order, germ = germ)
      maxP, maxQ = A.compute_demo_mismatch(case, res)
      results[germ] = res
      @printf("germ = %-10s converged = %-5s outer = %d   max |ΔP| = %.2e   max |ΔQ| = %.2e\n", germ, res.converged, res.outer_iters, maxP, maxQ)
   end
   @printf("max |V_deviation - V_noload| = %.2e pu\n", maximum(abs.(results[:deviation].V .- results[:noload].V)))
   res_flat = solve_case(case; order = opts.order, germ = :flat)
   maxP, maxQ = A.compute_demo_mismatch(case, res_flat)
   @printf("germ = :flat (legacy, full Y with flat germ) without NR polish: max |ΔP| = %.2e   max |ΔQ| = %.2e  <- not a solution\n", maxP, maxQ)

   res = results[:deviation]
   println("\nBus voltages (shift = $(opts.shift)°):")
   A.print_bus_voltages(res.V; labels = case.labels)
   println("\nBranch flows:")
   flows = branch_flows(res.V, case.branches; baseMVA = case.baseMVA)
   print_branch_flows(flows; labels = case.labels)

   # -------------------------------------------------------------------------
   section("3. Angle sweep: the PST moves active power between the two ring paths")
   # -------------------------------------------------------------------------
   println("shift°     P_45 MW     P_49 MW     P_56 MW    losses MW   |V| min   |V| max   st_dmin")
   for φ in -20.0:5.0:20.0
      c = demo_case_9bus_pst(shift_deg = φ, pst = pst, enforce_q_limits = false, sparse_output = opts.sparse)
      r = solve_pf_apslf(c; mode = :direct, order = opts.order, nr_polish = false, return_coeffs = true)
      fl = branch_flows(r.V, c.branches; baseMVA = c.baseMVA)
      st = A.stability_from_Vcoeff(r.Vcoeff; slack = c.slack, order = opts.order)
      @printf("%6.1f  %10.3f  %10.3f  %10.3f  %10.4f  %8.4f  %8.4f  %8.3f\n", φ,
         branch_active_power(r.V, c.branches, 4, 5) * c.baseMVA,
         branch_active_power(r.V, c.branches, 4, 9) * c.baseMVA,
         branch_active_power(r.V, c.branches, 5, 6) * c.baseMVA,
         sum(f.Ploss_MW for f in fl), minimum(abs.(r.V)), maximum(abs.(r.V)), st.dmin)
   end
   println("(P_45 + P_49 is the infeed from generator 1 at bus 4 and does not depend on the angle; the PST only changes the split.)")

   # -------------------------------------------------------------------------
   section("4. Regulated PST: outer loop on the angle until P_45 meets the setpoint")
   # -------------------------------------------------------------------------
   println("The angle enters Y through e^{jφ}, not polynomially in s, so it is iterated outside the recursion")
   println("(theory Section 6.5): secant on φ, one full APSLF solve per trial angle.\n")
   reg = solve_pf_pst_regulated(
      φ -> demo_case_9bus_pst(shift_deg = φ, pst = pst, enforce_q_limits = false, sparse_output = opts.sparse),
      4, 5, opts.target; order = opts.order, nr_polish = false, verbose = 1,
   )
   @printf("\nresult: shift = %.4f°   P_45 = %.6f pu (target %.3f)   converged = %s   solves = %d\n", reg.shift_deg, reg.P_pu, opts.target, reg.converged, length(reg.history))
   return nothing
end

if get(ENV, "APSLF_SUITE_NO_AUTORUN", "0") != "1"
   Base.invokelatest(main)
end
