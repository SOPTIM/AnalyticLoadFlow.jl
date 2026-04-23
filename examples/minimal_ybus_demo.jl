# =============================================================================
# File: examples/minimal_ybus_demo.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: 
# Minimal APSLF integration template using direct Y-bus input.
# Run with: julia --project=. examples/minimal_ybus_demo.jl [--case=9|lv400|118|all] [--inner=pq|direct_pv|both] [--order=40] [--no-nr-polish]
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

using APSLF
using LinearAlgebra
using Printf

# -----------------------------------------------------------------------------
# 1. Data contract comment
# -----------------------------------------------------------------------------
# APSLF applications can pass a small case object as a NamedTuple. This example
# expects the following fields:
#   Y        :: AbstractMatrix{<:Complex}   bus admittance matrix
#   bustype  :: Vector{Symbol}              :slack, :pv, or :pq per bus
#   Pspec    :: Vector{<:Real}              active power injections in pu
#   Qspec    :: Vector{<:Real}              reactive power injections in pu
#   Vm       :: Vector{<:Real}              voltage magnitude setpoints in pu
#   Qmin     :: Vector{<:Real}              PV reactive lower limits in pu
#   Qmax     :: Vector{<:Real}              PV reactive upper limits in pu
#   slack    :: Int                         slack-bus index
# Optional fields used only for reporting here are labels and baseMVA.

# -----------------------------------------------------------------------------
# 2. Presentation helpers
# -----------------------------------------------------------------------------
function bus_type_label(bt::Symbol)
   return bt == :slack ? "SL" : bt == :pv ? "PV" : bt == :pq ? "PQ" : string(bt)
end

function print_case_summary(case)
   println("APSLF minimal demo (Y-bus input)")
   println("buses=$(size(case.Y, 1))   slack=$(case.slack)   baseMVA=$(get(case, :baseMVA, 1.0))")
   return nothing
end

function print_bus_voltages_table(case, bustype_final, V; max_rows = nothing)
   labels = get(case, :labels, ["Bus$(i)" for i in axes(case.Y, 1)])
   nrow = max_rows === nothing ? length(V) : min(max_rows, length(V))
   println("bus  label        type   |V| (pu)   angle (deg)")
   @inbounds for i = 1:nrow
      @printf(
         "%3d  %-10s  %-4s  %8.6f  %10.4f°\n",
         i,
         labels[i],
         bus_type_label(bustype_final[i]),
         abs(V[i]),
         rad2deg(angle(V[i]))
      )
   end
   if nrow < length(V)
      println("... $(length(V) - nrow) additional buses omitted")
   end
   return nothing
end

function print_solution_summary(case, res; inner = :pq, order = 40)
   println("\nConverged   = ", res.converged)
   println("outer_iters = ", get(res, :outer_iters, 1))

   bt_final = get(res, :bustype, case.bustype)
   switched_pv = Int[]
   for i in eachindex(case.bustype)
      if case.bustype[i] == :pv && bt_final[i] == :pq
         push!(switched_pv, i)
      end
   end
   switch_log = get(res, :switch_log, ())
   if !isempty(switched_pv)
      println("\nNote: bus types changed due to Q-limits (PV→PQ switching).")
      println("Q-limit PV→PQ switches: ", length(switched_pv), " (buses: ", join(switched_pv, ", "), ")")
   elseif !isempty(switch_log)
      println("Q-limit PV→PQ switches: ", length(switch_log))
   else
      println("Q-limit PV→PQ switches: 0")
   end

   min_vm = minimum(abs.(res.V))
   max_vm = maximum(abs.(res.V))
   @printf("voltage range: min |V| = %.6f pu, max |V| = %.6f pu\n", min_vm, max_vm)

   println("\n--- Voltages ---")
   max_rows = length(res.V) <= 20 ? nothing : 20
   print_bus_voltages_table(case, bt_final, res.V; max_rows = max_rows)

   maxP, maxQ = APSLF.compute_demo_mismatch(case, res)
   println("\n--- Sanity check: mismatch on given Y-bus ---")
   @printf("Max |ΔP| (pu) = %.3e\n", maxP)
   @printf("Max |ΔQ| (pu) = %.3e\n", maxQ)
   return nothing
end

function print_stability_summary(res; slack, order)
   Vcoeff = get(res, :Vcoeff, nothing)
   println("\n--- APSLF stability (Padé poles) ---")
   if Vcoeff === nothing
      println("Not available (no Vcoeff returned).")
      return nothing
   end

   st = APSLF.stability_from_Vcoeff(Vcoeff; slack = slack, order = order, critical_poles = :auto)
   if isfinite(st.dmin)
      lvl = APSLF.st_level(st.dmin)
      @printf(
         "st_dmin = %.3e   st_lvl = %s   st_bus = %d   pole = %+.6f%+.6fi   [L/M]=[%d/%d]\n",
         st.dmin,
         lvl,
         st.bus,
         real(st.pole),
         imag(st.pole),
         st.L,
         st.M,
      )
      if !isempty(st.critical)
         println("critical poles (nearest to s=1):")
         for cp in st.critical
            @printf("  bus=%d  d=%.3e  pole=%+.6f%+.6fi\n", cp.bus, cp.distance, real(cp.pole), imag(cp.pole))
         end
      end
   else
      println("stability unavailable (no valid pole data).")
   end
   return nothing
end

function run_demo(case = APSLF.demo_case_9bus(); inner::Symbol, order::Int = 40, nr_polish::Bool = true)
   use_pade = true
   verbose = 0
   max_outer = 20

   println("\n", "="^80)
   print_case_summary(case)
   println("inner=$(inner)   order=$(order)   use_pade=$(use_pade)   nr_polish=$(nr_polish)")
   println("="^80)

   res = APSLF.solve_demo_case(
      case;
      inner = inner,
      order = order,
      use_pade = use_pade,
      nr_polish = nr_polish,
      verbose = verbose,
      max_outer = max_outer,
      return_coeffs = true,
   )
   print_solution_summary(case, res; inner = inner, order = order)
   print_stability_summary(res; slack = case.slack, order = order)
   return res
end

# -----------------------------------------------------------------------------
# 3. CLI argument parsing
# -----------------------------------------------------------------------------
function parse_demo_args(args = ARGS)
   opts = (case = "9", inner = :both, order = 40, nr_polish = true)
   for arg in args
      if arg == "--no-nr-polish"
         opts = merge(opts, (nr_polish = false,))
      elseif startswith(arg, "--case=")
         case = last(split(arg, "="; limit = 2))
         if !(case in ("9", "lv400", "118", "all"))
            error("Unsupported case: $(case). Use --case=9, --case=lv400, --case=118, or --case=all.")
         end
         opts = merge(opts, (case = case,))
      elseif startswith(arg, "--inner=")
         inner_text = last(split(arg, "="; limit = 2))
         inner = Symbol(inner_text)
         if !(inner in (:pq, :direct_pv, :both))
            error("Unsupported inner mode: $(inner_text). Use --inner=pq, --inner=direct_pv, or --inner=both.")
         end
         opts = merge(opts, (inner = inner,))
      elseif startswith(arg, "--order=")
         order_text = last(split(arg, "="; limit = 2))
         order = tryparse(Int, order_text)
         if order === nothing || order <= 0
            error("Invalid order: $(order_text). Use a positive integer, for example --order=40.")
         end
         opts = merge(opts, (order = order,))
      else
         error(
            "Unsupported option: $(arg). Supported options: --case=9|lv400|118|all, --inner=pq|direct_pv|both, --order=N, --no-nr-polish.",
         )
      end
   end
   return opts
end

function selected_cases(case_option)
   if case_option == "9"
      return [("9-bus teaching case", APSLF.demo_case_9bus())]
   elseif case_option == "lv400"
      return [("synthetic 400 V LV street-feeder case", APSLF.demo_case_lv_400v_streets())]
   elseif case_option == "118"
      return [("synthetic 118-bus integration case", APSLF.demo_case_118bus_synthetic())]
   elseif case_option == "all"
      return [
         ("9-bus teaching case", APSLF.demo_case_9bus()),
         ("synthetic 400 V LV street-feeder case", APSLF.demo_case_lv_400v_streets()),
         ("synthetic 118-bus integration case", APSLF.demo_case_118bus_synthetic()),
      ]
   else
      error("Unsupported case: $(case_option). Use --case=9, --case=lv400, --case=118, or --case=all.")
   end
end

function selected_inners(case_label, inner_option)
   if case_label in ("synthetic 400 V LV street-feeder case", "synthetic 118-bus integration case")
      if inner_option == :direct_pv
         error("--inner=direct_pv is demonstrated only for --case=9 in this compact example.")
      elseif inner_option == :both
         println("For LV/118 synthetic integration cases, --inner=both runs only :pq to keep the example compact.")
         return (:pq,)
      end
   end
   return inner_option == :both ? (:pq, :direct_pv) : (inner_option,)
end

# -----------------------------------------------------------------------------
# 4. main()
# -----------------------------------------------------------------------------
function main(args = ARGS)
   opts = parse_demo_args(args)
   cases = selected_cases(opts.case)
   results = Dict{Tuple{String,Symbol},Any}()

   for (case_label, case) in cases
      println("\nSelected case: $(case_label)")
      inners = selected_inners(case_label, opts.inner)
      for inner in inners
         label = inner == :pq ? "outer-loop PV handling" : "direct PV handling"
         println("\nRunning PV mode: $(label) (inner = :$(inner))")
         results[(case_label, inner)] = run_demo(case; inner = inner, order = opts.order, nr_polish = opts.nr_polish)
      end
   end

   if opts.inner == :both && opts.case == "9"
      case_label = "9-bus teaching case"
      println("\n", "-"^80)
      println("Comparison of both PV handling methods")
      println("-"^80)
      @printf(
         "outer-loop (:pq):      converged=%s, outer_iters=%d\n",
         string(results[(case_label, :pq)].converged),
         get(results[(case_label, :pq)], :outer_iters, 1)
      )
      @printf(
         "direct PV (:direct_pv): converged=%s, outer_iters=%d\n",
         string(results[(case_label, :direct_pv)].converged),
         get(results[(case_label, :direct_pv)], :outer_iters, 1)
      )
   end
   return nothing
end

if get(ENV, "APSLF_SUITE_NO_AUTORUN", "0") != "1"
   main_fn = getfield(@__MODULE__, :main)
   Base.invokelatest(main_fn)
end
