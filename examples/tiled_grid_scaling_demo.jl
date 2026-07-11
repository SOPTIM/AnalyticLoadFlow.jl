# =============================================================================
# File: examples/tiled_grid_scaling_demo.jl
# Date: 2026-07-09
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose:
# Parametric synthetic tiled-grid scaling example with compact timing output.
# Run with: julia --project=. examples/tiled_grid_scaling_demo.jl --buses=100 --samples=3
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
using Printf
using Statistics

const SUPPORTED_OPTIONS = "--buses=N, --order=N, --samples=N, --warmup=N, --aspect=FLOAT, --no-nr-polish"

function parse_positive_int_option(arg::AbstractString, name::AbstractString; minimum::Int, inclusive::Bool)
   value_text = last(split(arg, "="; limit = 2))
   value = tryparse(Int, value_text)
   if value === nothing || (inclusive ? value < minimum : value <= minimum)
      comparator = inclusive ? ">= $(minimum)" : "> $(minimum)"
      error("--$(name) must be $(comparator)")
   end
   return value
end

function parse_positive_float_option(arg::AbstractString, name::AbstractString)
   value_text = last(split(arg, "="; limit = 2))
   value = tryparse(Float64, value_text)
   if value === nothing || value <= 0.0
      error("--$(name) must be > 0")
   end
   return value
end

function parse_demo_args(args = ARGS)
   opts = (buses = 2000, order = 24, samples = 3, warmup = 1, aspect = 1.0, nr_polish = true)
   for arg in args
      if arg == "--no-nr-polish"
         opts = merge(opts, (nr_polish = false,))
      elseif startswith(arg, "--buses=")
         opts = merge(opts, (buses = parse_positive_int_option(arg, "buses"; minimum = 4, inclusive = true),))
      elseif startswith(arg, "--order=")
         opts = merge(opts, (order = parse_positive_int_option(arg, "order"; minimum = 0, inclusive = false),))
      elseif startswith(arg, "--samples=")
         opts = merge(opts, (samples = parse_positive_int_option(arg, "samples"; minimum = 1, inclusive = true),))
      elseif startswith(arg, "--warmup=")
         opts = merge(opts, (warmup = parse_positive_int_option(arg, "warmup"; minimum = 0, inclusive = true),))
      elseif startswith(arg, "--aspect=")
         opts = merge(opts, (aspect = parse_positive_float_option(arg, "aspect"),))
      else
         error("Unsupported option: $(arg). Supported options: $(SUPPORTED_OPTIONS).")
      end
   end
   return opts
end

function solve_case(case, opts)
   return AnalyticLoadFlow.solve_demo_case(
      case;
      inner = :pq,
      order = opts.order,
      use_pade = true,
      nr_polish = opts.nr_polish,
      verbose = 0,
      max_outer = 20,
      return_coeffs = false,
   )
end

function timed_solve_samples(case, opts)
   for _ = 1:opts.warmup
      solve_case(case, opts)
   end

   times_ms = Vector{Float64}(undef, opts.samples)
   results = Vector{Any}(undef, opts.samples)
   for i in eachindex(times_ms)
      t0 = time_ns()
      results[i] = solve_case(case, opts)
      times_ms[i] = (time_ns() - t0) / 1.0e6
   end
   return times_ms, results[end]
end

function main(args = ARGS)
   opts = parse_demo_args(args)
   cfg = (
      aspect_ratio = opts.aspect,
      r = 0.01,
      x = 0.08,
      g = 0.0,
      b = 0.0,
      base_mva = 100.0,
      load_mw_per_right_corner = 1.0,
      load_mvar_per_right_corner = 0.3,
      generation_balance = 0.98,
      vm_slack = 1.0,
      vm_flat = 1.0,
   )

   build_t0 = time_ns()
   case, meta = AnalyticLoadFlow.build_tiled_grid_spec(opts.buses, cfg)
   build_time_ms = (time_ns() - build_t0) / 1.0e6

   solve_times_ms, res = timed_solve_samples(case, opts)
   min_vm = minimum(abs.(res.V))
   max_vm = maximum(abs.(res.V))
   max_p_mismatch, max_q_mismatch = AnalyticLoadFlow.compute_demo_mismatch(case, res)
   max_mismatch = max(max_p_mismatch, max_q_mismatch)

   println("Synthetic tiled-grid APSLF scaling example")
   println("requested_buses = ", meta.requested_max_buses)
   println("actual_buses    = ", meta.actual_buses)
   println("rows x cols     = ", meta.rows, " x ", meta.cols)
   println("branches        = ", meta.branch_count)
   println("order           = ", opts.order)
   println("samples         = ", opts.samples)
   println("warmup          = ", opts.warmup)
   println("nr_polish       = ", opts.nr_polish)
   @printf("build_time_ms   = %.3f\n", build_time_ms)
   @printf(
      "solve_time_ms_min/median/max = %.3f / %.3f / %.3f\n",
      minimum(solve_times_ms),
      median(solve_times_ms),
      maximum(solve_times_ms),
   )
   println("converged       = ", res.converged)
   println("outer_iters     = ", get(res, :outer_iters, 1))
   @printf("voltage_range   = %.6f .. %.6f pu\n", min_vm, max_vm)
   @printf("max_mismatch    = %.3e pu\n", max_mismatch)
   return nothing
end

if get(ENV, "APSLF_SUITE_NO_AUTORUN", "0") != "1"
   main_fn = getfield(@__MODULE__, :main)
   Base.invokelatest(main_fn)
end
