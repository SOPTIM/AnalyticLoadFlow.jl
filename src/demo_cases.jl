# =============================================================================
# File: src/demo_cases.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Provides reusable synthetic demo cases and helper wrappers shared by examples and tests.
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

"""
    build_ybus_from_branches(nbus, branches)

Build a simple Y-bus from branch tuples `(i, j, r, x, b_total)`, where
`b_total` is the total line charging susceptance (split equally on both ends).
"""
function build_ybus_from_branches(nbus::Int, branches::Vector{NTuple{5,Float64}})
   Y = zeros(ComplexF64, nbus, nbus)
   for (i_f, j_f, r, x, b_total) in branches
      i = Int(i_f)
      j = Int(j_f)
      y_series = inv(complex(r, x))
      y_shunt_half = 0.5im * b_total
      Y[i, i] += y_series + y_shunt_half
      Y[j, j] += y_series + y_shunt_half
      Y[i, j] -= y_series
      Y[j, i] -= y_series
   end
   return Y
end

include(joinpath(@__DIR__, "..", "data", "synthetic_118_case.jl"))
include(joinpath(@__DIR__, "..", "data", "lv_400v_streets_case.jl"))

"""
    demo_case_9bus()

Self-contained 9-bus demo (inspired by common educational 9-bus layouts),
defined directly from branch parameters to build `Y`.
"""
function demo_case_9bus()
   labels = ["Bus1", "Bus2", "Bus3", "Bus4", "Bus5", "Bus6", "Bus7", "Bus8", "Bus9"]
   baseMVA = 100.0
   slack = 1
   bustype = [:slack, :pv, :pv, :pq, :pq, :pq, :pq, :pq, :pq]

   branches = NTuple{5,Float64}[
      (1.0, 4.0, 0.0000, 0.0576, 0.0000),
      (4.0, 5.0, 0.0170, 0.0920, 0.1580),
      (5.0, 6.0, 0.0390, 0.1700, 0.3580),
      (3.0, 6.0, 0.0000, 0.0586, 0.0000),
      (6.0, 7.0, 0.0119, 0.1008, 0.2090),
      (7.0, 8.0, 0.0085, 0.0720, 0.1490),
      (8.0, 2.0, 0.0000, 0.0625, 0.0000),
      (8.0, 9.0, 0.0320, 0.1610, 0.3060),
      (9.0, 4.0, 0.0100, 0.0850, 0.1760),
   ]
   Y = build_ybus_from_branches(length(labels), branches)

   Pspec = [0.0, 1.35, 0.85, 0.0, -0.90, 0.0, -1.00, 0.0, -1.25]
   Qspec = [0.0, 0.0, 0.0, 0.0, -0.30, 0.0, -0.35, 0.0, -0.50]
   Vm = [1.04, 1.025, 1.025, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]

   # Tight PV Q-limits to demonstrate PV→PQ switching in the output.
   Qmin = [-1e9, -0.10, -0.05, -1e9, -1e9, -1e9, -1e9, -1e9, -1e9]
   Qmax = [1e9, 0.15, 0.10, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9]

   return (
      Y = Y,
      labels = labels,
      baseMVA = baseMVA,
      bustype = bustype,
      slack = slack,
      Pspec = Pspec,
      Qspec = Qspec,
      Vm = Vm,
      Qmin = Qmin,
      Qmax = Qmax,
   )
end

"""
    solve_demo_case(case; inner=:pq, order=40, use_pade=true, nr_polish=true,
                    verbose=0, max_outer=20, return_coeffs=true, germ=:deviation, kwargs...)

Call APSLF for a case NamedTuple that follows the demo data contract.
"""
function solve_demo_case(
   case;
   inner::Symbol = :pq,
   order::Int = 40,
   use_pade::Bool = true,
   nr_polish::Bool = true,
   verbose::Int = 0,
   max_outer::Int = 20,
   return_coeffs::Bool = true,
   germ::Symbol = :deviation,
   kwargs...,
)
   return solve_pf_apslf_with_pv_q_limits(
      Matrix{ComplexF64}(case.Y),
      case.bustype,
      case.Pspec,
      case.Qspec,
      case.Vm,
      case.Qmin,
      case.Qmax;
      slack = case.slack,
      Vslack = ComplexF64(case.Vm[case.slack], 0.0),
      inner = inner,
      max_outer = max_outer,
      order = order,
      use_pade = use_pade,
      nr_polish = nr_polish,
      verbose = verbose,
      return_coeffs = return_coeffs,
      germ = germ,
      kwargs...,
   )
end

function compute_demo_mismatch(case, res)
   Q_final = get(res, :Q, case.Qspec)
   bt_final = get(res, :bustype, case.bustype)
   return max_mismatch(case.Y, bt_final, case.Pspec, Q_final, res.V; slack = case.slack)
end
