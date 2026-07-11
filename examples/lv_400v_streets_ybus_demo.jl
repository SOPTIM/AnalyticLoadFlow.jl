# =============================================================================
# File: examples/lv_400v_streets_ybus_demo.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: 
# Direct console entry point for the synthetic 400 V low-voltage street-feeder demo.
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

old_autorun = get(ENV, "APSLF_SUITE_NO_AUTORUN", nothing)
ENV["APSLF_SUITE_NO_AUTORUN"] = "1"
try
   include(joinpath(@__DIR__, "minimal_ybus_demo.jl"))
finally
   if old_autorun === nothing
      pop!(ENV, "APSLF_SUITE_NO_AUTORUN", nothing)
   else
      ENV["APSLF_SUITE_NO_AUTORUN"] = old_autorun
   end
end

function print_lv_network_diagram(case)
   labels = get(case, :labels, ["Bus$(i)" for i in axes(case.Y, 1)])
   lbl(i) = "$(labels[i])[$(i)]"

   println()
   println("--- Synthetic LV network topology ---")
   println("All non-slack buses are PQ loads.")
   println()
   println("$(lbl(1))  transformer / slack bus")
   println("├─ Street A: $(lbl(2)) ─ $(lbl(3)) ─ $(lbl(4)) ─ $(lbl(5))")
   println("│                         └─ $(lbl(14))")
   println("├─ Street B: $(lbl(6)) ─ $(lbl(7)) ─ $(lbl(8)) ─ $(lbl(9))")
   println("│                         └─ $(lbl(15))")
   println("└─ Street C: $(lbl(10)) ─ $(lbl(11)) ─ $(lbl(12)) ─ $(lbl(13))")
   println("                          └─ $(lbl(16))")
   println()
   return nothing
end

function main(args = ARGS)
   if !isempty(args)
      error("lv_400v_streets_ybus_demo.jl does not accept options; use minimal_ybus_demo.jl for CLI variants.")
   end

   println("Synthetic 400 V low-voltage street-feeder APSLF example")
   println("Generated data only: not a real grid model.")

   case = AnalyticLoadFlow.demo_case_lv_400v_streets()
   print_lv_network_diagram(case)

   run_demo(case; inner = :pq, order = 40, nr_polish = true)
   return nothing
end

if get(ENV, "APSLF_SUITE_NO_AUTORUN", "0") != "1"
   main_fn = getfield(@__MODULE__, :main)
   Base.invokelatest(main_fn)
end
