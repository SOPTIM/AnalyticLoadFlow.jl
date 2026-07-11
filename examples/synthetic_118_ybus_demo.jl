# =============================================================================
# File: examples/synthetic_118_ybus_demo.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: 
# Direct console entry point for the synthetic IEEE-118-sized Y-bus demo.
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

function main(args = ARGS)
   if !isempty(args)
      error("synthetic_118_ybus_demo.jl does not accept options; use minimal_ybus_demo.jl for CLI variants.")
   end

   println("Synthetic 118-bus-sized APSLF integration case")
   println("Generated data only: not official IEEE 118 benchmark data.")
   case = AnalyticLoadFlow.demo_case_118bus_synthetic()
   run_demo(case; inner = :pq, order = 40, nr_polish = true)
   return nothing
end

if get(ENV, "APSLF_SUITE_NO_AUTORUN", "0") != "1"
   main_fn = getfield(@__MODULE__, :main)
   Base.invokelatest(main_fn)
end
