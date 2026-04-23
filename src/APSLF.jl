# =============================================================================
# File: src/APSLF.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Defines the APSLF module, includes implementation files, and exports the supported reference API.
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

module APSLF

using LinearAlgebra
using Printf
using Logging
using SparseArrays

include("solver_core.jl")
include("yamlparams.jl")
include("utils.jl")
include("line_flows.jl")
include("demo_cases.jl")

# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------

# NOTE: Keep legacy naming temporarily for backwards compatibility.
# Prefer `solve_pf_series` going forward.

export
   # High-level APSLF interface
   solve_pf_series,
   APSLFEvaluationOptions,
   APSLFTimeoutError,
   evaluate_series,
   is_taylor_good,

   # Legacy alias (deprecated)
   solve_pf_apslf,
   # Line flows / losses
   line_flows_pi,
   total_line_losses,

   # Printing / formatting utilities
   polar_str,
   print_bus_voltages,
   print_bus_voltages_kv,
   print_line_flows,
   print_total_line_losses,
   safe_get,
   st_level,
   with_silent,
   diff_voltages,
   max_mismatch_on_specY,
   max_mismatch,
   any_pv_at_qlimit,
   choose_tiled_grid_dimensions,
   tiled_grid_bus_index,
   build_tiled_grid_ybus,
   build_tiled_grid_spec,
   build_ybus_from_branches,
   demo_case_9bus,
   build_synthetic_118_branches,
   demo_case_118bus_synthetic,
   demo_case_lv_400v_streets,
   solve_demo_case,
   compute_demo_mismatch

end
