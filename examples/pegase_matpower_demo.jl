# =============================================================================
# File: examples/pegase_matpower_demo.jl
# Date: 2026-09-17
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose:
# Large-network integration example: a PEGASE case (default case2869pegase,
# 2869 buses, 509 PV buses, transformers with ratio and phase shift) read from
# a MATPOWER case file, solved with the sparse direct PV kernel and compared
# with the solved state stored in the file. The case file is downloaded from
# the MATPOWER repository into data/_downloaded/ (git-ignored) when missing.
# Run with: julia --project=. examples/pegase_matpower_demo.jl [--case=case2869pegase|case1354pegase|path.m]
#           [--order=40] [--qlimits] [--polish] [--germ=deviation|noload]
#
# Data: PEGASE cases, CC BY 4.0, C. Josz, S. Fliscounakis, J. Maeght, P. Panciatici,
# "AC Power Flow Data in MATPOWER and QCQP Format: iTesla, RTE Snapshots, and PEGASE",
# https://arxiv.org/abs/1603.01533 (fictitious data, not for operation or planning).
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
using Downloads

const A = AnalyticLoadFlow
const MATPOWER_RAW = "https://raw.githubusercontent.com/MATPOWER/matpower/master/data/"
const DOWNLOAD_DIR = normpath(joinpath(@__DIR__, "..", "data", "_downloaded"))

function parse_args(args = ARGS)
   opts = (case = "case2869pegase", order = 40, qlimits = false, polish = false, germ = :deviation)
   for arg in args
      if startswith(arg, "--case=")
         opts = merge(opts, (case = String(last(split(arg, "="; limit = 2))),))
      elseif startswith(arg, "--order=")
         opts = merge(opts, (order = parse(Int, last(split(arg, "="; limit = 2))),))
      elseif arg == "--qlimits"
         opts = merge(opts, (qlimits = true,))
      elseif arg == "--polish"
         opts = merge(opts, (polish = true,))
      elseif startswith(arg, "--germ=")
         opts = merge(opts, (germ = Symbol(last(split(arg, "="; limit = 2))),))
      else
         error("Unsupported option: $(arg). Use --case=NAME|path.m, --order=N, --qlimits, --polish, --germ=deviation|noload.")
      end
   end
   return opts
end

"""
    case_file(name_or_path) -> path

Return the local path of a MATPOWER case file. A bare case name is looked up in
`data/_downloaded/` and downloaded from the MATPOWER repository when missing.
"""
function case_file(name_or_path::AbstractString)
   isfile(name_or_path) && return String(name_or_path)
   name = endswith(name_or_path, ".m") ? name_or_path : name_or_path * ".m"
   path = joinpath(DOWNLOAD_DIR, name)
   if !isfile(path)
      mkpath(DOWNLOAD_DIR)
      url = MATPOWER_RAW * name
      println("downloading $(url)\n        -> $(path)")
      Downloads.download(url, path)
   end
   return path
end

function main(args = ARGS)
   opts = parse_args(args)
   path = case_file(opts.case)

   println("="^88)
   println("APSLF on a MATPOWER case: ", basename(path))
   println("="^88)
   t_import = @elapsed case = matpower_case(path; verbose = 1)
   conv = case.conventions
   @printf("\nimport: %d buses, %d branches, %d PV buses, %d bus shunts, %.2f s\n", size(case.Y, 1), length(case.branches), count(==(:pv), case.bustype), length(case.bus_shunts), t_import)
   npst = count(A.is_phase_shifter, case.branches)
   ntr = count(A.is_transformer, case.branches)
   @printf("transformers: %d with off-nominal ratio or shift, %d of them phase shifters\n", ntr, npst)
   @printf("conventions chosen from the stored solution: angle in %s, sign %+d, ratio %s (max mismatch of the stored state %.2e pu, L1 %.2e pu)\n",
      conv.angle_unit, conv.angle_sign, conv.ratio_convention, conv.ref_mismatch_pu, conv.ref_mismatch_l1_pu)
   if npst > 0
      shifts = [b.shift_deg for b in case.branches if A.is_phase_shifter(b)]
      @printf("phase-shift angles after conversion: %d values in [%.2f°, %.2f°]\n", length(shifts), minimum(shifts), maximum(shifts))
   end

   println("\nsolver: mode = :direct (sparse direct PV kernel), order = $(opts.order), Padé, germ = :$(opts.germ), NR polish = $(opts.polish), Q limits = $(opts.qlimits)")
   # first call includes compilation; time the second one
   solve() = solve_pf_apslf(case; mode = :direct, order = opts.order, use_pade = true, nr_polish = opts.polish, enforce_q_limits = opts.qlimits, germ = opts.germ, return_coeffs = true)
   res = solve()
   t_solve = @elapsed res = solve()
   maxP, maxQ = A.compute_demo_mismatch(case, res)
   @printf("\nconverged = %s   effective mode = %s   outer iterations = %d   time = %.3f s (second call)\n", res.converged, res.effective_mode, res.outer_iters, t_solve)
   @printf("max |ΔP| = %.2e pu   max |ΔQ| (PQ buses) = %.2e pu   on the physical Y-bus\n", maxP, maxQ)
   @printf("|V| range = [%.4f, %.4f] pu\n", minimum(abs.(res.V)), maximum(abs.(res.V)))
   nsw = length(get(res, :switch_log, ()))
   opts.qlimits && @printf("PV→PQ switches due to Q limits: %d\n", nsw)

   dV = abs.(res.V .- case.V_ref)
   k = argmax(dV)
   @printf("\ncomparison with the solved state stored in the case file:\n  max |V - V_ref| = %.2e pu at bus %s, mean = %.2e pu\n", dV[k], case.labels[k], sum(dV) / length(dV))
   println("  (the stored state itself violates the equations by up to $(round(conv.ref_mismatch_pu; sigdigits = 2)) pu at the PST buses, so agreement to that level is expected)")

   Vcoeff = get(res, :Vcoeff, nothing)
   if Vcoeff !== nothing
      st = A.stability_from_Vcoeff(Vcoeff; slack = case.slack, order = opts.order)
      nrm = [maximum(abs.(Vcoeff[:, n+1])) for n = 0:opts.order]
      println("\nseries diagnostics:")
      @printf("  max |V^(n)| for n = 1, 5, 10, 20, %d: %.1e, %.1e, %.1e, %.1e, %.1e\n", opts.order, nrm[2], nrm[6], nrm[11], nrm[21], nrm[end])
      @printf("  nearest Padé pole to s = 1: distance %.3f at bus %s (%s)\n", st.dmin, case.labels[st.bus], A.st_level(st.dmin))
   end
   return nothing
end

if get(ENV, "APSLF_SUITE_NO_AUTORUN", "0") != "1"
   Base.invokelatest(main)
end
