# =============================================================================
# File: data/synthetic_118_case.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: 
#   Defines generated synthetic 118-bus-sized educational data for integration examples.
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
    build_synthetic_118_branches()

Build a connected, lightly meshed synthetic 118-bus branch set. Branch tuples
use the same `(i, j, r, x, b_total)` contract as `build_ybus_from_branches`.
"""
function build_synthetic_118_branches()
   branches = NTuple{5,Float64}[]

   # The backbone guarantees one connected component before extra meshing.
   for i = 1:117
      r = 0.018 + 0.0015 * mod(i, 5)
      x = 0.090 + 0.0040 * mod(i, 7)
      b = 0.006 + 0.0010 * mod(i, 3)
      push!(branches, (float(i), float(i + 1), r, x, b))
   end

   # Tie-lines make the synthetic network lightly meshed without storing data.
   for i = 1:9:100
      j = min(i + 17, 118)
      push!(branches, (float(i), float(j), 0.030, 0.150, 0.008))
   end

   for i = 5:13:105
      j = min(i + 23, 118)
      push!(branches, (float(i), float(j), 0.026, 0.130, 0.007))
   end

   for i = 12:16:108
      j = min(i + 31, 118)
      push!(branches, (float(i), float(j), 0.035, 0.180, 0.006))
   end

   return branches
end

"""
    demo_case_118bus_synthetic()

Return a synthetic 118-bus-sized integration case that follows the same APSLF
NamedTuple data contract as `demo_case_9bus()`. It is generated
algorithmically and is not official IEEE 118 benchmark data.
"""
function demo_case_118bus_synthetic()
   nbus = 118
   labels = ["Bus$(i)" for i = 1:nbus]
   baseMVA = 100.0
   slack = 1

   bustype = fill(:pq, nbus)
   bustype[slack] = :slack
   # Every 10th bus is PV so the case exercises generator-bus/Q-limit handling
   # with a simple deterministic pattern developers can inspect quickly.
   for i = 10:10:110
      bustype[i] = :pv
   end

   Pspec = zeros(Float64, nbus)
   Qspec = zeros(Float64, nbus)
   Vm = ones(Float64, nbus)
   Qmin = fill(-1.0, nbus)
   Qmax = fill(1.0, nbus)

   Vm[slack] = 1.02
   for i in eachindex(bustype)
      if bustype[i] == :pv
         Vm[i] = 1.01
         Pspec[i] = 0.085 + 0.005 * mod(i, 3)
         Qmin[i] = -0.35
         Qmax[i] = 0.45
      elseif bustype[i] == :pq
         Pspec[i] = -(0.0060 + 0.0005 * mod(i, 6))
         Qspec[i] = -(0.0020 + 0.0003 * mod(i, 5))
      end
   end

   branches = build_synthetic_118_branches()
   Y = build_ybus_from_branches(nbus, branches)

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
