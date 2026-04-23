# =============================================================================
# File: data/lv_400v_streets_case.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: 
#   Defines generated synthetic 400 V low-voltage radial street-feeder data for educational examples.
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


# Synthetic 400 V low-voltage radial street-feeder case.
# Parameters are compact educational data, not a real feeder model. Cable
# impedances use plausible low-voltage values converted to per-unit on
# baseMVA = 0.4 and nominal 0.4 kV line-to-line, where Zbase = 0.4 ohm.

"""
    build_lv_400v_streets_branches()

Build branch tuples `(i, j, r, x, b_total)` for a synthetic 400 V radial
street-feeder example. Resistances/reactances are already in per-unit.
"""
function build_lv_400v_streets_branches()
   zbase_ohm = 0.4
   r_ohm_per_km = 0.32
   x_ohm_per_km = 0.08

   cable_pu(length_km) = (r_ohm_per_km * length_km / zbase_ohm, x_ohm_per_km * length_km / zbase_ohm)

   sections = (
      (1, 2, 0.045),
      (2, 3, 0.050),
      (3, 4, 0.055),
      (4, 5, 0.050),
      (1, 6, 0.040),
      (6, 7, 0.050),
      (7, 8, 0.055),
      (8, 9, 0.045),
      (1, 10, 0.050),
      (10, 11, 0.050),
      (11, 12, 0.060),
      (12, 13, 0.050),
      (3, 14, 0.030),
      (7, 15, 0.030),
      (11, 16, 0.035),
   )

   branches = NTuple{5,Float64}[]
   for (i, j, length_km) in sections
      r, x = cable_pu(length_km)
      push!(branches, (float(i), float(j), r, x, 0.0))
   end
   return branches
end

"""
    demo_case_lv_400v_streets()

Return a compact synthetic 400 V low-voltage radial street-feeder case. Bus 1
is the transformer low-voltage/slack bus, and every other bus is a PQ load.
Negative `Pspec`/`Qspec` values are loads, matching the other demo cases.
"""
function demo_case_lv_400v_streets()
   nbus = 16
   labels = [
      "LV-Trafo",
      "A-01",
      "A-02",
      "A-03",
      "A-04",
      "B-01",
      "B-02",
      "B-03",
      "B-04",
      "C-01",
      "C-02",
      "C-03",
      "C-04",
      "A-02a",
      "B-02a",
      "C-02a",
   ]
   baseMVA = 0.4
   slack = 1

   bustype = fill(:pq, nbus)
   bustype[slack] = :slack

   # Small residential/commercial loads in kW/kvar on a 0.4 MVA base.
   loads_kw = [0.0, 10.0, 12.0, 9.0, 8.0, 11.0, 13.0, 9.0, 7.0, 12.0, 14.0, 10.0, 8.0, 5.0, 6.0, 5.0]
   loads_kvar = [0.0, 3.5, 4.0, 3.0, 2.6, 3.8, 4.4, 3.0, 2.3, 4.0, 4.8, 3.4, 2.7, 1.7, 2.0, 1.7]

   Pspec = -loads_kw ./ (baseMVA * 1000.0)
   Qspec = -loads_kvar ./ (baseMVA * 1000.0)
   Pspec[slack] = 0.0
   Qspec[slack] = 0.0

   Vm = ones(Float64, nbus)
   Qmin = fill(-1.0e9, nbus)
   Qmax = fill(1.0e9, nbus)

   branches = build_lv_400v_streets_branches()
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
