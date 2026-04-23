# =============================================================================
# File: src/line_flows.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Computes branch power flows and aggregate line-loss summaries for solved voltage states.
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
    line_flows_pi(V, edges, y_series_pu, y_sh_total_pu; Sbase_MVA=100.0)

Calculates power flows on lines of a network with identical Π-lines.
- V: Bus voltages in pu (ComplexF64)
- edges: Vector{Tuple{Int,Int}} of lines (i,j)
- y_series_pu: Series admittance per line in pu (ComplexF64)
- y_sh_total_pu: Total shunt admittance per line in pu (ComplexF64) (not halved!)
- Sbase_MVA: Base power (for MW/MVAr output)

Returns: Vector of NamedTuples per line:
(i, j, Sij_pu, Sji_pu, Sloss_pu, Pij_MW, Qij_MVAr, ...)
"""
function line_flows_pi(
   V::Vector{ComplexF64},
   edges::Vector{Tuple{Int,Int}},
   y_series_pu::ComplexF64,
   y_sh_total_pu::ComplexF64;
   Sbase_MVA::Float64 = 100.0,
)
   ysh_half = 0.5 * y_sh_total_pu
   flows = Vector{NamedTuple}(undef, length(edges))

   for (k, (i, j)) in enumerate(edges)
      Vi = V[i]
      Vj = V[j]

      Iij = y_series_pu * (Vi - Vj) + ysh_half * Vi
      Iji = y_series_pu * (Vj - Vi) + ysh_half * Vj

      Sij = Vi * conj(Iij)
      Sji = Vj * conj(Iji)
      Sloss = Sij + Sji

      flows[k] = (
         i = i,
         j = j,
         Sij_pu = Sij,
         Sji_pu = Sji,
         Sloss_pu = Sloss,
         Pij_MW = real(Sij) * Sbase_MVA,
         Qij_MVAr = imag(Sij) * Sbase_MVA,
         Pji_MW = real(Sji) * Sbase_MVA,
         Qji_MVAr = imag(Sji) * Sbase_MVA,
         Ploss_MW = real(Sloss) * Sbase_MVA,
         Qloss_MVAr = imag(Sloss) * Sbase_MVA,
      )
   end
   return flows
end


"""
    total_line_losses(flows)

Sums the net power losses over all lines (each line counted once).
`flows` comes from `line_flows_pi`.

Returns:
  (Ploss_MW, Qloss_MVAr, Sloss_pu)
"""
function total_line_losses(flows)
   Ploss_MW = 0.0
   Qloss_MVAr = 0.0
   Sloss_pu = 0.0 + 0.0im

   for f in flows
      Ploss_MW += f.Ploss_MW
      Qloss_MVAr += f.Qloss_MVAr
      Sloss_pu += f.Sloss_pu
   end

   return (Ploss_MW = Ploss_MW, Qloss_MVAr = Qloss_MVAr, Sloss_pu = Sloss_pu)
end
