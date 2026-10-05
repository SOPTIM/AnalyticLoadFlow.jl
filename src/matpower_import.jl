# =============================================================================
# File: src/matpower_import.jl
# Date: 2026-09-17
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Minimal MATPOWER case-file (.m) reader for large integration cases
#          such as the PEGASE networks, with the MATPOWER transformer
#          conventions (an experimental detection from the stored state is
#          available).
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
    parse_matpower_m(path) -> (name, baseMVA, bus, gen, branch)

Read the `mpc.baseMVA`, `mpc.bus`, `mpc.gen` and `mpc.branch` blocks of a
MATPOWER case file. Rows are returned as `Matrix{Float64}` in file order;
comments (`%`) and trailing semicolons are ignored. Nothing else of the
MATLAB syntax is interpreted.
"""
function parse_matpower_m(path::AbstractString)
   isfile(path) || throw(ArgumentError("MATPOWER case file not found: $(path)"))
   name = first(splitext(basename(path)))
   baseMVA = 100.0
   blocks = Dict{String,Vector{Vector{Float64}}}()
   current = nothing
   for raw in eachline(path)
      line = strip(first(split(raw, '%'; limit = 2)))
      isempty(line) && continue
      if current === nothing
         m = match(r"^mpc\.baseMVA\s*=\s*([0-9.eE+-]+)", line)
         if m !== nothing
            baseMVA = parse(Float64, m.captures[1])
            continue
         end
         m = match(r"^mpc\.(bus|gen|branch)\s*=\s*\[(.*)$", line)
         if m !== nothing
            current = m.captures[1]
            blocks[current] = Vector{Vector{Float64}}()
            line = strip(m.captures[2])
            isempty(line) && continue
         else
            continue
         end
      end
      if startswith(line, "]")
         current = nothing
         continue
      end
      row = strip(replace(line, ";" => " "))
      isempty(row) && continue
      vals = Float64[]
      for tok in split(row)
         v = tryparse(Float64, tok)
         v === nothing && break
         push!(vals, v)
      end
      isempty(vals) || push!(blocks[current], vals)
      endswith(line, "];") && (current = nothing)
   end
   for key in ("bus", "gen", "branch")
      haskey(blocks, key) || throw(ArgumentError("block mpc.$(key) not found in $(path)"))
   end
   tomat(rows) = begin
      ncol = minimum(length, rows)
      M = Matrix{Float64}(undef, length(rows), ncol)
      for (r, row) in enumerate(rows)
         M[r, :] .= row[1:ncol]
      end
      M
   end
   return (name = name, baseMVA = baseMVA, bus = tomat(blocks["bus"]), gen = tomat(blocks["gen"]), branch = tomat(blocks["branch"]))
end

"""
    matpower_case(path; angle_unit = :deg, angle_sign = 1, ratio_convention = :matpower,
                  sparse_output = true, verbose = 0)

Build the APSLF case NamedTuple (`Y`, `bustype`, `Pspec`, `Qspec`, `Vm`,
`Qmin`, `Qmax`, `slack`, `labels`, `baseMVA`, `branches`, `bus_shunts`,
`V_ref`, `conventions`) from a MATPOWER case file.

The file is read with the MATPOWER conventions: branch `angle` in degrees,
`ratio` as the tap on the from side. `V_ref` is the state `(Vm, Va)` stored in
the bus block; `conventions.ref_mismatch_pu` is its power mismatch on the
imported Y-bus. A stored state need not be a solution of the MATPOWER model
(the PEGASE files are not).

Experimental: `angle_unit = :auto`, `angle_sign = 0` and
`ratio_convention = :auto` try every combination (degrees/radians × sign ×
ratio/1/ratio) and keep the one whose Y-bus reproduces the stored state with
the smallest total (L1) power mismatch; all trials are returned in
`conventions.trials`. Because the stored state need not match the model, this
can pick a convention other than MATPOWER's.

Buses of MATPOWER type 4 (isolated) and out-of-service branches/generators are
dropped. Generators at a bus are aggregated (P, Q, Qmin, Qmax summed, `Vm`
taken from the first in-service unit). A PV bus without an in-service
generator becomes PQ.
"""
function matpower_case(
   path::AbstractString;
   angle_unit::Symbol = :deg,
   angle_sign::Int = 1,
   ratio_convention::Symbol = :matpower,
   sparse_output::Bool = true,
   verbose::Int = 0,
)
   angle_unit in (:auto, :deg, :rad) || throw(ArgumentError("angle_unit must be :auto, :deg or :rad"))
   angle_sign in (0, 1, -1) || throw(ArgumentError("angle_sign must be 0 (auto), 1 or -1"))
   ratio_convention in (:auto, :matpower, :inverse) || throw(ArgumentError("ratio_convention must be :auto, :matpower or :inverse"))
   mp = parse_matpower_m(path)
   baseMVA = mp.baseMVA
   bus = mp.bus
   gen = mp.gen
   br = mp.branch

   keep = findall(r -> bus[r, 2] != 4.0, 1:size(bus, 1))
   nbus = length(keep)
   index_of = Dict{Int,Int}()
   for (k, r) in enumerate(keep)
      index_of[Int(bus[r, 1])] = k
   end

   bustype = Vector{Symbol}(undef, nbus)
   Pspec = zeros(nbus)
   Qspec = zeros(nbus)
   Vm = ones(nbus)
   Qmin = fill(-1e9, nbus)
   Qmax = fill(1e9, nbus)
   labels = Vector{String}(undef, nbus)
   V_ref = Vector{ComplexF64}(undef, nbus)
   bus_shunts = NamedTuple{(:bus, :g, :b),Tuple{Int,Float64,Float64}}[]
   slack = 0
   for (k, r) in enumerate(keep)
      t = Int(bus[r, 2])
      bustype[k] = t == 3 ? :slack : (t == 2 ? :pv : :pq)
      t == 3 && slack == 0 && (slack = k)
      Pspec[k] = -bus[r, 3] / baseMVA
      Qspec[k] = -bus[r, 4] / baseMVA
      Vm[k] = bus[r, 8]
      labels[k] = string(Int(bus[r, 1]))
      V_ref[k] = bus[r, 8] * cis(deg2rad(bus[r, 9]))
      (bus[r, 5] != 0.0 || bus[r, 6] != 0.0) && push!(bus_shunts, (bus = k, g = bus[r, 5] / baseMVA, b = bus[r, 6] / baseMVA))
   end
   slack == 0 && throw(ArgumentError("no slack bus (type 3) in $(path)"))

   has_gen = falses(nbus)
   for r in axes(gen, 1)
      gen[r, 8] > 0.0 || continue
      k = get(index_of, Int(gen[r, 1]), 0)
      k == 0 && continue
      Pspec[k] += gen[r, 2] / baseMVA
      Qspec[k] += gen[r, 3] / baseMVA
      if !has_gen[k]
         Vm[k] = gen[r, 6]
         Qmin[k] = 0.0
         Qmax[k] = 0.0
         has_gen[k] = true
      end
      Qmin[k] += gen[r, 5] / baseMVA
      Qmax[k] += gen[r, 4] / baseMVA
   end
   for k = 1:nbus
      if bustype[k] == :pv && !has_gen[k]
         bustype[k] = :pq
      end
      if bustype[k] != :pv
         Qmin[k] = -1e9
         Qmax[k] = 1e9
      end
   end
   bustype[slack] = :slack

   function make_branches(unit::Symbol, sgn::Int, conv::Symbol)
      out = NamedTuple{(:i, :j, :r, :x, :b, :ratio, :shift_deg),Tuple{Int,Int,Float64,Float64,Float64,Float64,Float64}}[]
      for r in axes(br, 1)
         br[r, 11] > 0.0 || continue
         i = get(index_of, Int(br[r, 1]), 0)
         j = get(index_of, Int(br[r, 2]), 0)
         (i == 0 || j == 0) && continue
         ratio = br[r, 9] == 0.0 ? 1.0 : br[r, 9]
         conv == :inverse && (ratio = 1.0 / ratio)
         ang = sgn * br[r, 10]
         unit == :rad && (ang = rad2deg(ang))
         push!(out, transformer_branch(i, j; r = br[r, 3], x = br[r, 4], b = br[r, 5], ratio = ratio, shift_deg = ang))
      end
      return out
   end

   has_angle = any(r -> br[r, 10] != 0.0 && br[r, 11] > 0.0, 1:size(br, 1))
   has_ratio = any(r -> br[r, 9] != 0.0 && br[r, 9] != 1.0 && br[r, 11] > 0.0, 1:size(br, 1))
   units = angle_unit == :auto ? (has_angle ? (:deg, :rad) : (:deg,)) : (angle_unit,)
   signs = angle_sign == 0 ? (has_angle ? (1, -1) : (1,)) : (angle_sign,)
   convs = ratio_convention == :auto ? (has_ratio ? (:matpower, :inverse) : (:matpower,)) : (ratio_convention,)

   # (max, L1) mismatch of the stored solution: P at all non-slack buses, Q at PQ buses
   slack_bus = slack   # single assignment, so the closure does not box it
   function ref_mismatch(Y)
      Sc = calc_injections(Y, V_ref)
      mx = 0.0
      l1 = 0.0
      for k = 1:nbus
         k == slack_bus && continue
         dP = abs(real(Sc[k]) - Pspec[k])
         mx = max(mx, dP)
         l1 += dP
         if bustype[k] == :pq
            dQ = abs(imag(Sc[k]) - Qspec[k])
            mx = max(mx, dQ)
            l1 += dQ
         end
      end
      return mx, l1
   end

   trials = NamedTuple{(:angle_unit, :angle_sign, :ratio_convention, :ref_mismatch_pu, :ref_mismatch_l1_pu),Tuple{Symbol,Int,Symbol,Float64,Float64}}[]
   best = nothing
   for unit in units, sgn in signs, conv in convs
      branches = make_branches(unit, sgn, conv)
      Y = build_ybus(nbus, branches; bus_shunts = bus_shunts, sparse_output = true)
      mx, l1 = ref_mismatch(Y)
      push!(trials, (angle_unit = unit, angle_sign = sgn, ratio_convention = conv, ref_mismatch_pu = mx, ref_mismatch_l1_pu = l1))
      verbose >= 1 && @printf("matpower_case: angle=%-3s sign=%+d ratio=%-8s -> stored solution mismatch max = %.3e pu, L1 = %.3e pu\n", unit, sgn, conv, mx, l1)
      if best === nothing || l1 < best.l1
         best = (unit = unit, sgn = sgn, conv = conv, mis = mx, l1 = l1, branches = branches, Y = Y)
      end
   end
   Y = sparse_output ? best.Y : Matrix(best.Y)
   return (
      name = mp.name,
      Y = Y,
      bustype = bustype,
      Pspec = Pspec,
      Qspec = Qspec,
      Vm = Vm,
      Qmin = Qmin,
      Qmax = Qmax,
      slack = slack,
      labels = labels,
      baseMVA = baseMVA,
      branches = best.branches,
      bus_shunts = bus_shunts,
      V_ref = V_ref,
      conventions = (angle_unit = best.unit, angle_sign = best.sgn, ratio_convention = best.conv, ref_mismatch_pu = best.mis, ref_mismatch_l1_pu = best.l1, trials = trials),
   )
end
