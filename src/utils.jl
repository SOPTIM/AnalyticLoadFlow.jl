# =============================================================================
# File: src/utils.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Provides formatting, mismatch, stability, logging, and synthetic tiled-grid helper utilities.
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
    polar_str(V; digits_mag=6, digits_ang=4) -> String

Format a complex voltage phasor as `|V| ∠ ang°` (angle in degrees).
"""
function polar_str(V::ComplexF64; digits_mag::Int = 6, digits_ang::Int = 4)
   mag = abs(V)
   ang = angle(V) * 180 / pi
   return "|V|=$(round(mag, digits=digits_mag)) ∠ $(round(ang, digits=digits_ang))°"
end

"""
    print_bus_voltages(V; labels=nothing)

Print per-bus voltages in p.u. as `|V| ∠ ang°`.
`labels` may be a vector of bus names (length = length(V)).
"""
function print_bus_voltages(V::Vector{ComplexF64}; labels = nothing)
   n = length(V)
   for i = 1:n
      name = labels === nothing ? "Bus $i" : labels[i]
      println(rpad(name, 10), "  ", polar_str(V[i]))
   end
end

"""
    print_bus_voltages_kv(V, Vbase_LL_kV; labels=nothing)

Print per-bus voltages in kV (line-to-line base). Magnitude uses `abs(V)*Vbase_LL_kV`,
angle is in degrees.
"""
function print_bus_voltages_kv(V::Vector{ComplexF64}, Vbase_LL_kV::Float64; labels = nothing)
   n = length(V)
   for i = 1:n
      name = labels === nothing ? "Bus $i" : labels[i]
      mag_kV = abs(V[i]) * Vbase_LL_kV
      ang = angle(V[i]) * 180 / pi
      println(rpad(name, 10), "  |U|=$(round(mag_kV, digits=3)) kV ∠ $(round(ang, digits=4))°")
   end
end

"""
    print_total_line_losses(flows) -> tot

Compute and print total line losses. Requires `total_line_losses(flows)` returning
a struct/named-tuple with fields `Ploss_MW` and `Qloss_MVAr`.
"""
function print_total_line_losses(flows)
   tot = total_line_losses(flows)
   println("Total line losses: Ploss=$(round(tot.Ploss_MW, digits=3)) MW, Qloss=$(round(tot.Qloss_MVAr, digits=3)) MVAr")
   return tot
end

"""
    print_line_flows(flows; labels=nothing)

Print per-line directional flows and losses.
Each element `f` in `flows` must provide: `i`, `j`, `Pij_MW`, `Qij_MVAr`, `Pji_MW`, `Qji_MVAr`, `Ploss_MW`, `Qloss_MVAr`.
"""
function print_line_flows(flows; labels = nothing)
   for f in flows
      name_i = labels === nothing ? "Bus $(f.i)" : labels[f.i]
      name_j = labels === nothing ? "Bus $(f.j)" : labels[f.j]

      println("Line $(name_i) -> $(name_j):  P=$(round(f.Pij_MW, digits=3)) MW,  Q=$(round(f.Qij_MVAr, digits=3)) MVAr")
      println("Line $(name_j) -> $(name_i):  P=$(round(f.Pji_MW, digits=3)) MW,  Q=$(round(f.Qji_MVAr, digits=3)) MVAr")
      println(
         "Losses on $(name_i)-$(name_j):  Ploss=$(round(f.Ploss_MW, digits=3)) MW,  Qloss=$(round(f.Qloss_MVAr, digits=3)) MVAr",
      )
      println()
   end
end


# -----------------------------------------------------------------------------
# Logging control (to avoid I/O dominating benchmarks)
# -----------------------------------------------------------------------------
"""
    with_silent(f)

Execute function `f` with all logging suppressed to avoid I/O overhead in benchmarks.
"""
with_silent(f) = Logging.with_logger(NullLogger()) do
   f()
end

# -------------------------------------------------------------------------
# APSLF "stability" indicator from Padé poles (requires Vcoeff from solver)
# -------------------------------------------------------------------------
"""
    safe_get(x, sym::Symbol, default = nothing)

Safely get property `sym` from object `x`, returning `default` if property doesn't exist.
"""
safe_get(x, sym::Symbol, default = nothing) = (x !== nothing && sym in propertynames(x)) ? getproperty(x, sym) : default


"""
    st_level(dmin::Float64)

Convert stability distance to color-coded level (GRN/YEL/RED/NA).
"""
st_level(dmin::Float64) = isnan(dmin) ? "NA" : dmin < 0.10 ? "RED" : dmin < 0.30 ? "YEL" : "GRN"


"""
    diff_voltages(Va::Vector{ComplexF64}, Vb::Vector{ComplexF64})

Calculate maximum magnitude and angle differences between two voltage vectors.
"""
function diff_voltages(Va::Vector{ComplexF64}, Vb::Vector{ComplexF64})
   @assert length(Va) == length(Vb)
   mag_a = abs.(Va)
   mag_b = abs.(Vb)
   dV = maximum(abs.(mag_a .- mag_b))

   ang_a = angle.(Va)
   ang_b = angle.(Vb)  # radians
   d = ang_a .- ang_b
   d = mod.(d .+ pi, 2pi) .- pi            # wrap to [-pi, pi]
   dAng = maximum(abs.(rad2deg.(d)))
   return dV, dAng
end


"""
    max_mismatch_on_specY(spec, V::Vector{ComplexF64})

Calculate maximum P and Q mismatches for given voltages on specified admittance matrix.
"""
function max_mismatch_on_specY(spec, V::Vector{ComplexF64})
   S = APSLF.calc_injections(spec.Y, V)
   Pmis = real.(S) .- spec.Pspec
   Qmis = imag.(S) .- spec.Qspec

   Pmis[spec.slack] = 0.0
   Qmis[spec.slack] = 0.0
   @inbounds for i in eachindex(spec.bustype)
      spec.bustype[i] == :pv && (Qmis[i] = 0.0)
   end
   return maximum(abs.(Pmis)), maximum(abs.(Qmis))
end

"""
    max_mismatch(Y, bustype, Pspec, Qspec, V; slack)

Compute max active/reactive mismatch (p.u.) on PQ buses (and P mismatch on PV).
"""
function max_mismatch(
   Y::AbstractMatrix{ComplexF64},
   bustype::Vector{Symbol},
   Pspec::Vector{Float64},
   Qref::Vector{Float64},          # <-- reference Q (final schedule)
   V::Vector{ComplexF64};
   slack::Int,
)
   Sinj = APSLF.calc_injections(Matrix{ComplexF64}(Y), V)
   Pcalc = real.(Sinj)
   Qcalc = imag.(Sinj)

   maxP = 0.0
   maxQ = 0.0

   @inbounds for i in eachindex(bustype)
      i == slack && continue
      if bustype[i] == :pq
         maxP = max(maxP, abs(Pspec[i] - Pcalc[i]))
         maxQ = max(maxQ, abs(Qref[i] - Qcalc[i]))   # <-- compare to Qref
      elseif bustype[i] == :pv
         maxP = max(maxP, abs(Pspec[i] - Pcalc[i]))
      end
   end
   return (maxP, maxQ)
end




"""
    any_pv_at_qlimit(Y, bustype, Qmin, Qmax, V; eps = 1e-7)::Bool

Check if any PV bus is operating at its reactive power limit.
"""
function any_pv_at_qlimit(
   Y::AbstractMatrix{ComplexF64},
   bustype::Vector{Symbol},
   Qmin::Vector{Float64},
   Qmax::Vector{Float64},
   V::Vector{ComplexF64};
   eps::Float64 = 1e-7,
)::Bool
   Sinj = V .* conj.(Y * V)
   @inbounds for i in eachindex(bustype)
      bustype[i] == :pv || continue
      q = imag(Sinj[i])
      if (q <= Qmin[i] + eps) || (q >= Qmax[i] - eps)
         return true
      end
   end
   return false
end

# -----------------------------------------------------------------------------
# Synthetic tiled-grid demo network helpers
# -----------------------------------------------------------------------------

"""
    choose_tiled_grid_dimensions(max_buses::Int; aspect_ratio=1.0)

Choose `(rows, cols, nbus)` for the largest rectangular tiled grid with
`rows * cols <= max_buses`. The selected grid has at least two rows and two
columns; if multiple grids have the same bus count, the one closest to
`aspect_ratio = cols / rows` is chosen.

This helper makes scaling runs reproducible when a requested bus limit is not an
exact rectangle: the next lower feasible tiled size is used.
"""
function choose_tiled_grid_dimensions(max_buses::Int; aspect_ratio::Real = 1.0)
   max_buses >= 4 || error("max_buses must be >= 4")
   target_aspect = Float64(aspect_ratio)
   target_aspect > 0 || error("aspect_ratio must be positive")

   best = (rows = 2, cols = 2, nbus = 4, aspect_error = Inf)
   for rows = 2:max_buses
      max_cols = max_buses ÷ rows
      max_cols < 2 && break
      for cols = 2:max_cols
         nbus = rows * cols
         aspect = cols / rows
         aspect_error = abs(log(aspect / target_aspect))
         if nbus > best.nbus || (nbus == best.nbus && aspect_error < best.aspect_error)
            best = (rows = rows, cols = cols, nbus = nbus, aspect_error = aspect_error)
         end
      end
   end
   return (rows = best.rows, cols = best.cols, nbus = best.nbus)
end

@inline tiled_grid_bus_index(row::Int, col::Int, cols::Int)::Int = (row - 1) * cols + col

function _add_tiled_grid_line_stamp!(
   I::Vector{Int},
   J::Vector{Int},
   V::Vector{ComplexF64},
   from::Int,
   to::Int,
   z::ComplexF64,
   ysh::ComplexF64,
)
   y = inv(z)
   push!(I, from);
   push!(J, from);
   push!(V, y + ysh / 2)
   push!(I, to);
   push!(J, to);
   push!(V, y + ysh / 2)
   push!(I, from);
   push!(J, to);
   push!(V, -y)
   push!(I, to);
   push!(J, from);
   push!(V, -y)
   return nothing
end

"""
    build_tiled_grid_ybus(rows, cols; r, x, g=0.0, b=0.0)

Build the sparse Y-bus for a synthetic one-voltage-level demo grid.

Topology:
- every neighboring pair in horizontal direction is connected,
- every neighboring pair in vertical direction is connected,
- every rectangle has one diagonal from its upper-left to lower-right corner.

All connections are line PI equivalents with series impedance `r + im*x` and
total shunt admittance `g + im*b` split equally to both terminal buses. The
function returns `(Y, branch_count)`.
"""
function build_tiled_grid_ybus(rows::Int, cols::Int; r::Real, x::Real, g::Real = 0.0, b::Real = 0.0)
   rows >= 2 && cols >= 2 || error("rows and cols must be >= 2")
   r == 0.0 && x == 0.0 && error("r and x must not both be zero")

   nbus = rows * cols
   branch_count = rows * (cols - 1) + (rows - 1) * cols + (rows - 1) * (cols - 1)
   I = Vector{Int}(undef, 0)
   J = Vector{Int}(undef, 0)
   V = Vector{ComplexF64}(undef, 0)
   sizehint!(I, 4 * branch_count)
   sizehint!(J, 4 * branch_count)
   sizehint!(V, 4 * branch_count)

   z = ComplexF64(Float64(r), Float64(x))
   ysh = ComplexF64(Float64(g), Float64(b))

   # Horizontal tile edges.
   for row = 1:rows
      for col = 1:(cols-1)
         _add_tiled_grid_line_stamp!(
            I,
            J,
            V,
            tiled_grid_bus_index(row, col, cols),
            tiled_grid_bus_index(row, col + 1, cols),
            z,
            ysh,
         )
      end
   end
   # Vertical tile edges.
   for row = 1:(rows-1)
      for col = 1:cols
         _add_tiled_grid_line_stamp!(
            I,
            J,
            V,
            tiled_grid_bus_index(row, col, cols),
            tiled_grid_bus_index(row + 1, col, cols),
            z,
            ysh,
         )
      end
   end
   # One upper-left -> lower-right diagonal per rectangle.
   for row = 1:(rows-1)
      for col = 1:(cols-1)
         _add_tiled_grid_line_stamp!(
            I,
            J,
            V,
            tiled_grid_bus_index(row, col, cols),
            tiled_grid_bus_index(row + 1, col + 1, cols),
            z,
            ysh,
         )
      end
   end

   return sparse(I, J, V, nbus, nbus), branch_count
end

"""
    build_tiled_grid_spec(max_buses, cfg) -> spec, metadata

Create an APSLF-ready synthetic tiled-grid `spec` and metadata.

`cfg` is any object with the fields used by the example configuration:
`aspect_ratio`, `r`, `x`, `g`, `b`, `base_mva`,
`load_mw_per_right_corner`, `load_mvar_per_right_corner`,
`generation_balance`, `vm_slack`, and `vm_flat`.

The left side represents the generation side: the upper-left bus is slack and
the lower-left bus receives a scheduled PQ injection. The right side represents
the load side: the upper-right and lower-right buses receive equal PQ loads.
`generation_balance` is intentionally slightly below `1.0` by default so the
network is nearly balanced while the slack still covers a small residual and
line losses.
"""
function build_tiled_grid_spec(max_buses::Int, cfg)
   dims = choose_tiled_grid_dimensions(max_buses; aspect_ratio = cfg.aspect_ratio)
   Y, branch_count = build_tiled_grid_ybus(dims.rows, dims.cols; r = cfg.r, x = cfg.x, g = cfg.g, b = cfg.b)
   nbus = dims.nbus
   base_mva = cfg.base_mva

   top_left = tiled_grid_bus_index(1, 1, dims.cols)
   bottom_left = tiled_grid_bus_index(dims.rows, 1, dims.cols)
   top_right = tiled_grid_bus_index(1, dims.cols, dims.cols)
   bottom_right = tiled_grid_bus_index(dims.rows, dims.cols, dims.cols)

   bustype = fill(:pq, nbus)
   slack = top_left
   bustype[slack] = :slack

   Pspec = zeros(Float64, nbus)
   Qspec = zeros(Float64, nbus)
   load_p = cfg.load_mw_per_right_corner / base_mva
   load_q = cfg.load_mvar_per_right_corner / base_mva
   Pspec[top_right] = -load_p
   Pspec[bottom_right] = -load_p
   Qspec[top_right] = -load_q
   Qspec[bottom_right] = -load_q

   scheduled_gen_p = 2 * load_p * cfg.generation_balance
   scheduled_gen_q = 2 * load_q * cfg.generation_balance
   if bottom_left != slack
      Pspec[bottom_left] = scheduled_gen_p
      Qspec[bottom_left] = scheduled_gen_q
   end

   Vm = fill(Float64(cfg.vm_flat), nbus)
   Vm[slack] = Float64(cfg.vm_slack)
   Qmin = fill(-Inf, nbus)
   Qmax = fill(+Inf, nbus)
   labels = ["r$(row)c$(col)" for row = 1:dims.rows for col = 1:dims.cols]
   reference_V = ComplexF64.(Vm, zeros(Float64, nbus))

   spec = (
      Y = Y,
      bustype = bustype,
      slack = slack,
      Pspec = Pspec,
      Qspec = Qspec,
      Vm = Vm,
      Qmin = Qmin,
      Qmax = Qmax,
      labels = labels,
      baseMVA = base_mva,
      reference_V = reference_V,
      nbus = nbus,
   )
   metadata = (
      rows = dims.rows,
      cols = dims.cols,
      requested_max_buses = max_buses,
      actual_buses = nbus,
      branch_count = branch_count,
      generation_buses = (top_left, bottom_left),
      load_buses = (top_right, bottom_right),
      scheduled_generation_pu = scheduled_gen_p,
      scheduled_generation_q_pu = scheduled_gen_q,
      scheduled_load_pu = 2 * load_p,
      scheduled_load_q_pu = 2 * load_q,
   )
   return spec, metadata
end
