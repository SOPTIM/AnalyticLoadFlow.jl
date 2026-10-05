# =============================================================================
# File: src/transformers.jl
# Date: 2026-09-17
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Branch model with transformer ratio and phase shift (PST), Y-bus
#          stamping (dense or sparse), branch flows, a 9-bus PST demo case and
#          an outer loop for a regulated phase shifter (theory Section 6.5).
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
    pi_branch(i, j; r, x, b = 0.0)

π-model line between buses `i` and `j` with series impedance `r + jx` (pu) and
total line charging susceptance `b` (pu, split equally on both ends).
Returns the branch NamedTuple used by [`build_ybus`](@ref) and [`branch_flows`](@ref).
"""
pi_branch(i::Integer, j::Integer; r::Real, x::Real, b::Real = 0.0) =
   (i = Int(i), j = Int(j), r = Float64(r), x = Float64(x), b = Float64(b), ratio = 1.0, shift_deg = 0.0)

"""
    transformer_branch(i, j; r, x, b = 0.0, ratio = 1.0, shift_deg = 0.0)

Transformer between buses `i` and `j` with complex tap `t = ratio · e^{jφ}` on
side `i` (theory Section 6.5). `ratio ≠ 1` is an off-nominal ratio (OLTC),
`shift_deg ≠ 0` a phase shift (PST). The branch contributes

    I_i = (y + jb/2)/|t|² · V_i  −  y/conj(t) · V_j
    I_j =      −y/t · V_i        +  (y + jb/2) · V_j

with `y = 1/(r + jx)`. For `ratio = 1, shift_deg = 0` this is [`pi_branch`](@ref).
"""
transformer_branch(i::Integer, j::Integer; r::Real, x::Real, b::Real = 0.0, ratio::Real = 1.0, shift_deg::Real = 0.0) =
   (i = Int(i), j = Int(j), r = Float64(r), x = Float64(x), b = Float64(b), ratio = Float64(ratio), shift_deg = Float64(shift_deg))

is_transformer(br) = br.ratio != 1.0 || br.shift_deg != 0.0
is_phase_shifter(br) = br.shift_deg != 0.0

"""
    branch_admittances(br; nominal = false) -> (yii, yij, yji, yjj)

Entries of the 2×2 branch admittance matrix of a branch from
[`pi_branch`](@ref) / [`transformer_branch`](@ref). With `nominal = true` the
tap is forced to `1∠0°` and the shunt to zero, which gives the zero-row-sum
series part `Y0` used by the deviation embedding (theory Section 6.5, variant 1).
"""
function branch_admittances(br; nominal::Bool = false)
   y = inv(complex(br.r, br.x))
   ysh = nominal ? 0.0im : 0.5im * br.b
   t = nominal ? 1.0 + 0.0im : br.ratio * cis(deg2rad(br.shift_deg))
   yii = (y + ysh) / abs2(t)
   yij = -y / conj(t)
   yji = -y / t
   yjj = y + ysh
   return yii, yij, yji, yjj
end

"""
    build_ybus(nbus, branches; bus_shunts = (), nominal = false, sparse_output = false)

Bus admittance matrix from a vector of branches ([`pi_branch`](@ref),
[`transformer_branch`](@ref)) and optional bus shunts given as NamedTuples
`(bus, g, b)` in pu. A phase shifter makes `Y` non-symmetric (`Y[i,j] ≠ Y[j,i]`).

- `nominal = true` returns the series-only matrix `Y0` with every transformer at
  `1∠0°` and all shunts removed; `Y0` has zero row sums (theory Section 6.5).
- `sparse_output = true` returns a `SparseMatrixCSC{ComplexF64,Int}`, the
  preferred form for large networks and for the sparse direct PV kernel.
"""
function build_ybus(nbus::Integer, branches::AbstractVector; bus_shunts = (), nominal::Bool = false, sparse_output::Bool = false)
   nb = length(branches)
   I = Vector{Int}(undef, 4 * nb)
   J = Vector{Int}(undef, 4 * nb)
   Vv = Vector{ComplexF64}(undef, 4 * nb)
   k = 0
   for br in branches
      1 <= br.i <= nbus && 1 <= br.j <= nbus || throw(ArgumentError("branch ($(br.i),$(br.j)) outside 1:$(nbus)"))
      br.i != br.j || throw(ArgumentError("branch ($(br.i),$(br.j)) connects a bus to itself"))
      yii, yij, yji, yjj = branch_admittances(br; nominal = nominal)
      I[k+1] = br.i; J[k+1] = br.i; Vv[k+1] = yii
      I[k+2] = br.i; J[k+2] = br.j; Vv[k+2] = yij
      I[k+3] = br.j; J[k+3] = br.i; Vv[k+3] = yji
      I[k+4] = br.j; J[k+4] = br.j; Vv[k+4] = yjj
      k += 4
   end
   if !nominal
      for sh in bus_shunts
         push!(I, sh.bus); push!(J, sh.bus); push!(Vv, complex(Float64(sh.g), Float64(sh.b)))
      end
   end
   Ys = sparse(I, J, Vv, Int(nbus), Int(nbus))
   return sparse_output ? Ys : Matrix(Ys)
end

"""
    branch_flows(V, branches; baseMVA = 100.0)

Complex power flows on every branch (line or transformer) for the bus voltage
vector `V`. `Sij_pu` is measured at side `i` (positive into the branch),
`Sji_pu` at side `j`, `Sloss_pu = Sij + Sji`. MW/MVAr values use `baseMVA`.
"""
function branch_flows(V::AbstractVector{<:Complex}, branches::AbstractVector; baseMVA::Real = 100.0)
   flows = Vector{NamedTuple}(undef, length(branches))
   for (k, br) in enumerate(branches)
      yii, yij, yji, yjj = branch_admittances(br)
      Vi = ComplexF64(V[br.i])
      Vj = ComplexF64(V[br.j])
      Iij = yii * Vi + yij * Vj
      Iji = yji * Vi + yjj * Vj
      Sij = Vi * conj(Iij)
      Sji = Vj * conj(Iji)
      Sloss = Sij + Sji
      flows[k] = (
         i = br.i,
         j = br.j,
         ratio = br.ratio,
         shift_deg = br.shift_deg,
         Sij_pu = Sij,
         Sji_pu = Sji,
         Sloss_pu = Sloss,
         Pij_MW = real(Sij) * baseMVA,
         Qij_MVAr = imag(Sij) * baseMVA,
         Pji_MW = real(Sji) * baseMVA,
         Qji_MVAr = imag(Sji) * baseMVA,
         Ploss_MW = real(Sloss) * baseMVA,
         Qloss_MVAr = imag(Sloss) * baseMVA,
      )
   end
   return flows
end

"""
    branch_active_power(V, branches, i, j) -> P_ij (pu)

Active power entering the branch `(i, j)` at side `i`. The branch is looked up
in `branches` in the given orientation; the reversed orientation is accepted
and then the flow at side `j` of that branch is returned.
"""
function branch_active_power(V::AbstractVector{<:Complex}, branches::AbstractVector, i::Integer, j::Integer)
   for br in branches
      if br.i == i && br.j == j
         yii, yij, _, _ = branch_admittances(br)
         Vi = ComplexF64(V[br.i]); Vj = ComplexF64(V[br.j])
         return real(Vi * conj(yii * Vi + yij * Vj))
      elseif br.i == j && br.j == i
         _, _, yji, yjj = branch_admittances(br)
         Vi = ComplexF64(V[br.i]); Vj = ComplexF64(V[br.j])
         return real(Vj * conj(yji * Vi + yjj * Vj))
      end
   end
   throw(ArgumentError("no branch ($(i),$(j)) in branch list"))
end

"""
    print_branch_flows(flows; labels = nothing, io = stdout)

Console table for the result of [`branch_flows`](@ref).
"""
function print_branch_flows(flows; labels = nothing, io::IO = stdout)
   name(b) = labels === nothing ? string(b) : string(labels[b])
   println(io, "from        to          type   ratio   shift°     P_ij MW    Q_ij MVAr    P_ji MW    Q_ji MVAr   loss MW")
   for f in flows
      typ = f.shift_deg != 0.0 ? "PST" : (f.ratio != 1.0 ? "TR" : "line")
      @printf(io, "%-10s  %-10s  %-4s  %6.4f  %7.3f  %10.4f  %10.4f  %10.4f  %10.4f  %9.5f\n",
         name(f.i), name(f.j), typ, f.ratio, f.shift_deg, f.Pij_MW, f.Qij_MVAr, f.Pji_MW, f.Qji_MVAr, f.Ploss_MW)
   end
   return nothing
end

# -----------------------------------------------------------------------------
# 9-bus teaching case with a phase-shifting transformer
# -----------------------------------------------------------------------------

"""
    demo_9bus_branches()

Branch list of the 9-bus teaching case ([`demo_case_9bus`](@ref)) as
[`pi_branch`](@ref) NamedTuples.
"""
function demo_9bus_branches()
   return [
      pi_branch(1, 4; r = 0.0000, x = 0.0576, b = 0.0000),
      pi_branch(4, 5; r = 0.0170, x = 0.0920, b = 0.1580),
      pi_branch(5, 6; r = 0.0390, x = 0.1700, b = 0.3580),
      pi_branch(3, 6; r = 0.0000, x = 0.0586, b = 0.0000),
      pi_branch(6, 7; r = 0.0119, x = 0.1008, b = 0.2090),
      pi_branch(7, 8; r = 0.0085, x = 0.0720, b = 0.1490),
      pi_branch(8, 2; r = 0.0000, x = 0.0625, b = 0.0000),
      pi_branch(8, 9; r = 0.0320, x = 0.1610, b = 0.3060),
      pi_branch(9, 4; r = 0.0100, x = 0.0850, b = 0.1760),
   ]
end

"""
    demo_case_9bus_pst(; shift_deg = 0.0, ratio = 1.0, pst = (4, 5), sparse_output = false, enforce_q_limits = true)

The 9-bus teaching case with a phase-shifting transformer lumped into the
branch `pst` (default: the line 4–5, tap on bus 4). The ring 4-5-6-7-8-9-4 has
two parallel paths from the generator buses to each load, so the PST angle
redistributes the active power between the paths without changing the total.

Returns the demo data contract of [`demo_case_9bus`](@ref) plus the fields
`branches`, `pst` and `Y0` (series-only zero-row-sum matrix, theory Section 6.5).
With `enforce_q_limits = false` the PV reactive limits are opened so that the
angle sweep is not disturbed by PV→PQ switching.
"""
function demo_case_9bus_pst(;
   shift_deg::Real = 0.0,
   ratio::Real = 1.0,
   pst::Tuple{Int,Int} = (4, 5),
   sparse_output::Bool = false,
   enforce_q_limits::Bool = true,
)
   base = demo_case_9bus()
   branches = demo_9bus_branches()
   idx = findfirst(br -> (br.i, br.j) == pst || (br.j, br.i) == pst, branches)
   idx === nothing && throw(ArgumentError("branch $(pst) not found in the 9-bus case"))
   ln = branches[idx]
   branches[idx] = transformer_branch(pst[1], pst[2]; r = ln.r, x = ln.x, b = ln.b, ratio = ratio, shift_deg = shift_deg)
   Y = build_ybus(9, branches; sparse_output = sparse_output)
   Y0 = build_ybus(9, branches; nominal = true, sparse_output = sparse_output)
   Qmin = enforce_q_limits ? base.Qmin : fill(-1e9, 9)
   Qmax = enforce_q_limits ? base.Qmax : fill(1e9, 9)
   return merge(base, (Y = Y, Y0 = Y0, Qmin = Qmin, Qmax = Qmax, branches = branches, pst = (i = pst[1], j = pst[2], ratio = Float64(ratio), shift_deg = Float64(shift_deg))))
end

# -----------------------------------------------------------------------------
# Regulated phase shifter: outer loop on the PST angle (theory Section 6.5)
# -----------------------------------------------------------------------------

"""
    solve_pf_pst_regulated(case_of_shift, i, j, P_target_pu;
                           shift0 = 0.0, shift_min = -30.0, shift_max = 30.0,
                           p_tol = 1e-5, max_iter = 25, first_step_deg = 1.0,
                           verbose = 0, solver_kwargs...)

Outer loop for a **regulated** phase shifter. The angle enters the Y-bus
through `e^{jφ}`, i.e. not polynomially in the embedding parameter, so it is
adjusted outside the series recursion (theory Section 6.5): the secant method
iterates the angle `φ` until the active power `P_ij(φ)` entering branch `(i, j)`
at side `i` meets `P_target_pu`, and APSLF is re-run for every trial angle.

`case_of_shift(φ_deg)` must return a case NamedTuple with a `branches` field,
e.g. `φ -> demo_case_9bus_pst(shift_deg = φ)`. Remaining keyword arguments are
forwarded to [`solve_pf_apslf`](@ref).

Returns a NamedTuple `(shift_deg, P_pu, converged, iterations, history, res, case)`.
"""
function solve_pf_pst_regulated(
   case_of_shift,
   i::Integer,
   j::Integer,
   P_target_pu::Real;
   shift0::Real = 0.0,
   shift_min::Real = -30.0,
   shift_max::Real = 30.0,
   p_tol::Real = 1e-5,
   max_iter::Int = 25,
   first_step_deg::Real = 1.0,
   verbose::Int = 0,
   solver_kwargs...,
)
   history = NamedTuple{(:iter, :shift_deg, :P_pu, :converged),Tuple{Int,Float64,Float64,Bool}}[]
   function evaluate(φ::Float64, it::Int)
      local case = case_of_shift(φ)
      local res = solve_pf_apslf(case; solver_kwargs...)
      P = branch_active_power(res.V, case.branches, i, j)
      push!(history, (iter = it, shift_deg = φ, P_pu = P, converged = res.converged))
      verbose >= 1 && @printf("[pst %2d] shift = %9.4f°   P_%d%d = %10.6f pu   target = %10.6f pu   pf converged = %s\n", it, φ, i, j, P, P_target_pu, res.converged)
      return case, res, P
   end
   φ0 = clamp(Float64(shift0), shift_min, shift_max)
   case, res, P0 = evaluate(φ0, 0)
   f0 = P0 - P_target_pu
   if abs(f0) <= p_tol
      return (shift_deg = φ0, P_pu = P0, converged = res.converged, iterations = 0, history = history, res = res, case = case)
   end
   φ1 = clamp(φ0 + Float64(first_step_deg), shift_min, shift_max)
   φ1 == φ0 && (φ1 = clamp(φ0 - Float64(first_step_deg), shift_min, shift_max))
   case, res, P1 = evaluate(φ1, 1)
   f1 = P1 - P_target_pu
   it = 1
   while abs(f1) > p_tol && it < max_iter
      it += 1
      if abs(f1 - f0) < 1e-14
         break
      end
      φ2 = φ1 - f1 * (φ1 - φ0) / (f1 - f0)
      φ2 = clamp(φ2, shift_min, shift_max)
      case, res, P2 = evaluate(φ2, it)
      φ0, f0 = φ1, f1
      φ1, f1 = φ2, P2 - P_target_pu
      if (φ1 == shift_min || φ1 == shift_max) && φ0 == φ1
         break   # pinned at the angle limit
      end
   end
   return (shift_deg = φ1, P_pu = f1 + P_target_pu, converged = res.converged && abs(f1) <= p_tol, iterations = it, history = history, res = res, case = case)
end
