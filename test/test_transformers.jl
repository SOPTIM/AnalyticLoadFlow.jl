# =============================================================================
# File: test/test_transformers.jl
# Date: 2026-09-17
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Tests for transformer/PST branches, the deviation embedding in all
#          kernels, the phase ramp embedding, the regulated phase shifter loop
#          and the MATPOWER importer.
#          The germ and inner-solver tests are in test_germ_fix.jl.
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

using Test
using LinearAlgebra
using SparseArrays
using AnalyticLoadFlow

const A = AnalyticLoadFlow

_series_3bus() = build_ybus(3, [pi_branch(1, 2; r = 0.02, x = 0.10), pi_branch(2, 3; r = 0.03, x = 0.15), pi_branch(1, 3; r = 0.01, x = 0.08)])

@testset "branch builders and Y-bus stamping" begin
   br = demo_9bus_branches()
   @test length(br) == 9
   Y = build_ybus(9, br)
   @test maximum(abs.(Y .- A.demo_case_9bus().Y)) < 1e-14
   Ys = build_ybus(9, br; sparse_output = true)
   @test Ys isa SparseMatrixCSC{ComplexF64,Int}
   @test maximum(abs.(Matrix(Ys) .- Y)) < 1e-14
   Y0 = build_ybus(9, br; nominal = true)
   @test maximum(abs.(sum(Y0, dims = 2))) < 1e-12        # zero row sums
   # transformer with ratio and shift: non-symmetric, theory Section 6.5 entries
   t = 0.98 * cis(deg2rad(7.0))
   y = inv(complex(0.01, 0.1))
   tb = transformer_branch(1, 2; r = 0.01, x = 0.1, ratio = 0.98, shift_deg = 7.0)
   yii, yij, yji, yjj = A.branch_admittances(tb)
   @test yii ≈ y / abs2(t)
   @test yij ≈ -y / conj(t)
   @test yji ≈ -y / t
   @test yjj ≈ y
   Yt = build_ybus(2, [tb])
   @test Yt[1, 2] != Yt[2, 1]
   @test A.is_phase_shifter(tb) && A.is_transformer(tb) && !A.is_transformer(pi_branch(1, 2; r = 0.1, x = 0.2))
   sh = build_ybus(2, [tb]; bus_shunts = [(bus = 2, g = 0.0, b = 0.05)])
   @test sh[2, 2] ≈ Yt[2, 2] + 0.05im
   @test_throws ArgumentError build_ybus(2, [pi_branch(1, 3; r = 0.1, x = 0.2)])
end

@testset "PST case: exact solution, flows and angle sweep" begin
   c0 = demo_case_9bus_pst()
   @test maximum(abs.(c0.Y .- A.demo_case_9bus().Y)) < 1e-14
   P45 = Float64[]
   for φ in (-10.0, 0.0, 10.0)
      c = demo_case_9bus_pst(shift_deg = φ, enforce_q_limits = false, sparse_output = true)
      res = solve_pf_apslf(c; mode = :direct, order = 40, nr_polish = false)
      maxP, maxQ = A.compute_demo_mismatch(c, res)
      @test res.converged
      @test max(maxP, maxQ) < 1e-8
      flows = branch_flows(res.V, c.branches; baseMVA = c.baseMVA)
      @test length(flows) == 9
      # power balance: generation - load = losses
      Sinj = A.calc_injections(c.Y, res.V)
      @test abs(real(sum(Sinj)) - sum(f.Ploss_MW for f in flows) / c.baseMVA) < 1e-8
      push!(P45, branch_active_power(res.V, c.branches, 4, 5))
      @test branch_active_power(res.V, c.branches, 5, 4) ≈ real(flows[2].Sji_pu)
   end
   @test P45[1] > P45[2] > P45[3]   # negative shift on the 4-side pushes more power 4→5
   @test_throws ArgumentError branch_active_power(ones(ComplexF64, 9), c0.branches, 1, 9)
end

@testset "regulated PST outer loop" begin
   reg = solve_pf_pst_regulated(φ -> demo_case_9bus_pst(shift_deg = φ, enforce_q_limits = false), 4, 5, 0.6; order = 40, nr_polish = false)
   @test reg.converged
   @test abs(reg.P_pu - 0.6) < 1e-5
   @test reg.iterations <= 8
   @test reg.case.pst.shift_deg == reg.shift_deg
   # angle limit pins the result
   pinned = solve_pf_pst_regulated(φ -> demo_case_9bus_pst(shift_deg = φ, enforce_q_limits = false), 4, 5, 0.6; shift_min = -3.0, shift_max = 3.0, order = 40, nr_polish = false)
   @test !pinned.converged
   @test pinned.shift_deg == -3.0
end

@testset "MATPOWER importer with convention detection" begin
   # 4-bus case written in MATPOWER syntax: bus 2 fed through a PST (angle in radians, sign flipped),
   # bus 4 through an OLTC with the inverse ratio. The stored solution is produced by APSLF itself.
   base_branches = [
      pi_branch(1, 2; r = 0.01, x = 0.10, b = 0.02),
      transformer_branch(1, 3; r = 0.005, x = 0.05, shift_deg = 5.0),
      pi_branch(2, 3; r = 0.02, x = 0.15, b = 0.03),
      transformer_branch(3, 4; r = 0.0, x = 0.08, ratio = 0.95),
   ]
   shunts = [(bus = 2, g = 0.0, b = 0.05)]
   Y = build_ybus(4, base_branches; bus_shunts = shunts)
   spec = (Y = Y, bustype = [:slack, :pq, :pv, :pq], Pspec = [0.0, -0.5, 0.4, -0.3], Qspec = [0.0, -0.2, 0.0, -0.1],
      Vm = [1.02, 1.0, 1.01, 1.0], Qmin = [-1e9, -1e9, -1.0, -1e9], Qmax = [1e9, 1e9, 1.0, 1e9], slack = 1)
   res = solve_pf_apslf(spec; order = 40, nr_polish = true)
   @test res.converged
   Sinj = A.calc_injections(Y, res.V)
   path = joinpath(mktempdir(), "case4test.m")
   open(path, "w") do io
      println(io, "function mpc = case4test")
      println(io, "%   comment line ; with ] brackets")
      println(io, "mpc.version = '2';")
      println(io, "mpc.baseMVA = 100;")
      println(io, "%\tbus_i\ttype\tPd\tQd\tGs\tBs\tarea\tVm\tVa\tbaseKV\tzone\tVmax\tVmin")
      println(io, "mpc.bus = [")
      types = [3, 1, 2, 1, 4]
      for k = 1:4
         Pd = k == 3 ? 0.0 : -spec.Pspec[k] * 100
         Qd = k == 3 ? 0.0 : -spec.Qspec[k] * 100
         Bs = k == 2 ? 5.0 : 0.0
         println(io, "\t$(10k)\t$(types[k])\t$(Pd)\t$(Qd)\t0\t$(Bs)\t1\t$(abs(res.V[k]))\t$(rad2deg(angle(res.V[k])))\t110\t1\t1.1\t0.9;")
      end
      println(io, "\t99\t4\t0\t0\t0\t0\t1\t1\t0\t110\t1\t1.1\t0.9;   % isolated bus, dropped")
      println(io, "];")
      println(io, "mpc.gen = [")
      println(io, "\t10\t0\t0\t999\t-999\t1.02\t100\t1\t999\t-999\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;")
      println(io, "\t30\t40\t0\t100\t-100\t1.01\t100\t1\t999\t-999\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;")
      println(io, "\t30\t0\t0\t0\t0\t1.01\t100\t0\t999\t-999\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;   % out of service")
      println(io, "];")
      println(io, "mpc.branch = [")
      println(io, "\t10\t20\t0.01\t0.10\t0.02\t0\t0\t0\t0\t0\t1\t-360\t360;")
      println(io, "\t10\t30\t0.005\t0.05\t0\t0\t0\t0\t0\t$(-deg2rad(5.0))\t1\t-360\t360;")
      println(io, "\t20\t30\t0.02\t0.15\t0.03\t0\t0\t0\t0\t0\t1\t-360\t360;")
      println(io, "\t30\t40\t0\t0.08\t0\t0\t0\t0\t$(1/0.95)\t0\t1\t-360\t360;")
      println(io, "\t30\t40\t0\t0.08\t0\t0\t0\t0\t1\t0\t0\t-360\t360;   % out of service")
      println(io, "\t30\t99\t0\t0.08\t0\t0\t0\t0\t1\t0\t1\t-360\t360;   % to isolated bus, dropped")
      println(io, "];")
   end
   mp = parse_matpower_m(path)
   @test mp.baseMVA == 100.0
   @test size(mp.bus) == (5, 13) && size(mp.gen) == (3, 21) && size(mp.branch) == (6, 13)
   # experimental detection: the file is written in a non-MATPOWER convention
   c = matpower_case(path; angle_unit = :auto, angle_sign = 0, ratio_convention = :auto)
   @test size(c.Y) == (4, 4) && issparse(c.Y)
   @test c.bustype == [:slack, :pq, :pv, :pq]
   @test c.slack == 1
   @test c.labels == ["10", "20", "30", "40"]
   @test length(c.branches) == 4
   @test c.conventions.angle_unit == :rad
   @test c.conventions.angle_sign == -1
   @test c.conventions.ratio_convention == :inverse
   @test c.conventions.ref_mismatch_pu < 1e-6
   @test length(c.conventions.trials) == 8
   @test c.Pspec ≈ spec.Pspec && c.Qspec ≈ spec.Qspec
   @test c.Qmin[3] ≈ -1.0 && c.Qmax[3] ≈ 1.0
   @test c.Vm[3] ≈ 1.01 && c.Vm[1] ≈ 1.02
   @test c.bus_shunts == [(bus = 2, g = 0.0, b = 0.05)]
   @test maximum(abs.(Matrix(c.Y) .- Y)) < 1e-12
   res2 = solve_pf_apslf(c; order = 40, nr_polish = false)
   @test res2.converged
   @test maximum(abs.(res2.V .- res.V)) < 1e-6
   # default: MATPOWER convention, no detection (the stored state of this file does not fit it)
   c_deg = matpower_case(path)
   @test (c_deg.conventions.angle_unit, c_deg.conventions.angle_sign, c_deg.conventions.ratio_convention) == (:deg, 1, :matpower)
   @test length(c_deg.conventions.trials) == 1
   @test c_deg.conventions.ref_mismatch_pu > 1e-3
   @test_throws ArgumentError matpower_case(path; angle_unit = :grad)
   @test_throws ArgumentError parse_matpower_m(joinpath(dirname(path), "missing.m"))
end

@testset "deviation embedding (default germ) is exact in all kernels" begin
   case = A.demo_case_9bus()
   for Ymat in (case.Y, sparse(case.Y)), inner in (:pq, :direct_pv, :direct_pv_sparse)
      res = A.solve_pf_apslf_with_pv_q_limits(
         Ymat, case.bustype, case.Pspec, case.Qspec, case.Vm, case.Qmin, case.Qmax;
         slack = 1, Vslack = ComplexF64(case.Vm[1], 0), inner = inner, order = 40, nr_polish = false, max_outer = 30, use_sparse = issparse(Ymat),
      )
      @test res.converged
      @test max(A.compute_demo_mismatch(case, res)...) < 1e-8
      @test res.apslf_germ == :deviation
   end
   Y = _series_3bus()
   @test A.apslf_germ(Y, [2, 3], 1, 1.05 + 0im, :deviation) == fill(1.05 + 0im, 2)
   @test A.apslf_row_sums(Y) ≈ vec(sum(Y, dims = 2))
   ws = A.build_apslf_pq_workspace(Y; slack = 1, order = 20, germ = :deviation)
   @test_throws ArgumentError A.apslf_pq_solve!(ws, ComplexF64[0, -0.1, -0.1]; germ = :noload)
   c = demo_case_9bus_pst(shift_deg = 8.0, enforce_q_limits = false, sparse_output = true)
   res = solve_pf_apslf(c; order = 40, nr_polish = false)
   @test res.converged && max(A.compute_demo_mismatch(c, res)...) < 1e-10
   res_nl = solve_pf_apslf(c; order = 40, nr_polish = false, germ = :noload)
   @test maximum(abs.(res.V .- res_nl.V)) < 1e-8   # both variants describe the same function at s = 1
end

# Newton reference (rectangular NR of the package, flat start, no Q limits)
function _newton_ref(Y, bt, P, Q, Vm, slack)
   Vsl = ComplexF64(Vm[slack], 0.0)
   V = [i == slack ? Vsl : complex(lowercase(string(bt[i])) == "pv" ? Vm[i] : 1.0) for i in eachindex(bt)]
   cache = A.build_nr_rect_cache(copy(bt), slack)
   for _ = 1:30
      dx, _ = A.nr_refine_step_rect!(cache, A.sparse_Y(Y), V, P, Q, Vm; Vslack = Vsl)
      dx < 1e-14 && break
   end
   return V
end

# radial 1-2-3-4, bus 1 slack, shifter with series impedance z on 2-3 (tap on side 2)
_radial_pst(shift_deg, z) = build_ybus(
   4,
   [pi_branch(1, 2; r = 0.01, x = 0.1, b = 0.02), transformer_branch(2, 3; r = real(z), x = imag(z), shift_deg = shift_deg), pi_branch(3, 4; r = 0.01, x = 0.1, b = 0.02)];
   sparse_output = true,
)

@testset "phase ramp embedding (shift_embedding)" begin
   bt = [:slack, :pq, :pq, :pq]
   P = [0.0, -0.3, -0.5, -0.4]
   Q = [0.0, -0.1, -0.2, -0.1]
   Vm = ones(4)
   Qlim = fill(10.0, 4)
   solve(Y; kw...) = A.solve_pf_apslf_with_pv_q_limits(Y, bt, P, Q, Vm, -Qlim, Qlim; order = 24, nr_polish = false, enforce_q_limits = false, kw...)

   # small shifter in the 9-bus case: the default (:auto) ramps and matches Newton
   c = demo_case_9bus_pst(shift_deg = 30.0, enforce_q_limits = false, sparse_output = true)
   res = solve_pf_apslf(c; order = 24, nr_polish = false)
   @test res.shift_embedding == :phase_ramp && isempty(res.shift_skipped)
   @test res.converged
   Vn = _newton_ref(c.Y, c.bustype, c.Pspec, c.Qspec, c.Vm, c.slack)
   @test maximum(abs.(res.V .- Vn)) < 1e-10

   # shifts above 90 degrees (vector groups): the radial net only turns buses 3 and 4 by -φ,
   # for inductive and capacitive series elements with losses
   for z in (0.005 + 0.08im, 0.005 - 0.08im)
      V0 = solve(_radial_pst(0.0, z)).V
      for φ in (150.0, -150.0, 170.0)
         Y = _radial_pst(φ, z)
         Vref = copy(V0)
         Vref[3:4] .*= cis(-deg2rad(φ))
         res = solve(Y)
         @test res.shift_embedding == :phase_ramp
         @test res.converged
         @test maximum(abs.(res.V .- Vref)) < 1e-10
         sp = A.apslf_shift_split(Y, 24)
         @test real(sp.Ybase[2, 3]) < 0 && sp.Ybase[2, 3] == sp.Ybase[3, 2]   # physical sign kept
         ΔY1 = sparse(sp.rows, sp.cols, vec(sum(sp.coeff, dims = 2)), 4, 4)
         @test maximum(abs.(sp.Ybase .+ ΔY1 .- Y)) < 1e-12 * maximum(abs.(Y))
      end
   end

   # a lossless capacitive pair has no clear sign: not ramped, the result is the :off kernel
   Y = _radial_pst(150.0, 0.0 - 0.08im)
   res = solve(Y)
   @test res.shift_embedding == :off
   @test [s.reason for s in res.shift_skipped] == [:lossless_sign]
   @test res.V == solve(Y; shift_embedding = :off).V

   # no shifter: :auto is the unchanged kernel, bit for bit
   c9 = A.demo_case_9bus()
   r_auto = solve_pf_apslf(c9; order = 24, nr_polish = false, return_coeffs = true)
   r_off = solve_pf_apslf(c9; order = 24, nr_polish = false, return_coeffs = true, shift_embedding = :off)
   @test r_auto.shift_embedding == :off
   @test r_auto.V == r_off.V && r_auto.Vcoeff == r_off.Vcoeff

   # an entry without structural partner: :auto falls back to :off instead of throwing
   Yd = Matrix(_radial_pst(30.0, 0.005 + 0.08im) + build_ybus(4, [pi_branch(1, 4; r = 0.02, x = 0.2)]; sparse_output = true))
   Yd[3, 2] = 0
   Yd = dropzeros!(sparse(Yd))
   res = solve(Yd)
   @test res.shift_embedding == :off
   @test [(s.i, s.j, s.reason) for s in res.shift_skipped] == [(2, 3, :no_partner)]
   @test res.V == solve(Yd; shift_embedding = :off).V

   # the outer PQ loop has no ramp: :auto runs :off there, an explicit :phase_ramp is an error
   Y = _radial_pst(30.0, 0.005 + 0.08im)
   @test solve(Y; inner = :pq).shift_embedding == :off
   @test_throws ArgumentError solve(Y; inner = :pq, shift_embedding = :phase_ramp)
   @test_throws ArgumentError solve(Y; germ = :noload, shift_embedding = :phase_ramp)
   @test solve(Y; germ = :noload).shift_embedding == :off
   @test_throws ArgumentError solve(Y; shift_embedding = :sym_rowsum)
end

const _PEGASE_1354 = normpath(joinpath(@__DIR__, "..", "data", "_downloaded", "case1354pegase.m"))
if isfile(_PEGASE_1354)
   @testset "PEGASE 1354 (local download): MATPOWER conventions and pure APSLF solution" begin
      c = matpower_case(_PEGASE_1354)
      @test size(c.Y) == (1354, 1354)
      @test (c.conventions.angle_unit, c.conventions.angle_sign, c.conventions.ratio_convention) == (:deg, 1, :matpower)
      res = solve_pf_apslf(c; order = 40, nr_polish = false, enforce_q_limits = false)
      @test res.converged
      @test res.effective_mode == :direct
      @test max(A.compute_demo_mismatch(c, res)...) < 1e-8
   end
else
   @info "PEGASE test skipped (no data/_downloaded/case1354pegase.m); run examples/pegase_matpower_demo.jl to download it"
end
