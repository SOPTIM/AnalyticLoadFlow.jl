# =============================================================================
# File: src/double_double.jl
# Date: 2026-11-11
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose:
# Double-double arithmetic (about 32 significant digits) for the direct PV
# kernel with `precision = :double_double`: real and complex numbers as
# unevaluated sums of two Float64, a dense LU for the Padé systems and the
# recursion with iterative refinement on the Float64 LU of the kernel.
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

# Error-free transformations (Dekker, Knuth); _two_prod relies on fma.
@inline function _two_sum(a::Float64, b::Float64)
   s = a + b
   bb = s - a
   return s, (a - (s - bb)) + (b - bb)
end

@inline function _quick_two_sum(a::Float64, b::Float64)
   s = a + b
   return s, b - (s - a)
end

@inline function _two_prod(a::Float64, b::Float64)
   p = a * b
   return p, fma(a, b, -p)
end

"""
    DD(hi, lo)

Double-double number `hi + lo` with `|lo| <= ulp(hi)/2`, about 32 significant
digits. Supports `+`, `-`, `*`, `/` with `DD` and `Float64`; `Float64(x)`
rounds back.
"""
struct DD
   hi::Float64
   lo::Float64
end

DD(x::Float64) = DD(x, 0.0)
Base.Float64(x::DD) = x.hi + x.lo
Base.zero(::Type{DD}) = DD(0.0, 0.0)
Base.BigFloat(x::DD) = BigFloat(x.hi) + BigFloat(x.lo)

@inline function Base.:+(a::DD, b::DD)
   s1, s2 = _two_sum(a.hi, b.hi)
   t1, t2 = _two_sum(a.lo, b.lo)
   s2 += t1
   s1, s2 = _quick_two_sum(s1, s2)
   s2 += t2
   s1, s2 = _quick_two_sum(s1, s2)
   return DD(s1, s2)
end

@inline Base.:-(a::DD) = DD(-a.hi, -a.lo)
@inline Base.:-(a::DD, b::DD) = a + (-b)

@inline function Base.:*(a::DD, b::DD)
   p1, p2 = _two_prod(a.hi, b.hi)
   p2 += a.hi * b.lo + a.lo * b.hi
   p1, p2 = _quick_two_sum(p1, p2)
   return DD(p1, p2)
end

@inline function Base.:*(a::DD, b::Float64)
   p1, p2 = _two_prod(a.hi, b)
   p2 += a.lo * b
   p1, p2 = _quick_two_sum(p1, p2)
   return DD(p1, p2)
end
@inline Base.:*(a::Float64, b::DD) = b * a

@inline function Base.:/(a::DD, b::DD)
   q1 = a.hi / b.hi
   r = a - b * q1
   q2 = r.hi / b.hi
   r = r - b * q2
   q3 = r.hi / b.hi
   q1, q2 = _quick_two_sum(q1, q2)
   return DD(q1, q2) + DD(q3)
end

"""
    CDD(re, im)

Complex double-double number. Supports `+`, `-`, `*`, `/`, `conj`;
`ComplexF64(z)` rounds back.
"""
struct CDD
   re::DD
   im::DD
end

CDD(z::ComplexF64) = CDD(DD(real(z)), DD(imag(z)))
CDD(x::Float64) = CDD(DD(x), DD(0.0))
Base.ComplexF64(z::CDD) = complex(Float64(z.re), Float64(z.im))
Base.zero(::Type{CDD}) = CDD(zero(DD), zero(DD))
Base.conj(z::CDD) = CDD(z.re, -z.im)

@inline Base.:+(a::CDD, b::CDD) = CDD(a.re + b.re, a.im + b.im)
@inline Base.:-(a::CDD, b::CDD) = CDD(a.re - b.re, a.im - b.im)
@inline Base.:-(a::CDD) = CDD(-a.re, -a.im)
@inline Base.:*(a::CDD, b::CDD) = CDD(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
@inline Base.:*(a::CDD, b::DD) = CDD(a.re * b, a.im * b)
@inline Base.:*(a::DD, b::CDD) = b * a
@inline Base.:*(a::CDD, b::ComplexF64) = a * CDD(b)

@inline function Base.:/(a::CDD, b::CDD)
   den = b.re * b.re + b.im * b.im
   return CDD((a.re * b.re + a.im * b.im) / den, (a.im * b.re - a.re * b.im) / den)
end

@inline _abs2_hi(z::CDD) = z.re.hi * z.re.hi + z.im.hi * z.im.hi
@inline _isfinite_cdd(z::CDD) = isfinite(z.re.hi) && isfinite(z.re.lo) && isfinite(z.im.hi) && isfinite(z.im.lo)

# Solve A*x = b in place (LU with partial pivoting); the solution is left in
# b. Returns false on a zero or non-finite pivot.
function _lu_solve_dd!(A::Matrix{CDD}, b::Vector{CDD})
   n = length(b)
   @inbounds for k = 1:n
      piv = k
      pmax = _abs2_hi(A[k, k])
      for i = (k+1):n
         v = _abs2_hi(A[i, k])
         if v > pmax
            pmax = v
            piv = i
         end
      end
      (pmax > 0.0 && isfinite(pmax)) || return false
      if piv != k
         for j = k:n
            A[k, j], A[piv, j] = A[piv, j], A[k, j]
         end
         b[k], b[piv] = b[piv], b[k]
      end
      akk = A[k, k]
      for i = (k+1):n
         f = A[i, k] / akk
         for j = (k+1):n
            A[i, j] -= f * A[k, j]
         end
         b[i] -= f * b[k]
      end
   end
   @inbounds for k = n:-1:1
      acc = b[k]
      for j = (k+1):n
         acc -= A[k, j] * b[j]
      end
      b[k] = acc / A[k, k]
   end
   return true
end

# Padé [L/M] value at s = 1 of the series c (c[n+1] = c_n), same construction
# as pade_eval. Returns nothing where the system is singular or the value is
# not finite.
function _pade_value_dd(c::AbstractVector{CDD}, L::Int, M::Int)
   if M == 0
      num = zero(CDD)
      @inbounds for n = 0:L
         num += c[n+1]
      end
      return num
   end
   A = Matrix{CDD}(undef, M, M)
   b = Vector{CDD}(undef, M)
   @inbounds for row = 1:M
      n = L + row
      b[row] = -c[n+1]
      for j = 1:M
         A[row, j] = c[(n-j)+1]
      end
   end
   _lu_solve_dd!(A, b) || return nothing
   num = zero(CDD)
   @inbounds for n = 0:L
      acc = c[n+1]
      for j = 1:min(n, M)
         acc += b[j] * c[(n-j)+1]
      end
      num += acc
   end
   den = CDD(1.0)
   @inbounds for j = 1:M
      den += b[j]
   end
   v = num / den
   return _isfinite_cdd(v) ? v : nothing
end

# Series value at s = 1: Padé [N-N÷2 / N÷2] as in evaluate_series, or the
# plain sum for taylor. Falls back to the Float64 evaluation of the rounded
# coefficients where the double-double Padé fails.
function _series_value_dd(c::AbstractVector{CDD}, taylor::Bool, eval_opts::APSLFEvaluationOptions)
   N = length(c) - 1
   M = taylor ? 0 : div(N, 2)
   v = _pade_value_dd(c, N - M, M)
   v === nothing || return ComplexF64(v)
   return evaluate_series(ComplexF64[ComplexF64(x) for x in c], eval_opts; collect_logs = false).voltage
end

# r = rhs - A*x with x in double-double
function _residual_dd!(r::Vector{DD}, A::SparseMatrixCSC{Float64,Int}, x::Vector{DD}, rhs::Vector{DD})
   copyto!(r, rhs)
   rv = rowvals(A)
   nz = nonzeros(A)
   @inbounds for j in axes(A, 2)
      xj = -x[j]
      for p in nzrange(A, j)
         r[rv[p]] += xj * nz[p]
      end
   end
   return r
end

# Solve A*x = rhs to double-double accuracy with the Float64 factorization F
# of A: iterative refinement, residual in double-double. Stops when the next
# correction would fall below the double-double rounding level or the
# corrections stop shrinking.
function _refined_solve_dd!(
   x::Vector{DD},
   A::SparseMatrixCSC{Float64,Int},
   F,
   rhs::Vector{DD},
   r::Vector{DD},
   r64::Vector{Float64},
   dx::Vector{Float64};
   maxit::Int = 20,
)
   n = length(x)
   @inbounds for i = 1:n
      r64[i] = Float64(rhs[i])
   end
   ldiv!(dx, F, r64)
   xmax = 0.0
   @inbounds for i = 1:n
      x[i] = DD(dx[i])
      xmax = max(xmax, abs(dx[i]))
   end
   xmax > 0.0 || return x
   prev = xmax
   for _ = 1:maxit
      _residual_dd!(r, A, x, rhs)
      @inbounds for i = 1:n
         r64[i] = Float64(r[i])
      end
      ldiv!(dx, F, r64)
      dmax = 0.0
      @inbounds for i = 1:n
         x[i] += DD(dx[i])
         dmax = max(dmax, abs(dx[i]))
      end
      isfinite(dmax) || break
      # contraction per step is dmax/prev; the error left is about dmax times that
      rate = dmax / prev
      (rate > 0.5 || dmax * rate <= 1e-31 * xmax) && break
      prev = dmax
   end
   return x
end

# Recursion and evaluation of apslf_pf_pv_direct in double-double. The
# system matrix Asys and its Float64 LU Fsys are those of the Float64 kernel;
# every order is solved with iterative refinement, the convolutions and the
# Padé evaluation run in double-double. Vcoeff, Wcoeff and Qcoeff receive the
# coefficients rounded to Float64. Returns (V, Qpv).
function _apslf_direct_dd!(
   Vcoeff::Matrix{ComplexF64},
   Wcoeff::Matrix{ComplexF64},
   Qcoeff::Matrix{Float64},
   Asys::SparseMatrixCSC{Float64,Int},
   Fsys,
   bustype::Vector{Symbol},
   Pspec::Vector{Float64},
   Sstar::Vector{ComplexF64},
   Vm::Vector{Float64},
   V0::Vector{ComplexF64},
   drow::Vector{ComplexF64},
   nonslack::Vector{Int},
   pos_ns::Vector{Int},
   pv::Vector{Int},
   pv_pos::Vector{Int},
   slack::Int,
   Vslack::ComplexF64,
   order::Int,
   deviation::Bool,
   sh_rows::Vector{Int},
   sh_cols::Vector{Int},
   sh_coeff::Matrix{ComplexF64},
   taylor::Bool,
   eval_opts::APSLFEvaluationOptions,
   timeout_check,
)
   nbus = length(bustype)
   nn = length(nonslack)
   npv = length(pv)
   nunk = 2 * nn + npv
   shifted = !isempty(sh_rows)

   Vd = fill(zero(CDD), nbus, order + 1)
   Wd = fill(zero(CDD), nbus, order + 1)
   Qd = fill(zero(DD), npv, order)
   one_c = CDD(1.0)
   @inbounds for i = 1:nbus
      Vd[i, 1] = CDD(V0[i])
      Wd[i, 1] = one_c / Vd[i, 1]
   end

   rhs = fill(zero(DD), nunk)
   x = fill(zero(DD), nunk)
   res = fill(zero(DD), nunk)
   r64 = zeros(Float64, nunk)
   dx = zeros(Float64, nunk)
   sh_corr = fill(zero(CDD), shifted ? nbus : 0)
   half = DD(0.5)

   for n = 1:order
      _maybe_check_timeout(timeout_check, :apslf_direct_order)

      if shifted
         fill!(sh_corr, zero(CDD))
         @inbounds for e in eachindex(sh_rows)
            acc = zero(CDD)
            for k = 1:n
               acc += Vd[sh_cols[e], n-k+1] * sh_coeff[e, k]
            end
            sh_corr[sh_rows[e]] += acc
         end
      end

      @inbounds for i in nonslack
         ti = pos_ns[i]
         r_re = 2 * (ti - 1) + 1
         wc = conj(Wd[i, n])
         if _is_pq(bustype[i])
            f = wc * Sstar[i]
         elseif _is_pv(bustype[i])
            kpv = pv_pos[i]
            acc = zero(CDD)
            for m = 0:(n-2)
               acc += conj(Wd[i, (n-1-m)+1]) * Qd[kpv, m+1]
            end
            # P*conj(W^(n-1)) - j*acc
            f = wc * DD(Pspec[i]) + CDD(acc.im, -acc.re)
         else
            error("apslf_pf_pv_direct: unsupported bustype at bus $i: $(bustype[i])")
         end
         deviation && (f -= Vd[i, n] * drow[i])
         shifted && (f -= sh_corr[i])
         rhs[r_re] = f.re
         rhs[r_re+1] = f.im
      end

      @inbounds for (k, bus) in enumerate(pv)
         acc = zero(CDD)
         for m = 1:(n-1)
            acc += Vd[bus, m+1] * conj(Vd[bus, (n-m)+1])
         end
         e = -(half * acc.re)
         if n == 1
            vm = DD(Vm[bus])
            e += half * (vm * vm - (Vd[bus, 1].re * Vd[bus, 1].re + Vd[bus, 1].im * Vd[bus, 1].im))
         end
         rhs[2*nn+k] = e
      end

      _refined_solve_dd!(x, Asys, Fsys, rhs, res, r64, dx)

      @inbounds for (t, bus) in enumerate(nonslack)
         Vd[bus, n+1] = CDD(x[t], x[nn+t])
      end
      @inbounds for k = 1:npv
         Qd[k, n] = x[2*nn+k]
      end

      @inbounds for i = 1:nbus
         i == slack && continue
         acc = zero(CDD)
         for m = 1:n
            acc += Vd[i, m+1] * Wd[i, (n-m)+1]
         end
         Wd[i, n+1] = -(acc / Vd[i, 1])
      end
   end

   @inbounds for n = 1:(order+1), i = 1:nbus
      Vcoeff[i, n] = ComplexF64(Vd[i, n])
      Wcoeff[i, n] = ComplexF64(Wd[i, n])
   end
   @inbounds for n = 1:order, k = 1:npv
      Qcoeff[k, n] = Float64(Qd[k, n])
   end

   V = zeros(ComplexF64, nbus)
   V[slack] = Vslack
   @inbounds for i = 1:nbus
      i == slack && continue
      V[i] = _series_value_dd(view(Vd, i, :), taylor, eval_opts)
   end
   Qpv = zeros(Float64, npv)
   cq = Vector{CDD}(undef, order)
   @inbounds for k = 1:npv
      for n = 1:order
         cq[n] = CDD(Qd[k, n], zero(DD))
      end
      Qpv[k] = real(_series_value_dd(cq, taylor, eval_opts))
   end
   return V, Qpv
end
