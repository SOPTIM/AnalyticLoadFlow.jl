# =============================================================================
# File: src/solver_core.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
#
# Copyright 2026 SOPTIM AG
# Purpose: Implements the APSLF solver core, Padé evaluation, PV/PQ handling, Q-limit switching, and optional Newton polish helpers.
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
    solver_core.jl

APSLF (Analytic Power Series Load Flow) core routines + optional NR polishing,
plus an outer PV/Q-limit handling loop.

This file is intentionally self-contained (no dependencies on utils) to avoid
world-age / precompile issues in interactive workflows (REPL/Revise).

# Notation / conventions
- Network model: nodal admittance matrix `Y ∈ ℂ^{n×n}`
- Bus voltages:  `V ∈ ℂ^n` (phasors)
- Bus currents:  `I = Y*V`
- Complex power injection at bus `i`:
  `S_i = P_i + jQ_i = V_i * conj(I_i)`
  Vector form: `S = V .* conj(Y*V)`

- Slack bus: one reference bus with fixed voltage `V_slack`, not solved for.
- PQ bus: specified `P` and `Q` (`Pspec`, `Qspec`). Unknown: voltage magnitude/angle.
- PV bus: specified `P` and `|V|` (`Vm`). Unknown: `Q` and angle (`|V|` fixed).

# APSLF basics (high-level)
APSLF constructs an analytic continuation in an embedding parameter `s`, and
expresses voltages as a power series:

`V_i(s) = Σ_{n=0..N} V_i^(n) s^n`

and its reciprocal:

`W_i(s) = 1 / V_i(s) = Σ_{n=0..N} W_i^(n) s^n`

Power flow equations are enforced coefficient-by-coefficient.
A "germ" `V^(0)` is chosen (here often a flat start).

Evaluation at the physical operating point corresponds to `s = 1.0`.
To improve convergence radius, a Padé approximant `[L/M]` can be used
instead of simple series summation.
"""

using LinearAlgebra
using Printf
using SparseArrays

struct APSLFTimeoutError <: Exception
   where::Symbol
   elapsed_s::Float64
   timeout_s::Float64
end

function Base.showerror(io::IO, err::APSLFTimeoutError)
   print(io, "APSLF solver timeout at ", err.where, ": elapsed_s=", err.elapsed_s, ", timeout_s=", err.timeout_s)
end

_maybe_check_timeout(timeout_check, where::Symbol) = timeout_check === nothing ? nothing : timeout_check(where)

# -----------------------------------------------------------------------------
# Bus type normalization helpers (case-insensitive handling of :pv/:PV etc.)
# -----------------------------------------------------------------------------
@inline _bt_norm(bt::Symbol) = Symbol(uppercase(String(bt)))
@inline _is_pv(bt::Symbol) = (_bt_norm(bt) == :PV)
@inline _is_pq(bt::Symbol) = (_bt_norm(bt) == :PQ)
@inline _is_sl(bt::Symbol) = begin
   b = _bt_norm(bt)
   (b == :SLACK) || (b == :SL)
end

# -----------------------------------------------------------------------------
# Sparse policy helpers
# -----------------------------------------------------------------------------

"""
    maybe_sparse_Y(Y; nbus=size(Y,1), use_sparse=:auto, sparse_nbus_min=110)

Return either `Y` (as-is) or `sparse(Y)` depending on policy.

- use_sparse = :auto  => use sparse when nbus >= sparse_nbus_min
- use_sparse = true   => always convert dense -> sparse
- use_sparse = false  => never convert (but if Y already sparse, keep it sparse)
"""
function maybe_sparse_Y(
   Y::AbstractMatrix{ComplexF64};
   nbus::Int = size(Y, 1),
   use_sparse::Union{Bool,Symbol} = :auto,
   sparse_nbus_min::Int = 110,
)
   @assert size(Y, 1) == size(Y, 2)
   @assert use_sparse === :auto || use_sparse === true || use_sparse === false

   if use_sparse === false
      return Y
   elseif use_sparse === true
      return issparse(Y) ? Y : sparse(Y)
   else
      # :auto
      if nbus >= sparse_nbus_min
         return issparse(Y) ? Y : sparse(Y)
      else
         return Y
      end
   end
end

"""
    calc_injections(Y::Matrix{ComplexF64}, V::Vector{ComplexF64}) -> Vector{ComplexF64}

Compute complex power injections `S` from nodal admittance matrix `Y` and bus voltages `V`.

Definitions:
- `I = Y*V`
- `S = V .* conj(I)`

Per-bus form:
`S_i = V_i * conj( Σ_k Y_{ik} V_k )`
"""
function calc_injections(Y::AbstractMatrix{ComplexF64}, V::Vector{ComplexF64})
   I = Y * V
   return V .* conj.(I)
end

"""
    calc_injections!(Sinj, I, Y, V) -> Sinj

In-place version:
- I := Y*V
- Sinj := V .* conj(I)

`Sinj` and `I` must be preallocated ComplexF64 vectors of length nbus.
"""
function calc_injections!(
   Sinj::Vector{ComplexF64},
   I::Vector{ComplexF64},
   Y::AbstractMatrix{ComplexF64},
   V::Vector{ComplexF64},
)
   mul!(I, Y, V)
   @inbounds @simd for i in eachindex(V)
      Sinj[i] = V[i] * conj(I[i])
   end
   return Sinj
end

@inline _isfinite_complex(z) = isfinite(real(z)) && isfinite(imag(z))

@inline function _all_finite_complex(x::AbstractVector)
   @inbounds for k in eachindex(x)
      _isfinite_complex(x[k]) || return false
   end
   return true
end

@inline function _all_finite_complex_matrix(A::AbstractMatrix)
   @inbounds for j in axes(A, 2), i in axes(A, 1)
      _isfinite_complex(A[i, j]) || return false
   end
   return true
end

@inline function _first_nonfinite_complex(x::AbstractVector)
   @inbounds for k in eachindex(x)
      _isfinite_complex(x[k]) || return (k, x[k])
   end
   return nothing
end

@inline function _first_nonfinite_complex_matrix(A::AbstractMatrix)
   @inbounds for j in axes(A, 2), i in axes(A, 1)
      _isfinite_complex(A[i, j]) || return ((i, j), A[i, j])
   end
   return nothing
end

@inline function _assert_finite_vec(x::AbstractVector, name::AbstractString)
   bad = _first_nonfinite_complex(x)
   bad === nothing || error("Non-finite entry in $name at index $(bad[1]): $(bad[2])")
   return nothing
end

@inline function _assert_finite_mat(A::AbstractMatrix, name::AbstractString)
   bad = _first_nonfinite_complex_matrix(A)
   bad === nothing || error("Non-finite entry in $name at ($(bad[1][1]),$(bad[1][2])): $(bad[2])")
   return nothing
end

@inline function _assert_finite_coeff_matrix(A::AbstractMatrix, name::AbstractString)
   bad = _first_nonfinite_complex_matrix(A)
   bad === nothing || error("First non-finite coefficient in $name at bus=$(bad[1][1]), order=$(bad[1][2] - 1): $(bad[2])")
   return nothing
end

@inline function _series_eval_finite_prefix(c::AbstractVector{ComplexF64}; s::ComplexF64 = 1.0 + 0.0im)
   last_finite = 0
   @inbounds for k in eachindex(c)
      _isfinite_complex(c[k]) || break
      last_finite = k
   end

   last_finite == 0 && return ComplexF64(NaN, NaN)

   acc = 0.0 + 0.0im
   @inbounds for k = last_finite:-1:1
      acc = acc * s + c[k]
   end
   return acc
end

@inline function _is_recoverable_linear_solve_error(err)
   return err isa SingularException || err isa DomainError || err isa ArgumentError
end

@inline function _is_nonfinite_inner_error(err)::Bool
   return (err isa ErrorException) && occursin("Non-finite", sprint(showerror, err))
end

"""
    APSLFEvaluationOptions(; mode=:auto, auto=(...), pade=(...))

Configuration for APSLF series evaluation.

`mode` selects how voltage series are evaluated:

- `:taylor`: direct series summation
- `:pade`: Padé evaluation
- `:auto`: choose Taylor or Padé based on simple coefficient checks

The `auto` and `pade` named tuples contain conservative numerical thresholds
and fallback settings used by `evaluate_series`.
"""
struct APSLFEvaluationOptions
   mode::Symbol
   auto::NamedTuple
   pade::NamedTuple
end

function APSLFEvaluationOptions(;
   mode::Symbol = :auto,
   auto::NamedTuple = (
      enabled = true,
      prefer_pade = true,
      min_order_before_check = 6,
      coeff_decay_tol = 1.0e-6,
      partial_sum_tol = 1.0e-8,
      max_last_coeff_ratio = 0.85,
      pade_min_order = 6,
      pade_max_order = 20,
      pade_tolerance = 1.0e-7,
   ),
   pade::NamedTuple = (order_strategy = :balanced, fallback_to_taylor = true),
)
   @assert mode in (:auto, :taylor, :pade) "evaluation mode must be :auto, :taylor, or :pade"
   return APSLFEvaluationOptions(mode, auto, pade)
end

"""
    is_taylor_good(coeffs, opts) -> Bool

Return `true` when the tail of a voltage power series appears sufficiently
small for direct Taylor summation under the supplied `opts` thresholds.

This is a lightweight heuristic used by automatic series evaluation. It is not
a formal convergence proof.
"""
@inline function is_taylor_good(coeffs::AbstractVector{ComplexF64}, opts::NamedTuple)
   N = length(coeffs) - 1
   N < opts.min_order_before_check && return false
   last = norm(coeffs[end])
   base = max(norm(coeffs[1]), eps())
   rel = last / base
   if rel >= opts.coeff_decay_tol || rel > opts.max_last_coeff_ratio
      return false
   end
   return true
end

"""
    evaluate_series(coeffs, options; collect_logs=true) -> NamedTuple

Evaluate a complex power series using Taylor summation, Padé approximation, or
automatic selection according to `options`.

Returns a named tuple containing at least:

- `voltage`
- `method_used`
- `reason`
- `pade_triggered`
- `logs`
"""
function evaluate_series(
   coeffs::AbstractVector{ComplexF64},
   options::APSLFEvaluationOptions;
   collect_logs::Bool = true,
)::NamedTuple
   logs = collect_logs ? String["APSLF Evaluation", "----------------", "mode_requested     : $(options.mode)"] : String[]
   N = length(coeffs) - 1
   do_pade = false
   reason = :forced_taylor
   method_used = :taylor
   if options.mode == :taylor
      do_pade = false
      reason = :forced_taylor
   elseif options.mode == :pade
      do_pade = true
      reason = :forced_pade
   else
      taylor_ok = options.auto.enabled && is_taylor_good(coeffs, options.auto)
      do_pade = !taylor_ok && options.auto.prefer_pade
      reason = do_pade ? :weak_taylor_convergence : :taylor_good
   end

   v = 0.0 + 0.0im
   pade_triggered = (options.mode == :auto && do_pade)
   if do_pade
      M = clamp(div(N, 2), 0, N)
      L = N - M
      try
         v = _pade_eval_or_series(coeffs, L, M; s = 1.0 + 0.0im)
         method_used = :pade
      catch err
         if options.pade.fallback_to_taylor && _is_recoverable_linear_solve_error(err)
            v = sum(coeffs)
            method_used = :taylor
            reason = :pade_failed_fallback_to_taylor
            collect_logs && push!(logs, "warning            : pade_failed_fallback_to_taylor")
         else
            rethrow()
         end
      end
   else
      v = sum(coeffs)
      method_used = :taylor
   end
   collect_logs && push!(logs, "method_used        : $(method_used)")
   collect_logs && push!(logs, "reason             : $(reason)")
   collect_logs && push!(logs, "pade_triggered     : $(pade_triggered)")
   return (voltage = v, method_used = method_used, reason = reason, pade_triggered = pade_triggered, logs = logs)
end

@inline function _pade_eval_or_series(c::AbstractVector{ComplexF64}, L::Int, M::Int; s::ComplexF64 = 1.0 + 0.0im)
   if _all_finite_complex(c)
      try
         v = pade_eval(c, L, M; s = s)
         _isfinite_complex(v) && return v
      catch err
         if !_is_recoverable_linear_solve_error(err)
            rethrow()
         end
      end
   end

   return _series_eval_finite_prefix(c; s = s)
end

"""
    pade_eval(c, L, M; s=1+0im) -> ComplexF64

Evaluate the rational Padé approximant `[L/M]` of a power series at point `s` (default `s=1`).

# Arguments
- `c::AbstractVector{ComplexF64}`: Power series coefficients `[c₀, c₁, ..., cₙ]`
- `L::Int`: Degree of numerator polynomial (≥ 0)
- `M::Int`: Degree of denominator polynomial (≥ 0)
- `s::ComplexF64 = 1+0im`: Evaluation point in complex plane

# Returns
- `ComplexF64`: Value of Padé approximant P(s)/Q(s) at the given point

# Mathematical Formulation
Given a power series:
```
f(s) = c₀ + c₁s + c₂s² + ... + cₙsⁿ
```

Constructs and evaluates the Padé [L/M] rational approximation:
```
f(s) ≈ P(s)/Q(s) = (a₀ + a₁s + ... + aₗsᴸ) / (1 + b₁s + ... + bₘsᴹ)
```

# Algorithm Steps
1. **Handle M=0 case**: Pure polynomial evaluation using Horner's method
2. **Build linear system**: Solve Toeplitz-like system for denominator coefficients `b₁, ..., bₘ`
3. **Compute numerator**: Calculate `a₀, ..., aₗ` via convolution with denominator
4. **Horner evaluation**: Evaluate both polynomials efficiently at point `s`
5. **Return ratio**: `P(s) / Q(s)`

# Linear System for Denominator
The denominator coefficients satisfy:
```
Σⱼ₌₁ᴹ bⱼ c_{L+i-j} = -c_{L+i}  for i = 1, ..., M
```
This ensures the rational function matches the series expansion through order L+M.

# Numerator Construction
Once denominator is known:
```
aₙ = cₙ + Σⱼ₌₁^{min(n,M)} bⱼ c_{n-j}  for n = 0, ..., L
```

# Constraints and Special Cases
- Requires `length(c) - 1 ≥ L + M` (sufficient series coefficients)
- **M = 0**: Returns pure polynomial `P(s) = Σₙ₌₀ᴸ cₙ sⁿ` (no rational approximation)
- **Denominator normalization**: `Q(s) = 1 + b₁s + ... + bₘsᴹ` (constant term = 1)
- **Horner evaluation**: O(L + M) operations for numerical stability

# Series Coefficient Indexing
Input vector `c` uses 1-based Julia indexing:
- `c[1]` → coefficient c₀ (constant term)
- `c[n+1]` → coefficient cₙ (sⁿ term)

# Applications in APSLF
- **Series acceleration**: Better approximation quality with fewer terms
- **Convergence extension**: Padé extends radius of convergence beyond direct summation
- **Voltage evaluation**: `V(s=1)` via Padé of voltage series coefficients
- **Stability analysis**: Denominator poles indicate system stability margins

# Numerical Considerations
- Uses direct linear solve for denominator system (efficient for moderate M)
- Horner's method provides numerical stability for polynomial evaluation
- Handles general complex coefficients and evaluation points
- Avoid evaluation near denominator zeros (poles of approximant)

# Example
```julia
# Series: f(s) = 1 + 0.5s + 0.25s² + 0.125s³ + ...  (geometric series)
c = ComplexF64[1.0, 0.5, 0.25, 0.125, 0.0625]

# Evaluate [2/1] Padé at s=0.8
result = pade_eval(c, 2, 1; s=0.8+0.0im)

# Compare to exact: f(0.8) = 1/(1-0.5*0.8) = 1/0.6 ≈ 1.667
# Padé should be very accurate for this geometric series
```

# Performance Notes
- Computational cost: O(M²) for linear system + O(L+M) for evaluation
- More efficient than `pade_build` + separate evaluation for single points
- Recomputes denominator coefficients each call (use `pade_build` for multiple evaluations)

# See Also
- `pade_build`: Build numerator/denominator coefficients separately
- `pade_poles`: Extract poles of Padé approximant for stability analysis
- `apslf_pq`: Main APSLF solver using Padé evaluation for voltage series
"""
function pade_eval(c::AbstractVector{ComplexF64}, L::Int, M::Int; s::ComplexF64 = 1.0 + 0.0im)
   N = length(c) - 1
   @assert N >= L + M "Need at least N >= L+M coefficients for [L/M] Padé."
   _assert_finite_vec(c, "Padé coefficient vector c")

   if M == 0
      # [L/0] => pure polynomial, denominator = 1
      num = 0.0 + 0.0im
      for n = L:-1:0
         num = num * s + c[n+1]
      end
      return num
   end

   # Build Toeplitz-like system A*b = rhs for denominator coefficients b
   # such that series of (num/den) matches series c up to order L+M.
   A = zeros(ComplexF64, M, M)
   rhs = zeros(ComplexF64, M)
   for row = 1:M
      n = L + row
      rhs[row] = -c[n+1]                 # -c_{L+row}
      for j = 1:M
         A[row, j] = c[(n-j)+1]       # c_{L+row-j}
      end
   end

   _assert_finite_mat(A, "Padé matrix A")
   _assert_finite_vec(rhs, "Padé rhs b")

   b = A \ rhs  # b1..bM
   _assert_finite_vec(b, "Padé denominator q")

   # Numerator coefficients a_0..a_L:
   #   a_n = c_n + Σ_{j=1..min(n,M)} b_j c_{n-j}
   a = zeros(ComplexF64, L + 1)  # a0..aL
   for n = 0:L
      acc = c[n+1]
      for j = 1:min(n, M)
         acc += b[j] * c[(n-j)+1]
      end
      a[n+1] = acc
   end

   # Evaluate numerator polynomial a(s) via Horner
   num = 0.0 + 0.0im
   for n = L:-1:0
      num = num * s + a[n+1]
   end

   # Evaluate denominator polynomial 1 + b1*s + ... + bM*s^M via Horner
   den = b[M]
   for j = (M-1):-1:1
      den = den * s + b[j]
   end
   den = den * s + (1.0 + 0.0im)

   return num / den
end


"""
    pade_build(c, L, M) -> (a, b)

Build Padé [L/M] numerator and denominator coefficients from power series coefficients.

# Arguments
- `c::AbstractVector{ComplexF64}`: Power series coefficients `[c₀, c₁, ..., cₙ]`
- `L::Int`: Degree of numerator polynomial (≥ 0)
- `M::Int`: Degree of denominator polynomial (≥ 0)

# Returns
- `a::Vector{ComplexF64}`: Numerator coefficients `[a₀, a₁, ..., aₗ]` (length L+1)
- `b::Vector{ComplexF64}`: Denominator coefficients `[b₀, b₁, ..., bₘ]` (length M+1)

where `b₀ = 1` (normalized form).

# Mathematical Formulation
Given a power series:
```
f(s) = c₀ + c₁s + c₂s² + ... + cₙsⁿ
```

Constructs the Padé [L/M] rational approximation:
```
f(s) ≈ P(s)/Q(s) = (a₀ + a₁s + ... + aₗsᴸ) / (1 + b₁s + ... + bₘsᴹ)
```

# Algorithm
1. **Denominator coefficients**: Solve the linear system for `b₁, ..., bₘ`:
   ```
   Σⱼ₌₁ᴹ bⱼ c_{L+i-j} = -c_{L+i}  for i = 1, ..., M
   ```
   This ensures the series expansion of P(s)/Q(s) matches the original series
   through order L+M.

2. **Numerator coefficients**: Compute via convolution:
   ```
   aₙ = cₙ + Σⱼ₌₁^{min(n,M)} bⱼ c_{n-j}  for n = 0, ..., L
   ```

# Constraints
- Requires `length(c) - 1 ≥ L + M` (i.e., `N ≥ L + M`)
- For `M = 0`: returns pure polynomial (denominator = [1])
- Denominator is normalized so `b₀ = 1`

# Series Coefficient Storage
Input vector `c` uses 1-based Julia indexing:
- `c[1]` contains `c₀` (constant term)
- `c[n+1]` contains `cₙ` (coefficient of sⁿ)

# Applications in APSLF
- Voltage series evaluation: improves convergence radius over direct summation
- Stability analysis: denominator roots give Padé pole locations
- Series acceleration: better approximation with fewer terms

# Example
```julia
# Power series: f(s) = 1 + 2s + 3s² + 4s³ + ...
c = ComplexF64[1, 2, 3, 4, 5, 6]  # c₀ through c₅

# Build [2/2] Padé approximant
a, b = pade_build(c, 2, 2)

# Result: P(s)/Q(s) where
# P(s) = a[1] + a[2]s + a[3]s²
# Q(s) = b[1] + b[2]s + b[3]s² = 1 + b[2]s + b[3]s²
```

# See Also
- `pade_eval`: Evaluate Padé approximant at specific point
- `pade_poles`: Extract denominator poles for stability analysis
- `poly_roots`: General polynomial root finding
"""
function pade_build(c::AbstractVector{ComplexF64}, L::Int, M::Int)
   N = length(c) - 1
   @assert N >= L + M "Need at least N >= L+M coefficients for [L/M] Padé."

   A = zeros(ComplexF64, M, M)
   rhs = zeros(ComplexF64, M)
   for row = 1:M
      n = L + row
      rhs[row] = -c[n+1]                 # -c_{L+row}
      for j = 1:M
         A[row, j] = c[(n-j)+1]       # c_{L+row-j}
      end
   end

   b1 = A \ rhs  # b1..bM

   a = zeros(ComplexF64, L + 1)
   for n = 0:L
      acc = c[n+1]
      for j = 1:min(n, M)
         acc += b1[j] * c[(n-j)+1]
      end
      a[n+1] = acc
   end

   b = Vector{ComplexF64}(undef, M + 1)
   b[1] = 1.0 + 0.0im
   @inbounds for j = 1:M
      b[j+1] = b1[j]
   end

   return a, b
end

"""
    poly_roots(coeffs) -> Vector{ComplexF64}

Compute roots of a polynomial using companion-matrix eigenvalues.

# Arguments
- `coeffs::AbstractVector{ComplexF64}`: Polynomial coefficients in ascending order of powers

# Returns
- `Vector{ComplexF64}`: All roots of the polynomial (length = degree)

# Mathematical Formulation
Given coefficients `[c₀, c₁, ..., cₘ]`, finds roots of:
```
q(s) = c₀ + c₁s + c₂s² + ... + cₘsᴹ
```


# Algorithm
Uses the companion matrix method:
1. **Normalization**: Convert to monic form `sᴹ + aₘ₋₁ sᴹ⁻¹ + ... + a₀`
2. **Companion matrix**: Construct the M×M companion matrix
3. **Eigenvalues**: Roots are the eigenvalues of that matrix

# Coefficient Storage Convention
- `coeffs[1]` → constant term `c₀`
- `coeffs[k+1]` → coefficient `cₖ` of `sᵏ`
- `coeffs[end]` → leading coefficient `cₘ` (must be nonzero)

# Constraints
- Requires degree ≥ 1 (`length(coeffs) ≥ 2`)
- Leading coefficient must be nonzero

# Numerical Considerations
- Companion-matrix root finding is general but can be ill-conditioned for higher degrees
  or poorly scaled coefficients.
- Accuracy depends strongly on polynomial conditioning and coefficient scaling.
- For degrees above ~20–40 (problem-dependent), expect potential sensitivity; rescaling
  or using specialized root finders may be beneficial.

# Applications in APSLF
- Extracting Padé denominator roots (poles) for heuristic diagnostics.

# See Also
- `pade_poles`
- `pade_build`
"""
function poly_roots(coeffs::AbstractVector{ComplexF64})
   M = length(coeffs) - 1
   @assert M >= 1 "Need degree >= 1."
   qM = coeffs[end]
   @assert abs(qM) > 0 "Leading coefficient must be nonzero."

   # Make monic: s^M + a_{M-1} s^{M-1} + ... + a0
   a = coeffs[1:(end-1)] ./ qM  # a0..a_{M-1}

   C = zeros(ComplexF64, M, M)
   @inbounds for i = 2:M
      C[i, i-1] = 1.0 + 0.0im
   end
   @inbounds for i = 1:M
      C[i, M] = -a[i]  # last column: -a0..-a_{M-1}
   end
   return eigvals(C)
end

"""
    pade_poles(c, L, M) -> Vector{ComplexF64}

Extract poles of the Padé [L/M] approximant, i.e. the roots of its denominator polynomial.

# Arguments
- `c::AbstractVector{ComplexF64}`: Power series coefficients `[c₀, c₁, ..., cₙ]`
- `L::Int`: Numerator degree (≥ 0)
- `M::Int`: Denominator degree (≥ 0)

# Returns
- `Vector{ComplexF64}`: The `M` poles (roots of the Padé denominator), length `M`.
  For `M = 0`, returns an empty vector.

# Background
The Padé approximant is constructed in the normalized form:
```
f(s) ≈ P(s)/Q(s) = (a₀ + a₁ s + ... + a_L s^L) / (1 + b₁ s + ... + b_M s^M)
```

So the denominator is normalized with constant term `b₀ = 1` (not monic in the
leading coefficient). Root finding (`poly_roots`) internally normalizes to monic form.

# Notes on Interpretation (APSLF context)
- Poles are primarily a *heuristic diagnostic* derived from the Padé denominator.
- A small distance `|pole - 1|` may correlate with reduced analytic continuation margin,
  but it is not a universal or rigorous stability certificate.
- Spurious pole-zero pairs can occur, especially at higher orders or with noisy coefficients.

# Constraints
- Requires `length(c) - 1 ≥ L + M`.

# See Also
- `pade_build`
- `poly_roots`
- `stability_from_Vcoeff`
"""
function pade_poles(c::AbstractVector{ComplexF64}, L::Int, M::Int)
   _, b = pade_build(c, L, M)   # b0..bM with b0=1
   return poly_roots(b)
end

"""
    stability_from_Vcoeff(Vcoeff; slack=1, order=size(Vcoeff,2)-1, max_buses=0, critical_poles=:auto)
        -> NamedTuple

Compute a heuristic APSLF margin indicator from Padé poles of per-bus voltage series.

# Arguments
- `Vcoeff::Matrix{ComplexF64}`: Voltage series coefficients (nbus × order+1)
  - Row `i` contains `[V_i^(0), V_i^(1), ..., V_i^(order)]`
- `slack::Int = 1`: Slack bus index (excluded from pole analysis)
- `order::Int = size(Vcoeff,2)-1`: Series order used for Padé construction
- `max_buses::Int = 0`: If >0, analyze only the first `max_buses` non-slack buses (0 = all)
- `critical_poles::Union{Symbol,Int} = :auto`: Number of critical poles to report.
  - `:auto`: scale with network size as `ceil(sqrt(n_non_slack))`
  - `k::Int`: return the `k` poles with smallest `|pole - 1|`

# Returns
NamedTuple with:
- `dmin::Float64`: Minimum distance `min |pole - 1|` over all analyzed buses/poles
- `pole::ComplexF64`: Pole that attains `dmin`
- `bus::Int`: Bus index where that pole was found
- `L::Int`, `M::Int`: Degrees used for Padé [L/M]
- `critical::Vector{NamedTuple}`: Critical poles sorted by ascending distance:
  each entry is `(bus::Int, pole::ComplexF64, distance::Float64)`

# Definition
For each non-slack bus `i`, form a Padé approximant of:
```
V_i(s) = Σ_{n=0}^{order} V_i^(n) s^n
```
using:
- `M = order ÷ 2`
- `L = order - M`

Then compute:
```
dmin = min_{i,k} |p_{i,k} - 1|
```
where `p_{i,k}` are the Padé denominator roots (poles).

# Interpretation (Important)
- This is a *heuristic analytic-continuation margin* based on Padé poles.
- A smaller `dmin` can indicate that the Padé approximation has a pole closer to the
  physical evaluation point `s = 1`, which may correlate with reduced continuation margin.
- It is **not** a rigorous voltage-stability proof and should not be equated with
  eigenvalue-based small-signal stability or classic V–Q/PV margin certificates.
- Spurious poles may occur; compare multiple [L/M] choices or orders if you rely on it.

# Notes
- `max_buses` is purely a speed knob; it may miss the globally closest pole.

# See Also
- `pade_poles`
- `pade_build`
"""
function stability_from_Vcoeff(
   Vcoeff::Matrix{ComplexF64};
   slack::Int = 1,
   order::Int = size(Vcoeff, 2) - 1,
   max_buses::Int = 0,
   critical_poles::Union{Symbol,Int} = :auto,
)
   Nv = order
   Mv = Nv ÷ 2
   Lv = Nv - Mv

   nbus = size(Vcoeff, 1)
   buses = [i for i = 1:nbus if i != slack]
   if max_buses > 0 && length(buses) > max_buses
      buses = buses[1:max_buses]
   end

   best_d = Inf
   best_p = 0.0 + 0.0im
   best_bus = 0
   candidates = NamedTuple{(:bus, :pole, :distance),Tuple{Int,ComplexF64,Float64}}[]

   @inbounds for i in buses
      cV = @view Vcoeff[i, :]  # V^(0..order), no alloc
      poles = pade_poles(cV, Lv, Mv)
      for p in poles
         d = abs(p - (1.0 + 0.0im))
         push!(candidates, (bus = i, pole = p, distance = d))
         if d < best_d
            best_d = d
            best_p = p
            best_bus = i
         end
      end
   end

   if critical_poles === :auto
      k = isempty(buses) ? 0 : ceil(Int, sqrt(length(buses)))
   elseif critical_poles isa Int
      critical_poles >= 0 || throw(ArgumentError("critical_poles must be :auto or a non-negative Int."))
      k = critical_poles
   else
      throw(ArgumentError("critical_poles must be :auto or a non-negative Int."))
   end

   sort!(candidates, by = x -> x.distance)
   k = min(k, length(candidates))
   critical = candidates[1:k]

   return (dmin = best_d, pole = best_p, bus = best_bus, L = Lv, M = Mv, critical = critical)
end

function _handle_deprecated_apslf_germ_kwargs(kwargs; context::AbstractString, strict::Bool = true)
   for (key, value) in pairs(kwargs)
      if key === :flatstart
         @warn "`flatstart` is deprecated for the APSLF solver. APSLF always uses the canonical analytic germ V(s=0)=1∠0. This option will be removed." context =
            context flatstart = value
      elseif key === :V0_germ
         @warn "`V0_germ` is not a Newton start value. Passing a custom APSLF germ changes the analytic embedding and is no longer supported in the main APSLF solver path." context =
            context
      else
         strict && throw(ArgumentError("Unknown keyword argument `$key` in $context."))
      end
   end
   return nothing
end

"""
    apslf_germ(Y, nonslack, slack, Vslack, germ; F=nothing) -> Vector{ComplexF64}

Return the APSLF germ `V^(0)` for the non-slack buses in the order of `nonslack`.

- `germ = :flat`   → canonical flat germ `1∠0`. Exact at order 0 only if the
  constant matrix has zero row sums (pure series network) and `Vslack = 1`.
- `germ = :noload` → solution of the linear no-load problem
  `Y_red V^(0) = -Y[red, slack] Vslack`. This is the exact `s = 0` state for
  any `Y`: line shunts, off-nominal transformer ratios, phase shifters (PST)
  and `Vslack ≠ 1` are all absorbed into the germ, and the recursion runs
  unchanged with a non-uniform `W^(0) = 1 ./ V^(0)` (theory Section 6.5,
  variant 2). Costs one additional solve with the same factorization.
- `germ = :deviation` → flat germ `Vslack·1` combined with the deviation
  embedding `Y(s) = Y0 + s (Y - Y0)`, `Y0 = Y - diag(Y·1)` (theory Section
  6.5, variant 1). `Y0` has zero row sums by construction, so the constant
  germ is exact at order 0 for any `Y`, and the deviation `diag(Y·1)` (bus
  shunts, transformer and PST row sums) is ramped up with `s` on the
  right-hand side of the recursion. The germ stays at nominal voltage even
  when the no-load state is far from the operating point (large networks
  with strong Ferranti rise), which usually gives the larger convergence
  radius of the two variants.

`F` may be a factorization of `Y[nonslack, nonslack]` to avoid refactoring.
"""
function apslf_germ(
   Y::AbstractMatrix{ComplexF64},
   nonslack::Vector{Int},
   slack::Int,
   Vslack::ComplexF64,
   germ::Symbol;
   F = nothing,
)
   germ in (:flat, :noload, :deviation) || throw(ArgumentError("germ must be :noload, :deviation or :flat, got :$(germ)"))
   n = length(nonslack)
   germ == :flat && return fill(1.0 + 0.0im, n)
   germ == :deviation && return fill(Vslack, n)
   Yslack = Vector{ComplexF64}(Y[nonslack, slack])
   return apslf_noload_germ(F === nothing ? lu(issparse(Y) ? sparse(Y[nonslack, nonslack]) : Matrix(Y[nonslack, nonslack])) : F, Yslack, Vslack)
end

"""
    apslf_row_sums(Y) -> Vector{ComplexF64}

Row sums `Y·1` of the bus admittance matrix; the diagonal deviation of the
`:deviation` embedding (see [`apslf_germ`](@ref)).
"""
apslf_row_sums(Y::AbstractMatrix{ComplexF64}) = Vector{ComplexF64}(vec(sum(Y, dims = 2)))

"""
    apslf_noload_germ(F, Yslack, Vslack) -> Vector{ComplexF64}

No-load germ from a factorization `F` of the reduced matrix and the reduced
slack column `Yslack = Y[nonslack, slack]`: solves `Y_red V0 = -Yslack * Vslack`.
"""
function apslf_noload_germ(F, Yslack::Vector{ComplexF64}, Vslack::ComplexF64)
   V0 = F \ (-Yslack .* Vslack)
   all(isfinite, V0) || error("APSLF no-load germ is not finite; the reduced Y-bus may be singular.")
   vmin = minimum(abs, V0)
   vmin > 1e-6 || error("APSLF no-load germ has a (near-)zero bus voltage (min |V0| = $(vmin)); W = 1/V is undefined.")
   return V0
end

function _warn_degenerate_q_limits(
   bustype::Vector{Symbol},
   Qmin::Vector{Float64},
   Qmax::Vector{Float64};
   atol::Float64 = 1e-12,
)
   buses_equal = Int[]
   @inbounds for i in eachindex(bustype)
      _is_pv(bustype[i]) || continue
      qmin = Qmin[i]
      qmax = Qmax[i]
      if isfinite(qmin) && isfinite(qmax)
         if qmin > qmax + atol
            throw(ArgumentError("Invalid reactive limits at bus $i: Qmin ($qmin) > Qmax ($qmax)."))
         end
         if abs(qmax - qmin) <= atol
            push!(buses_equal, i)
         end
      end
   end

   if !isempty(buses_equal)
      @warn "Detected PV buses with Qmin ≈ Qmax. These buses are effectively fixed-Q and may destabilize direct APSLF convergence." buses =
         buses_equal
   end
   return nothing
end

function _canonicalize_bustype_lower(bustype::Vector{Symbol})
   bt = similar(bustype)
   @inbounds for i in eachindex(bustype)
      b = bustype[i]
      if _is_pv(b)
         bt[i] = :pv
      elseif _is_pq(b)
         bt[i] = :pq
      elseif _is_sl(b)
         bt[i] = :slack
      else
         bt[i] = b
      end
   end
   return bt
end

function _demote_degenerate_pv_buses!(
   bt::Vector{Symbol},
   Q::Vector{Float64},
   Qmin::Vector{Float64},
   Qmax::Vector{Float64};
   atol::Float64 = 1e-12,
)
   buses_demoted = Int[]
   @inbounds for i in eachindex(bt)
      _is_pv(bt[i]) || continue
      qmin = Qmin[i]
      qmax = Qmax[i]
      if isfinite(qmin) && isfinite(qmax) && abs(qmax - qmin) <= atol
         bt[i] = :pq
         Q[i] = qmin
         push!(buses_demoted, i)
      end
   end
   return buses_demoted
end

"""
    apslf_pq(Y, S; slack=1, Vslack=1+0im, order=24, use_pade=true)
        -> (V, Vcoeff, Wcoeff)

APSLF for a PQ-only network with a single slack bus fixed at `Vslack`.

Model:
- Slack bus voltage is fixed (not part of unknowns)
- All other buses treated as PQ (specified `S = P + jQ`)

Coefficient recursion (standard PQ embedding):
Let `V(s) = Σ V^(n) s^n` and `W(s) = 1/V(s) = Σ W^(n) s^n`.

For `n ≥ 1`:
`Y_red * V_red^(n) = conj(S_red) .* conj(W_red^(n-1))`

where "red" excludes the slack bus.

Inverse series recursion for `W(s)=1/V(s)` with canonical germ `v0 = V^(0)`:
- `w_0 = 1/v0`
- `w_n = -(1/v0) * Σ_{m=1..n} v_m w_{n-m}`

Evaluation at `s=1`:
- by direct series sum `Σ V^(n)` or
- by Padé `[L/M]` approximation of the voltage series coefficients.
"""
function apslf_pq(
   Y::AbstractMatrix{ComplexF64},
   S::Vector{ComplexF64};
   slack::Int = 1,
   Vslack::ComplexF64 = 1.0 + 0.0im,
   order::Int = 24,
   use_pade::Bool = true,
   germ::Symbol = :deviation,
   evaluation_options::Union{Nothing,APSLFEvaluationOptions} = nothing,
   timeout_check = nothing,
   kwargs...,
)
   _handle_deprecated_apslf_germ_kwargs(kwargs; context = "apslf_pq")

   nbus = size(Y, 1)
   @assert size(Y, 2) == nbus
   @assert length(S) == nbus

   pq = [i for i = 1:nbus if i != slack]
   npq = length(pq)

   # Reduced linear system for coefficient solves:
   # If Y is sparse, keep the reduced matrix sparse (UMFPACK LU)
   Yred = issparse(Y) ? sparse(Y[pq, pq]) : Matrix(Y[pq, pq])
   # :deviation embedding: constant matrix Y0 = Y - diag(Y·1), deviation d = Y·1 on the RHS
   drow = germ == :deviation ? apslf_row_sums(Y)[pq] : ComplexF64[]
   germ == :deviation && (Yred = Yred - Diagonal(drow))
   F = lu(Yred)

   # Coefficient storage (reduced ordering)
   # Vcoeff[:, n+1] holds V^(n), because column 1 is n=0.
   Vcoeff = zeros(ComplexF64, npq, order + 1)
   Wcoeff = zeros(ComplexF64, npq, order + 1)

   # ---- APSLF germ V^(0): canonical flat 1∠0, or the no-load solution (germ=:noload).
   V0 = apslf_germ(Y, pq, slack, Vslack, germ; F = F)

   Vcoeff[:, 1] .= V0
   Wcoeff[:, 1] .= 1.0 ./ V0

   # conj(S) is used in the standard APSLF coefficient equation for PQ
   Sstar = conj.(S[pq])

   for n = 1:order
      _maybe_check_timeout(timeout_check, :apslf_pq_order)
      # At order n:
      #   rhs = conj(S) .* conj(W^(n-1))   (reflection condition, theory Sec. 2.4/4.2)
      # Note: Wcoeff[:, n] corresponds to W^(n-1) due to shift.
      rhs = Sstar .* conj.(Wcoeff[:, n])
      germ == :deviation && (rhs .-= drow .* @view(Vcoeff[:, n]))   # - d .* V^(n-1)
      Vcoeff[:, n+1] = F \ rhs

      # Inverse series recursion:
      #   w_n = -(1/v0) * Σ_{m=1..n} v_m w_{n-m}
      @inbounds for i = 1:npq
         acc = 0.0 + 0.0im
         for m = 1:n
            acc += Vcoeff[i, m+1] * Wcoeff[i, (n-m)+1]
         end
         Wcoeff[i, n+1] = -acc / Vcoeff[i, 1]
      end
   end

   # Assemble full voltage vector (PF ordering)
   V = zeros(ComplexF64, nbus)
   V[slack] = Vslack

   # Evaluate at s=1
   eval_opts = isnothing(evaluation_options) ? APSLFEvaluationOptions(mode = use_pade ? :pade : :taylor) : evaluation_options
   @inbounds for (idx, bus) in enumerate(pq)
      c = @view Vcoeff[idx, :]
      V[bus] = evaluate_series(c, eval_opts).voltage
   end

   return V, Vcoeff, Wcoeff
end


# =============================================================================
# APSLF PQ workspace (LU reuse)
# =============================================================================
"""
    APSLFPQWorkspace{TY, TF}

Pre-allocated workspace for efficient repeated APSLF PQ-only power flow solutions.

# Type Parameters
- `TY <: AbstractMatrix{ComplexF64}`: Type of reduced admittance matrix (sparse or dense)
- `TF <: LinearAlgebra.Factorization`: Type of LU factorization of reduced system

# Fields
- `nbus::Int`: Total number of buses in the system
- `slack::Int`: Slack bus index (fixed voltage reference)
- `order::Int`: APSLF series expansion order

## Bus Indexing and Reduction
- `pq::Vector{Int}`: Non-slack bus indices in original PF ordering
- `npq::Int`: Number of non-slack buses (`length(pq)`)
- `Yred::TY`: Reduced admittance matrix `Y[pq, pq]` (excludes slack bus)
- `F::TF`: LU factorization of `Yred` for efficient linear solves

## Coefficient Storage (Reduced Ordering)
- `Vcoeff::Matrix{ComplexF64}`: Voltage series coefficients (npq × order+1)
  - Row `k` corresponds to bus `pq[k]`
  - Column `n+1` stores coefficient `V^(n)` for all non-slack buses
- `Wcoeff::Matrix{ComplexF64}`: Inverse series coefficients (npq × order+1)
  - `W(s) = 1/V(s)` where `W^(n)` stored in column `n+1`

## Work Arrays (Reduced Ordering)
- `rhs::Vector{ComplexF64}`: Right-hand side for linear system (length npq)
- `Sstar::Vector{ComplexF64}`: Conjugated power specifications (length npq)

## Output Buffer (Full PF Ordering)
- `Vfull::Vector{ComplexF64}`: Complete voltage solution vector (length nbus)

# Purpose and Usage
This workspace eliminates repeated allocations and expensive LU factorizations
when solving multiple APSLF PQ power flows with the same system structure:

1. **One-time setup**: Build workspace with `build_apslf_pq_workspace(Y, ...)`
2. **Repeated solves**: Call `apslf_pq_solve!(ws, S, ...)` multiple times
3. **Efficiency gains**: Reuses `F = lu(Yred)` and all coefficient storage

# Mathematical Context
The workspace supports the APSLF PQ recursion:
```
Y_red * V^(n) = conj(S_red) .* conj(W^(n-1))    for n = 1, 2, ..., order
```

where:
- `V^(0) = germ` (typically flat start or provided initialization)
- `W^(n)` follows inverse series recursion: `W^(n) = -(1/V^(0)) * Σ_{m=1}^n V^(m) W^(n-m)`

# Storage Conventions
- **Coefficient indexing**: Column `n+1` stores order-`n` coefficient (0-based math, 1-based Julia)
- **Reduced ordering**: Non-slack buses only, indexed by position in `pq` vector
- **Full ordering**: Complete nbus-length vectors using original bus numbering

# Performance Benefits
- **LU reuse**: Most expensive operation (factorization) done once
- **Memory reuse**: All work arrays pre-allocated, eliminating GC pressure
- **Cache locality**: Contiguous coefficient storage improves memory access patterns
- **Type stability**: Concrete types `TY` and `TF` enable compiler optimizations

# Example
```julia
# Setup: build workspace once
Y = create_admittance_matrix(...)
ws = build_apslf_pq_workspace(Y; slack=1, order=24)

# Repeated solves: efficient iterations
for load_level in [0.8, 0.9, 1.0, 1.1, 1.2]
    S = base_load * load_level
    V, Vcoeff, Wcoeff = apslf_pq_solve!(ws, S; use_pade=true)

    # Process results...
    voltages[load_level] = copy(V)  # copy since ws.Vfull is reused
end

# Workspace automatically manages all internal storage
```

# Thread Safety
**Not thread-safe**: Each workspace should be used by only one thread.
For parallel execution, create separate workspace instances per thread.

# See Also
- `build_apslf_pq_workspace`: Constructor function
- `apslf_pqq_solve!`: Main solving function using workspace
- `apslf_pq`: Allocating version for single-use solves
"""
mutable struct APSLFPQWorkspace{TY<:AbstractMatrix{ComplexF64},TF<:LinearAlgebra.Factorization}
   nbus::Int
   slack::Int
   order::Int

   pq::Vector{Int}
   npq::Int
   Yred::TY
   F::TF
   Yslack::Vector{ComplexF64}   # Y[pq, slack], needed for the no-load germ
   germ::Symbol                 # embedding the factorization F was built for
   drow::Vector{ComplexF64}     # (Y·1)[pq] for the :deviation embedding, empty otherwise

   Vcoeff::Matrix{ComplexF64}
   Wcoeff::Matrix{ComplexF64}
   rhs::Vector{ComplexF64}
   Sstar::Vector{ComplexF64}

   Vfull::Vector{ComplexF64}
end

"""
    build_apslf_pq_workspace(Y; slack=1, order=24, use_sparse=:auto, sparse_nbus_min=110)
        -> APSLFPQWorkspace

Construct a pre-allocated workspace for efficient repeated APSLF PQ-only power flow solutions.

# Arguments
- `Y::AbstractMatrix{ComplexF64}`: Nodal admittance matrix (nbus × nbus)

# Keyword Arguments
- `slack::Int = 1`: Slack bus index (fixed voltage reference, excluded from unknowns)
- `order::Int = 24`: APSLF series expansion order for voltage coefficients
- `use_sparse::Union{Bool,Symbol} = :auto`: Sparsity policy for reduced admittance matrix
  - `:auto` → Use sparse format when `nbus >= sparse_nbus_min`
  - `true` → Always convert to sparse format
  - `false` → Keep original matrix format (dense or sparse)
- `sparse_nbus_min::Int = 110`: Threshold for automatic sparse conversion

# Returns
- `APSLFPQWorkspace{TY,TF}`: Pre-allocated workspace with concrete types for efficiency

# Purpose
Eliminates expensive repeated operations when solving multiple APSLF PQ power flows:
- **LU factorization**: Most costly operation (`lu(Yred)`) computed once and reused
- **Memory allocation**: All coefficient matrices and work arrays pre-allocated
- **Type stability**: Concrete parameterized types enable compiler optimizations

# Workspace Components

## System Structure
- Extracts non-slack buses and builds reduced admittance matrix `Yred = Y[pq,pq]`
- Computes and stores LU factorization `F = lu(Yred)` for efficient linear solves
- Sets up indexing maps between full and reduced bus orderings

## Coefficient Storage
- `Vcoeff`: Voltage series coefficients (npq × order+1) in reduced ordering
- `Wcoeff`: Inverse voltage series coefficients (npq × order+1) for recursion
- Column `n+1` stores coefficient of order `n` (0-based math, 1-based Julia)

## Work Arrays
- `rhs`: Right-hand side vector for linear systems (length npq)
- `Sstar`: Conjugated power specifications (length npq)
- `Vfull`: Complete voltage solution vector (length nbus) in full PF ordering

# APSLF Recursion Support
The workspace enables the standard PQ recursion:
```
Y_red * V^(n) = conj(S_red) .* conj(W^(n-1))    for n = 1, 2, ..., order
```

where the inverse series follows: `W^(n) = -(1/V^(0)) * Σ_{m=1}^n V^(m) W^(n-m)`

# Sparsity Policy
- **Automatic detection**: Uses system size to determine optimal matrix format
- **Performance tuning**: `sparse_nbus_min` threshold balances sparse overhead vs benefits
- **Memory efficiency**: Sparse format reduces storage and improves solve times for large systems

# Performance Benefits
- Can be significantly faster for parameter sweeps than repeatedly calling an allocating solver,
  because the LU factorization and work buffers are reused.
- Actual speedup depends on system size, sparsity, factorization cost, and the number of repeated solves.

# Constraints and Notes
- **Fixed structure**: System topology (Y sparsity pattern) and slack bus cannot change
- **PQ buses only**: All non-slack buses treated as PQ with specified complex power
- **Thread safety**: Each thread requires its own workspace instance
- **Order consistency**: All solves with workspace must use same series order

# See Also
- `apslf_pq_solve!`: Main solving function using this workspace
- `apslf_pq`: Allocating version for single-use solutions
- `APSLFPQWorkspace`: Workspace type documentation
"""
function build_apslf_pq_workspace(
   Y::AbstractMatrix{ComplexF64};
   slack::Int = 1,
   order::Int = 24,
   use_sparse::Union{Bool,Symbol} = :auto,
   sparse_nbus_min::Int = 110,
   germ::Symbol = :deviation,
)
   nbus = size(Y, 1)
   @assert size(Y, 2) == nbus
   @assert 1 <= slack <= nbus
   @assert order >= 1
   germ in (:flat, :noload, :deviation) || throw(ArgumentError("germ must be :noload, :deviation or :flat, got :$(germ)"))

   Yuse = maybe_sparse_Y(Y; nbus = nbus, use_sparse = use_sparse, sparse_nbus_min = sparse_nbus_min)

   pq = [i for i = 1:nbus if i != slack]
   npq = length(pq)

   Yred = issparse(Yuse) ? sparse(Yuse[pq, pq]) : Matrix(Yuse[pq, pq])
   drow = germ == :deviation ? apslf_row_sums(Yuse)[pq] : ComplexF64[]
   germ == :deviation && (Yred = Yred - Diagonal(drow))
   F = lu(Yred)

   Vcoeff = zeros(ComplexF64, npq, order + 1)
   Wcoeff = zeros(ComplexF64, npq, order + 1)
   rhs = zeros(ComplexF64, npq)
   Sstar = zeros(ComplexF64, npq)
   Vfull = zeros(ComplexF64, nbus)

   Yslack = Vector{ComplexF64}(Yuse[pq, slack])

   return APSLFPQWorkspace(nbus, slack, order, pq, npq, Yred, F, Yslack, germ, drow, Vcoeff, Wcoeff, rhs, Sstar, Vfull)
end

"""
    apslf_pq_solve!(ws, Y, S; Vslack=1+0im, use_pade=true)
        -> (V, Vcoeff, Wcoeff)

Same semantics as apslf_pq, but:
- reuses ws.F = lu(Yred) across calls (critical for PV secant loop)
- computes rhs in-place (no rhs allocation per order)
- stores coefficients into ws.Vcoeff/ws.Wcoeff buffers
- returns V as a reference to ws.Vfull (copy if you need persistence)
"""
function apslf_pq_solve!(
   ws::APSLFPQWorkspace,
   S::Vector{ComplexF64};
   Vslack::ComplexF64 = 1.0 + 0.0im,
   use_pade::Bool = true,
   germ::Union{Nothing,Symbol} = nothing,
   evaluation_options::Union{Nothing,APSLFEvaluationOptions} = nothing,
   timeout_check = nothing,
   kwargs...,
)
   _handle_deprecated_apslf_germ_kwargs(kwargs; context = "apslf_pq_solve!")
   # The factorization in the workspace belongs to one embedding; :flat and :noload share it.
   germ = germ === nothing ? ws.germ : germ
   if (germ == :deviation) != (ws.germ == :deviation)
      throw(ArgumentError("workspace was built with germ=:$(ws.germ); rebuild it for germ=:$(germ)"))
   end

   @assert length(S) == ws.nbus
   @assert ws.order >= 1

   pq = ws.pq
   npq = ws.npq
   order = ws.order

   Vcoeff = ws.Vcoeff
   Wcoeff = ws.Wcoeff
   rhs = ws.rhs
   Sstar = ws.Sstar

   # ---- APSLF germ V^(0): canonical flat 1∠0, or the no-load solution (germ=:noload)
   if germ == :flat
      @inbounds for k = 1:npq
         Vcoeff[k, 1] = 1.0 + 0.0im
      end
   elseif germ == :noload
      Vcoeff[:, 1] .= apslf_noload_germ(ws.F, ws.Yslack, Vslack)
   elseif germ == :deviation
      @inbounds for k = 1:npq
         Vcoeff[k, 1] = Vslack
      end
   else
      throw(ArgumentError("germ must be :noload, :deviation or :flat, got :$(germ)"))
   end

   # W^(0) = 1/V^(0)
   @inbounds for k = 1:npq
      Wcoeff[k, 1] = 1.0 / Vcoeff[k, 1]
   end

   # Sstar := conj(S[pq])
   @inbounds for k = 1:npq
      Sstar[k] = conj(S[pq[k]])
   end

   # Recursion
   for n = 1:order
      _maybe_check_timeout(timeout_check, :apslf_pq_order)
      # rhs := Sstar .* conj(W^(n-1))  (Wcoeff[:, n] is W^(n-1); reflection condition)
      @inbounds for k = 1:npq
         rhs[k] = Sstar[k] * conj(Wcoeff[k, n])
      end
      if germ == :deviation
         @inbounds for k = 1:npq
            rhs[k] -= ws.drow[k] * Vcoeff[k, n]   # - d .* V^(n-1)
         end
      end

      # Solve Yred * V^(n) = rhs  (store into column n+1)
      ldiv!(@view(Vcoeff[:, n+1]), ws.F, rhs)

      # Inverse series recursion:
      # w_n = -(1/v0) * Σ_{m=1..n} v_m w_{n-m}
      @inbounds for i = 1:npq
         acc = 0.0 + 0.0im
         for m = 1:n
            acc += Vcoeff[i, m+1] * Wcoeff[i, (n-m)+1]
         end
         Wcoeff[i, n+1] = -acc / Vcoeff[i, 1]
      end
   end

   # Assemble full V (PF ordering)
   V = ws.Vfull
   fill!(V, 0.0 + 0.0im)
   V[ws.slack] = Vslack

   eval_opts = isnothing(evaluation_options) ? APSLFEvaluationOptions(mode = use_pade ? :pade : :taylor) : evaluation_options
   @inbounds for (idx, bus) in enumerate(pq)
      c = @view Vcoeff[idx, :]
      V[bus] = evaluate_series(c, eval_opts).voltage
   end

   return V, Vcoeff, Wcoeff
end

# =============================================================================
# Rectangular NR polish (analytic Jacobian, PQ/PV)
# State: x = [Vr(non-slack); Vi(non-slack)]
# Residual per non-slack bus i:
#   PQ: [ΔP_i; ΔQ_i]   with ΔP=Pcalc-Pspec, ΔQ=Qcalc-Qspec
#   PV: [ΔP_i; ΔV_i]   with ΔV=|V|-Vm
# =============================================================================


@inline function _bt_sig(bustype::Vector{Symbol})::UInt64
   h = UInt64(0)
   @inbounds for bt in bustype
      b = _bt_norm(bt)
      # mix: small stable hash combiner
      h ⊻= (UInt64(hash(b)) + 0x9e3779b97f4a7c15) + (h << 6) + (h >> 2)
   end
   return h
end

"""
    NRRectCache

Pre-allocated workspace for efficient Newton-Raphson power flow iterations in rectangular coordinates.

# Fields

## System Structure
- `nbus::Int`: Total number of buses in the system
- `slack::Int`: Slack bus index (fixed voltage reference)

## Bus Indexing and Mapping
- `non_slack::Vector{Int}`: Non-slack bus indices in original PF ordering
- `pos_ns::Vector{Int}`: Mapping from bus index to position in `non_slack` (0 if slack)
- `rowP::Vector{Int}`: Mapping from bus index to ΔP equation row in Jacobian (0 if slack)

## Bus Type Masks (Normalized)
- `pv_mask::BitVector`: Mask for PV buses (true for :PV buses, length nbus)
- `pq_mask::BitVector`: Mask for PQ buses (true for :PQ buses, length nbus)

## Pre-allocated Work Arrays
- `F::Vector{Float64}`: Mismatch vector of length 2×(nbus-1)
- `dx::Vector{Float64}`: Newton step vector of length 2×(nbus-1)
- `I::Vector{ComplexF64}`: Current injection work array (length nbus)
- `Sinj::Vector{ComplexF64}`: Complex power injection work array (length nbus)

# Purpose and Usage
This cache eliminates repeated allocations during Newton-Raphson iterations:

1. **One-time setup**: Build cache with `build_nr_rect_cache(bustype, slack)`
2. **Repeated iterations**: Reuse for mismatch computation and Jacobian assembly
3. **Memory efficiency**: All work arrays pre-allocated to exact required sizes

# State Vector and Equation Structure
The Newton-Raphson method uses rectangular voltage coordinates:

**State vector** (length 2×(nbus-1)):
```
x = [Vr₁, Vr₂, ..., Vrₙ₋₁, Vi₁, Vi₂, ..., Viₙ₋₁]
```

**Mismatch vector** (length 2×(nbus-1)):
```
F = [ΔP₁, ΔQ₁/ΔV₁, ΔP₂, ΔQ₂/ΔV₂, ...]
```

where each non-slack bus contributes two equations:
- **PQ buses**: `[ΔP, ΔQ]` (active and reactive power mismatches)
- **PV buses**: `[ΔP, ΔV]` (active power and voltage magnitude mismatches)

# Indexing Conventions

## Bus-to-Position Mapping
- `pos_ns[bus]`: Position of `bus` in `non_slack` vector (0 if slack)
- `non_slack[k]`: Bus index at position `k` in non-slack ordering

## Equation Row Assignment
- `rowP[bus]`: Row index for ΔP equation of `bus` (0 if slack)
- Row `rowP[bus] + 1`: ΔQ (PQ) or ΔV (PV) equation for same bus

## State Variable Columns
For non-slack bus at position `k = pos_ns[bus]`:
- Column `k`: Real voltage component `Vr`
- Column `k + (nbus-1)`: Imaginary voltage component `Vi`

# Bus Type Normalization
The cache normalizes bus type symbols to uppercase:
- `:pq` or `:PQ` → stored as `:PQ` in masks
- `:pv` or `:PV` → stored as `:PV` in masks
- `:slack` or `:Slack` → excluded from non_slack list

# Memory Layout Optimization
- Contiguous F and dx arrays improve cache locality
- BitVector masks provide memory-efficient bus type queries
- Integer indexing arrays avoid hash table lookups in hot paths

# Example Usage
```julia
# Build cache once for given system structure
cache = build_nr_rect_cache(bustype, slack)

# Newton-Raphson iteration loop
for iter = 1:max_iter
    # Compute mismatch (reuses cache.I, cache.Sinj)
    F = mismatch_rectangular!(cache, Y, V, Pspec, Qspec, Vm)

    # Build Jacobian (uses same work arrays)
    J = build_rect_jac_sparse(cache, Y, V, Vm)

    # Newton step (reuses cache.dx)
    dx_norm, F_norm = nr_refine_step_rect!(cache, Y, V, Pspec, Qspec, Vm; ...)

    # Check convergence
    max(dx_norm, F_norm) < tol && break
end
```

# Thread Safety
**Not thread-safe**: Each thread requires its own cache instance since work arrays
are modified in-place during computations.

# See Also
- `build_nr_rect_cache`: Constructor function
- `mismatch_rectangular!`: Mismatch computation using this cache
- `build_rect_jac_sparse`: Sparse Jacobian assembly using this cache
- `nr_refine_step_rect!`: Complete Newton step using this cache
"""
mutable struct NRRectCache
   nbus::Int
   slack::Int

   non_slack::Vector{Int}     # buses excluding slack
   pos_ns::Vector{Int}        # bus -> position in non_slack (0 if slack)
   rowP::Vector{Int}          # bus -> row index (ΔP) in F/J (0 if slack)
   pv_mask::BitVector         # PV buses mask (normalized)
   pq_mask::BitVector         # PQ buses mask (normalized)

   F::Vector{Float64}         # mismatch vector, length 2*(n-1)
   dx::Vector{Float64}        # step vector, length 2*(n-1)

   I::Vector{ComplexF64}      # work: I = Y*V
   Sinj::Vector{ComplexF64}   # work: S = V .* conj(I)
end

"""
    build_nr_rect_cache(bustype, slack) -> NRRectCache

Construct a pre-allocated cache for efficient Newton-Raphson power flow iterations in rectangular coordinates.

# Arguments
- `bustype::Vector{Symbol}`: Bus type specification for each bus
  - `:PQ` or `:pq`: Active and reactive power specified
  - `:PV` or `:pv`: Active power and voltage magnitude specified
  - `:Slack` or `:slack`: Reference bus with fixed voltage (excluded from equations)
- `slack::Int`: Slack bus index (1-based, must correspond to slack bus in `bustype`)

# Returns
- `NRRectCache`: Pre-allocated workspace containing system structure and work arrays

# Cache Construction
The function builds efficient indexing structures and pre-allocates all work arrays:

## System Structure Analysis
- Identifies non-slack buses and creates position mappings
- Normalizes bus type symbols to uppercase (`:pq` → `:PQ`, etc.)
- Builds BitVector masks for fast bus type queries during iterations

## Indexing Maps
- `non_slack`: List of buses excluding slack (used for state variables)
- `pos_ns[bus]`: Position of `bus` in `non_slack` vector (0 if slack)
- `rowP[bus]`: Row index for ΔP equation of `bus` in mismatch vector (0 if slack)

## Equation Structure
For each non-slack bus, assigns two consecutive rows in the mismatch vector:
- **PQ buses**: `[ΔP, ΔQ]` (active and reactive power mismatches)
- **PV buses**: `[ΔP, ΔV]` (active power and voltage magnitude mismatches)

## Pre-allocated Arrays
- Mismatch vector `F` and Newton step `dx`: length 2×(nbus-1)
- Work arrays `I` and `Sinj` for current and power computations: length nbus
- All arrays sized exactly to avoid runtime allocations

# State Vector Structure
The Newton-Raphson state vector uses rectangular voltage coordinates:
```
x = [Vr₁, Vr₂, ..., Vrₙ₋₁, Vi₁, Vi₂, ..., Viₙ₋₁]
```
where subscripts refer to positions in the `non_slack` bus list.

# Bus Type Normalization
Accepts both lowercase and uppercase bus type symbols:
- `:pq`, `:PQ` → normalized to `:PQ`
- `:pv`, `:PV` → normalized to `:PV`
- `:slack`, `:Slack` → identified but excluded from `non_slack`

# Memory Layout Optimization
- Contiguous arrays improve cache locality during matrix operations
- BitVector masks provide O(1) bus type queries without hash lookups
- Integer indexing arrays enable direct array access in hot computation paths


# Performance Benefits
- **Eliminates allocations**: All work arrays pre-sized and reused
- **Fast indexing**: Direct array access via pre-computed integer maps
- **Cache efficiency**: Contiguous memory layout improves performance
- **Type stability**: Concrete array types enable compiler optimizations

# Constraints
- **Fixed topology**: Bus count, slack bus, and type pattern cannot change after construction
- **Thread safety**: Each thread requires its own cache instance
- **Bus numbering**: Assumes buses numbered 1 to nbus with valid slack index

# See Also
- `NRRectCache`: Cache type documentation with field descriptions
- `mismatch_rectangular!`: Mismatch computation using this cache
- `build_rect_jac_sparse`: Sparse Jacobian assembly using this cache
- `nr_refine_step_rect!`: Complete Newton step using this cache
"""
function build_nr_rect_cache(bustype::Vector{Symbol}, slack::Int)
   n = length(bustype)
   @assert 1 <= slack <= n

   non_slack = Int[]
   sizehint!(non_slack, n - 1)
   for i = 1:n
      i == slack && continue
      push!(non_slack, i)
   end

   pos_ns = zeros(Int, n)
   for (k, bus) in enumerate(non_slack)
      pos_ns[bus] = k
   end

   rowP = zeros(Int, n)
   row = 0
   for i = 1:n
      i == slack && continue
      row += 2
      rowP[i] = row - 1  # ΔP row
   end

   pv_mask = falses(n)
   pq_mask = falses(n)
   @inbounds for i = 1:n
      bt = _bt_norm(bustype[i])
      pv_mask[i] = (bt == :PV)
      pq_mask[i] = (bt == :PQ)
   end

   m = 2 * (n - 1)
   return NRRectCache(
      n,
      slack,
      non_slack,
      pos_ns,
      rowP,
      pv_mask,
      pq_mask,
      zeros(Float64, m),
      zeros(Float64, m),
      zeros(ComplexF64, n),
      zeros(ComplexF64, n),
   )
end

# In-place mismatch computation
"""
    mismatch_rectangular!(cache, Y, V, Pspec, Qspec, Vm) -> Vector{Float64}

Compute power flow mismatch vector in rectangular coordinates for Newton-Raphson iterations.

# Arguments
- `cache::NRRectCache`: Pre-allocated cache containing system structure and work arrays
- `Y::AbstractMatrix{ComplexF64}`: Nodal admittance matrix (nbus × nbus)
- `V::Vector{ComplexF64}`: Current voltage phasor estimates (length nbus)
- `Pspec::Vector{Float64}`: Specified active power injections (p.u., length nbus)
- `Qspec::Vector{Float64}`: Specified reactive power injections (p.u., length nbus)
- `Vm::Vector{Float64}`: Specified voltage magnitudes for PV buses (p.u., length nbus)

# Returns
- `Vector{Float64}`: Mismatch vector F of length 2×(nbus-1) stored in cache.F

# Mismatch Formulation
For each non-slack bus i, computes two residual equations:

**PQ buses:**
- `F[2k-1] = Pi_calc - Pspec[i]`  (active power mismatch)
- `F[2k]   = Qi_calc - Qspec[i]`  (reactive power mismatch)

**PV buses:**
- `F[2k-1] = Pi_calc - Pspec[i]`  (active power mismatch)
- `F[2k]   = |Vi| - Vm[i]`        (voltage magnitude mismatch)

where `Pi_calc = real(Si)` and `Qi_calc = imag(Si)` are computed from:
- `I = Y * V` (nodal current injections)
- `S = V .* conj(I)` (complex power injections)

# State Vector Ordering
The state vector corresponds to rectangular voltage components of non-slack buses:
`x = [Vr₁, Vr₂, ..., Vrₙ₋₁, Vi₁, Vi₂, ..., Viₙ₋₁]`

where subscripts refer to non-slack bus positions.

# Notes
- Slack bus is excluded from equations (voltage fixed)
- Bus types must be :PQ or :PV (normalized internally)
- Modifies `cache.F` in-place and uses `cache.I`, `cache.Sinj` as work arrays
- Function assumes cache was built with same bus type pattern and slack bus

# Example
```julia
cache = build_nr_rect_cache(bustype, slack)
F = mismatch_rectangular!(cache, Y, V, Pspec, Qspec, Vm)
# F now contains power flow mismatches for NR step
```
"""
function mismatch_rectangular!(
   cache::NRRectCache,
   Y::AbstractMatrix{ComplexF64},
   V::Vector{ComplexF64},
   Pspec::Vector{Float64},
   Qspec::Vector{Float64},
   Vm::Vector{Float64},
)
   n = cache.nbus
   @assert length(V) == n
   @assert length(Pspec) == n && length(Qspec) == n && length(Vm) == n

   # I := Y*V ; Sinj := V .* conj(I)
   mul!(cache.I, Y, V)
   @inbounds @simd for i = 1:n
      cache.Sinj[i] = V[i] * conj(cache.I[i])
   end

   F = cache.F
   slack = cache.slack
   rowP = cache.rowP
   pv = cache.pv_mask
   pq = cache.pq_mask

   @inbounds for i = 1:n
      i == slack && continue
      rP = rowP[i]       # ΔP row
      r2 = rP + 1        # ΔQ (PQ) or ΔV (PV) row

      Pi = real(cache.Sinj[i])
      Qi = imag(cache.Sinj[i])

      F[rP] = Pi - Pspec[i]

      if pq[i]
         F[r2] = Qi - Qspec[i]
      elseif pv[i]
         F[r2] = abs(V[i]) - Vm[i]
      else
         error("mismatch_rectangular!: unsupported bus type at bus $i (neither PQ nor PV).")
      end
   end

   return F
end

"""
    build_rect_jac_sparse(cache, Y, V, Vm) -> SparseMatrixCSC{Float64}

Build sparse analytic Jacobian matrix for rectangular Newton-Raphson power flow iterations.

# Arguments
- `cache::NRRectCache`: Pre-allocated cache containing system structure and work arrays
- `Y::SparseMatrixCSC{ComplexF64}`: Sparse nodal admittance matrix (nbus × nbus)
- `V::Vector{ComplexF64}`: Current voltage phasor estimates (length nbus)
- `Vm::Vector{Float64}`: Specified voltage magnitudes for PV buses (p.u., length nbus). Not used by the Jacobian itself (the derivative of |V_i| - Vm_i does not depend on Vm); kept for a uniform call signature with the residual functions.

# Returns
- `SparseMatrixCSC{Float64}`: Sparse Jacobian matrix J of size [2×(nbus-1)] × [2×(nbus-1)]

# Jacobian Structure
The Jacobian relates mismatch vector F to rectangular voltage state vector x:
- **State vector:** `x = [Vr₁, Vr₂, ..., Vrₙ₋₁, Vi₁, Vi₂, ..., Viₙ₋₁]`
- **Mismatch vector:** `F = [ΔP₁, ΔQ₁/ΔV₁, ΔP₂, ΔQ₂/ΔV₂, ...]`

# Derivative Formulations
For complex power injection `Si = Vi * conj(Ii)` where `I = Y*V`:

**Power derivatives (PQ buses):**
- `∂Si/∂Vrⱼ = conj(Ii)δᵢⱼ + Vi*conj(Yᵢⱼ)`
- `∂Si/∂Viⱼ = j*(conj(Ii)δᵢⱼ - Vi*conj(Yᵢⱼ))`

**Voltage magnitude derivatives (PV buses):**
- `∂|Vi|/∂Vri = Vri/|Vi|`  (local derivative only)
- `∂|Vi|/∂Vii = Vii/|Vi|`  (local derivative only)

where δᵢⱼ is the Kronecker delta.

# Sparsity Pattern
Exploits the sparsity pattern of Y matrix:
- Only computes derivatives where Yᵢⱼ ≠ 0
- PV voltage magnitude constraints have only local (diagonal-block) derivatives
- Results in sparse Jacobian with structure inherited from Y

# Bus Type Handling
- **PQ buses:** Both ΔP and ΔQ equations contribute to Jacobian
- **PV buses:** ΔP equation uses power derivatives, ΔV equation uses magnitude derivatives
- **Slack bus:** Excluded from both equations and state variables

# Notes
- Requires Y to be `SparseMatrixCSC{ComplexF64}` for efficient sparse traversal
- Uses `cache.I` as work array for current injections I = Y*V
- Modifies cache.I in-place but preserves input arguments
- For dense matrices, use `build_rect_jac_dense` instead

# Performance
Computational complexity scales with nnz(Y) rather than nbus², making it
efficient for large sparse power systems.

# Example
```julia
cache = build_nr_rect_cache(bustype, slack)
Y_sparse = sparse(Y)  # ensure sparse format
J = build_rect_jac_sparse(cache, Y_sparse, V, Vm)
# J is now ready for Newton step: dx = -(J \\ F)
```
"""
function build_rect_jac_sparse(cache::NRRectCache, Y::SparseMatrixCSC{ComplexF64}, V::Vector{ComplexF64}, ::Vector{Float64})
   n = cache.nbus
   slack = cache.slack
   non_slack = cache.non_slack
   pos_ns = cache.pos_ns
   rowP = cache.rowP
   pv = cache.pv_mask
   pq = cache.pq_mask

   # I = Y*V needed for diagonal term
   mul!(cache.I, Y, V)

   m = 2 * (n - 1)
   nvar = 2 * (n - 1)

   Iidx = Int[]
   Jidx = Int[]
   Vals = Float64[]
   # heuristic sizehint
   sizehint!(Iidx, 16 * nnz(Y))
   sizehint!(Jidx, 16 * nnz(Y))
   sizehint!(Vals, 16 * nnz(Y))

   rv = rowvals(Y)
   nzval = nonzeros(Y)

   # Columns correspond to state variables of non-slack buses only:
   # colVr = pos_ns[j], colVi = (n-1) + pos_ns[j]
   for j = 1:n
      col_pos = pos_ns[j]
      col_pos == 0 && continue  # slack: no state variable

      colVr = col_pos
      colVi = (n - 1) + col_pos

      for ptr in nzrange(Y, j)
         i = rv[ptr]
         i == slack && continue   # slack has no equations

         rP = rowP[i]
         rP == 0 && continue
         r2 = rP + 1

         Yij = nzval[ptr]

         # ∂S_i/∂Vrⱼ = conj(I_i)*δᵢⱼ + V_i*conj(Yᵢⱼ)
         # ∂S_i/∂Viⱼ = j*(conj(I_i)*δᵢⱼ - V_i*conj(Yᵢⱼ))
         dS_dVr = V[i] * conj(Yij)
         dS_dVi = -im * (V[i] * conj(Yij))

         if i == j
            Ii = cache.I[i]
            dS_dVr += conj(Ii)
            dS_dVi += im * conj(Ii)
         end

         dP_Vr = real(dS_dVr)
         dP_Vi = real(dS_dVi)

         # ΔP row always present (PQ & PV)
         if dP_Vr != 0.0
            push!(Iidx, rP)
            push!(Jidx, colVr)
            push!(Vals, dP_Vr)
         end
         if dP_Vi != 0.0
            push!(Iidx, rP)
            push!(Jidx, colVi)
            push!(Vals, dP_Vi)
         end

         # Second row:
         if pq[i]
            dQ_Vr = imag(dS_dVr)
            dQ_Vi = imag(dS_dVi)
            if dQ_Vr != 0.0
               push!(Iidx, r2)
               push!(Jidx, colVr)
               push!(Vals, dQ_Vr)
            end
            if dQ_Vi != 0.0
               push!(Iidx, r2)
               push!(Jidx, colVi)
               push!(Vals, dQ_Vi)
            end
         elseif pv[i]
            # ΔV row handled after loop (local derivative only)
         else
            error("build_rect_jac_sparse: unsupported bus type at bus $i")
         end
      end
   end

   # PV ΔV_i = |V_i| - Vm[i]:
   # ∂|V|/∂Vr = Vr/|V|, ∂|V|/∂Vi = Vi/|V|
   @inbounds for bus in non_slack
      pv[bus] || continue
      rV = rowP[bus] + 1
      pos = pos_ns[bus]
      pos == 0 && continue

      vm = abs(V[bus])
      vm == 0.0 && continue
      dVr = real(V[bus]) / vm
      dVi = imag(V[bus]) / vm

      colVr = pos
      colVi = (n - 1) + pos

      if dVr != 0.0
         push!(Iidx, rV)
         push!(Jidx, colVr)
         push!(Vals, dVr)
      end
      if dVi != 0.0
         push!(Iidx, rV)
         push!(Jidx, colVi)
         push!(Vals, dVi)
      end
   end

   return sparse(Iidx, Jidx, Vals, m, nvar)
end

"""
    build_rect_jac_dense(cache, Y, V, Vm) -> Matrix{Float64}

Build dense analytic Jacobian matrix for rectangular Newton-Raphson power flow iterations.

# Arguments
- `cache::NRRectCache`: Pre-allocated cache containing system structure and work arrays
- `Y::AbstractMatrix{ComplexF64}`: Nodal admittance matrix (nbus × nbus), any format
- `V::Vector{ComplexF64}`: Current voltage phasor estimates (length nbus)
- `Vm::Vector{Float64}`: Specified voltage magnitudes for PV buses (p.u., length nbus). Not used by the Jacobian itself (the derivative of |V_i| - Vm_i does not depend on Vm); kept for a uniform call signature with the residual functions.

# Returns
- `Matrix{Float64}`: Dense Jacobian matrix J of size [2×(nbus-1)] × [2×(nbus-1)]

# Jacobian Structure
The Jacobian relates mismatch vector F to rectangular voltage state vector x:
- **State vector:** `x = [Vr₁, Vr₂, ..., Vrₙ₋₁, Vi₁, Vi₂, ..., Viₙ₋₁]`
- **Mismatch vector:** `F = [ΔP₁, ΔQ₁/ΔV₁, ΔP₂, ΔQ₂/ΔV₂, ...]`

# Derivative Formulations
For complex power injection `Si = Vi * conj(Ii)` where `I = Y*V`:

**Power derivatives (PQ buses):**
- `∂Si/∂Vrⱼ = conj(I_i)*δᵢⱼ + Vi*conj(Yᵢⱼ)`
- `∂Si/∂Viⱼ = j*(conj(I_i)*δᵢⱼ - V_i*conj(Yᵢⱼ))`

**Voltage magnitude derivatives (PV buses):**
- `∂|Vi|/∂Vri = Vri/|Vi|`  (local derivative only)
- `∂|Vi|/∂Vii = Vii/|Vi|`  (local derivative only)

where δᵢⱼ is the Kronecker delta.

# Bus Type Handling
- **PQ buses:** Both ΔP and ΔQ equations contribute to Jacobian
- **PV buses:** ΔP equation uses power derivatives, ΔV equation uses magnitude derivatives
- **Slack bus:** Excluded from both equations and state variables

# Dense vs Sparse
This function computes all possible derivatives regardless of Y sparsity:
- Suitable for small to medium systems where Y is dense or nearly dense
- For large sparse systems, use `build_rect_jac_sparse` for better performance
- Automatically handles any Y matrix format (dense, sparse, structured)

# Notes
- Uses `cache.I` as work array for current injections I = Y*V
- Modifies cache.I in-place but preserves input arguments
- Computational complexity O(nbus²) regardless of Y sparsity
- Returns standard dense Float64 matrix suitable for direct factorization

# Performance
For systems with nbus < ~200 or dense Y matrices, this may be faster than
sparse methods due to reduced overhead. For large sparse systems, prefer
`build_rect_jac_sparse`.

# Example
```julia
cache = build_nr_rect_cache(bustype, slack)
J = build_rect_jac_dense(cache, Y, V, Vm)  # Y can be any matrix type
# J is now ready for Newton step: dx = -(J \\ F)
```
"""
function build_rect_jac_dense(cache::NRRectCache, Y::AbstractMatrix{ComplexF64}, V::Vector{ComplexF64}, ::Vector{Float64})
   n = cache.nbus
   non_slack = cache.non_slack
   pos_ns = cache.pos_ns
   rowP = cache.rowP
   pv = cache.pv_mask
   pq = cache.pq_mask

   # I = Y*V
   mul!(cache.I, Y, V)

   m = 2 * (n - 1)
   nvar = 2 * (n - 1)
   J = zeros(Float64, m, nvar)

   # For each equation bus i (non-slack), for each variable bus j (non-slack)
   @inbounds for j in non_slack
      col_pos = pos_ns[j]
      colVr = col_pos
      colVi = (n - 1) + col_pos

      for i in non_slack
         rP = rowP[i]
         r2 = rP + 1

         Yij = Y[i, j]

         dS_dVr = V[i] * conj(Yij)
         dS_dVi = -im * (V[i] * conj(Yij))

         if i == j
            Ii = cache.I[i]
            dS_dVr += conj(Ii)
            dS_dVi += im * conj(Ii)
         end

         # ΔP row
         J[rP, colVr] += real(dS_dVr)
         J[rP, colVi] += real(dS_dVi)

         # second row
         if pq[i]
            J[r2, colVr] += imag(dS_dVr)
            J[r2, colVi] += imag(dS_dVi)
         elseif pv[i]
            # ΔV later
         else
            error("build_rect_jac_dense: unsupported bus type at bus $i")
         end
      end
   end

   # PV ΔV local derivative
   @inbounds for i in non_slack
      pv[i] || continue
      rV = rowP[i] + 1
      pos = pos_ns[i]
      pos == 0 && continue

      vmi = abs(V[i])
      vmi == 0.0 && continue
      dVr = real(V[i]) / vmi
      dVi = imag(V[i]) / vmi

      J[rV, pos] = dVr
      J[rV, (n-1)+pos] = dVi
   end

   return J
end

"""
    nr_refine_step_rect!(cache, Y, V, Pspec, Qspec, Vm; Vslack, damping=1.0, use_sparse=true)
        -> (maxabs_dx, maxabs_F)

Perform one Newton-Raphson refinement step in rectangular coordinates with analytic Jacobian.

# Arguments
- `cache::NRRectCache`: Pre-allocated cache containing system structure and work arrays
- `Y::AbstractMatrix{ComplexF64}`: Nodal admittance matrix (nbus × nbus)
- `V::Vector{ComplexF64}`: Current voltage phasor estimates (length nbus), modified in-place
- `Pspec::Vector{Float64}`: Specified active power injections (p.u., length nbus)
- `Qspec::Vector{Float64}`: Specified reactive power injections (p.u., length nbus)
- `Vm::Vector{Float64}`: Specified voltage magnitudes for PV buses (p.u., length nbus)

# Keyword Arguments
- `Vslack::ComplexF64`: Fixed slack bus voltage (applied after update)
- `damping::Float64 = 1.0`: Step damping factor ∈ (0,1] for convergence control
- `use_sparse::Bool = true`: Use sparse Jacobian if Y is SparseMatrixCSC, else dense

# Returns
- `Tuple{Float64, Float64}`: Maximum absolute step size and residual magnitude
  - `maxabs_dx`: `maximum(abs.(dx))` - largest voltage component change
  - `maxabs_F`: `maximum(abs.(F))` - largest power flow mismatch

# Algorithm Steps
1. **Mismatch computation**: Compute residual vector F using `mismatch_rectangular!`
2. **Jacobian assembly**: Build analytic Jacobian J (sparse or dense based on Y type)
3. **Linear solve**: Compute Newton step `dx = -(J \\ F)`
4. **Voltage update**: Apply damped step to non-slack bus voltages in rectangular form
5. **Slack enforcement**: Restore slack bus voltage to specified value

# State Vector Structure
Updates voltage components of non-slack buses in rectangular coordinates:
- `dx = [dVr₁, dVr₂, ..., dVrₙ₋₁, dVi₁, dVi₂, ..., dViₙ₋₁]`
- Voltage update: `V[bus] += damping * complex(dVr, dVi)` for non-slack buses

# Convergence Monitoring
The returned values are typically used for convergence assessment:
- `maxabs_dx < tolerance`: Step size convergence criterion
- `maxabs_F < tolerance`: Residual convergence criterion

# Bus Type Handling
- **PQ buses**: Both active and reactive power equations enforced
- **PV buses**: Active power + voltage magnitude constraint equations
- **Slack bus**: Voltage held fixed (excluded from state vector and equations)

# Performance Notes
- Jacobian sparsity automatically detected from Y matrix type
- For `SparseMatrixCSC{ComplexF64}` and `use_sparse=true`: uses efficient sparse algorithms
- For dense matrices or `use_sparse=false`: uses dense linear algebra
- Cache reuse eliminates allocations for mismatch vector and work arrays

# Example
```julia
cache = build_nr_rect_cache(bustype, slack)
for iter = 1:max_iter
    dx_norm, F_norm = nr_refine_step_rect!(
        cache, Y, V, Pspec, Qspec, Vm;
        Vslack=1.0+0.0im, damping=0.8
    )
    max(dx_norm, F_norm) < 1e-10 && break  # converged
end
```
"""
function nr_refine_step_rect!(
   cache::NRRectCache,
   Y::AbstractMatrix{ComplexF64},
   V::Vector{ComplexF64},
   Pspec::Vector{Float64},
   Qspec::Vector{Float64},
   Vm::Vector{Float64};
   Vslack::ComplexF64,
   damping::Float64 = 1.0,
   use_sparse::Bool = true,
)
   # Residual
   F = mismatch_rectangular!(cache, Y, V, Pspec, Qspec, Vm)

   # Jacobian
   J = if use_sparse && (Y isa SparseMatrixCSC{ComplexF64})
      build_rect_jac_sparse(cache, Y, V, Vm)
   else
      build_rect_jac_dense(cache, Y, V, Vm)
   end

   # Solve J*dx = -F
   cache.dx .= -(J \ F)

   # Apply update to V (non-slack only)
   non_slack = cache.non_slack
   n1 = cache.nbus - 1
   @inbounds for (k, bus) in enumerate(non_slack)
      dVr = cache.dx[k]
      dVi = cache.dx[n1+k]
      V[bus] += damping * complex(dVr, dVi)
   end

   # Keep slack fixed
   V[cache.slack] = Vslack

   return (maximum(abs.(cache.dx)), maximum(abs.(F)))
end



# =============================================================================
# src/solver_core.jl  (DROP-IN ADDITION)
# =============================================================================
# Adds a sparse variant of the direct PV APSLF kernel:
#   - Builds the augmented real system Asys as SparseMatrixCSC{Float64}
#   - Assembles the Y-block by iterating nnz(Y) (no O(n^2) dense fill)
#   - Uses sparse LU once, then ldiv! per series order
#
# Notes:
# - Expects Y to be sparse for best performance; will convert dense Y to sparse.
# - Uses the same output contract as apslf_pf_pv_direct.
# - Docstrings/comments are in English (as per your preference).

"""
    apslf_pf_pv_direct_sparse(Y, bustype, Pspec, Qspec, Vm; kwargs...) -> (V, Qpv, Vcoeff, Wcoeff, Qcoeff)

Sparse variant of `apslf_pf_pv_direct`.

Key differences vs dense version:
- Builds the constant augmented real system `Asys` as sparse CSC (Float64).
- Assembles network block from `Y` using nnz iteration (sparse traversal).
- Uses sparse LU (`lu(Asys)`) once and reuses it across all series orders.

This is the right kernel for batch processing on large sparse networks and is
also the correct starting point for a future GPU sparse-direct implementation.

See `apslf_pf_pv_direct` docstring for algorithmic background and semantics.
"""
function apslf_pf_pv_direct_sparse(
   Y::AbstractMatrix{ComplexF64},
   bustype::Vector{Symbol},
   Pspec::Vector{Float64},
   Qspec::Vector{Float64},
   Vm::Vector{Float64};
   slack::Int = 1,
   Vslack::ComplexF64 = 1.0 + 0.0im,
   order::Int = 24,
   use_pade::Bool = true,
   germ::Symbol = :deviation,
   evaluation_options::Union{Nothing,APSLFEvaluationOptions} = nothing,
   debug::Bool = false,
   debug_every::Int = 1,
   debug_maxn::Int = order,
   self_check::Bool = true,
   mis_tol_p::Float64 = 1e-6,
   mis_tol_q::Float64 = 1e-6,
   v_tol::Float64 = 1e-6,
   qpv_tol::Float64 = 1e-6,
   qpv_from_inj::Bool = true,
   timeout_check = nothing,
)
   nbus = length(bustype)
   @assert size(Y, 1) == nbus && size(Y, 2) == nbus
   @assert length(Pspec) == nbus && length(Qspec) == nbus && length(Vm) == nbus
   @assert 1 <= slack <= nbus
   @assert order >= 1

   # Force sparse Y for nnz traversal (cheap if already sparse)
   Yuse = issparse(Y) ? (Y::SparseMatrixCSC{ComplexF64,Int}) : sparse(Y)

   # Identify PV buses (accept :pv or :PV etc.)
   pv = findall(i -> _is_pv(bustype[i]), 1:nbus)
   npv = length(pv)

   # Non-slack buses and position maps
   nonslack = [i for i = 1:nbus if i != slack]
   nn = length(nonslack)

   pos_ns = zeros(Int, nbus)      # bus -> position in nonslack (0 if slack)
   @inbounds for (t, bus) in enumerate(nonslack)
      pos_ns[bus] = t
   end

   pv_pos = zeros(Int, nbus)      # bus -> position in pv list (0 if not pv)
   @inbounds for (k, bus) in enumerate(pv)
      pv_pos[bus] = k
   end

   # -----------------------
   # Coefficient arrays (PF ordering)
   # -----------------------
   Vcoeff = zeros(ComplexF64, nbus, order + 1)   # col 1 => V^(0)
   Wcoeff = zeros(ComplexF64, nbus, order + 1)   # col 1 => W^(0)
   Qcoeff = zeros(Float64, npv, order)           # col n => Q^(n-1)  (SHIFTED)

   # Germ: canonical flat 1∠0 (germ=:flat) or no-load solution (germ=:noload)
   V0 = fill(1.0 + 0.0im, nbus)
   V0[nonslack] .= apslf_germ(Yuse, nonslack, slack, Vslack, germ)
   V0[slack] = Vslack
   # :deviation embedding: constant matrix Y0 = Y - diag(d), d = Y·1 ramped with s on the RHS
   drow = germ == :deviation ? apslf_row_sums(Yuse) : zeros(ComplexF64, nbus)
   Vcoeff[:, 1] .= V0
   Wcoeff[:, 1] .= 1.0 ./ V0
   Wcoeff[slack, 1] = 1.0 / Vslack
   npv > 0 && (Qcoeff[:, 1] .= 0.0)

   # -----------------------
   # Build sparse augmented real system Asys (constant for all n>0)
   # Unknowns: x = [Vr(nonslack); Vi(nonslack); Q^(n-1)(pv)]
   # -----------------------
   nunk = 2 * nn + npv

   Iidx = Int[]
   Jidx = Int[]
   Vals = Float64[]
   # Rough size hint: ~4*nnz(Y) for Y block + small PV additions
   sizehint!(Iidx, 4 * nnz(Yuse) + 8 * npv + 64)
   sizehint!(Jidx, 4 * nnz(Yuse) + 8 * npv + 64)
   sizehint!(Vals, 4 * nnz(Yuse) + 8 * npv + 64)

   # conj(W^(0)) for PV coupling term
   W0c = conj.(Wcoeff[:, 1])

   # Assemble Y-block by iterating over nnz(Y):
   # For each nonzero Y[i,j] where both i and j are non-slack,
   # add contributions to real/imag split equations.
   rv = rowvals(Yuse)
   nzval = nonzeros(Yuse)

   @inbounds for j = 1:nbus
      tj = pos_ns[j]
      tj == 0 && continue  # slack column => no state variable
      colVr = tj
      colVi = nn + tj

      for ptr in nzrange(Yuse, j)
         i = rv[ptr]
         ti = pos_ns[i]
         ti == 0 && continue  # slack row => no equations

         yij = nzval[ptr]
         g = real(yij)
         b = imag(yij)

         r_re = 2 * (ti - 1) + 1
         r_im = r_re + 1

         # Real row:
         #   Asys[r_re, Vr_j] += g
         #   Asys[r_re, Vi_j] += -b
         push!(Iidx, r_re)
         push!(Jidx, colVr)
         push!(Vals, g)
         push!(Iidx, r_re)
         push!(Jidx, colVi)
         push!(Vals, -b)

         # Imag row:
         #   Asys[r_im, Vr_j] += b
         #   Asys[r_im, Vi_j] += g
         push!(Iidx, r_im)
         push!(Jidx, colVr)
         push!(Vals, b)
         push!(Iidx, r_im)
         push!(Jidx, colVi)
         push!(Vals, g)
      end
   end

   # :deviation embedding: subtract the row sums d from the diagonal (Y0 = Y - diag(d))
   if germ == :deviation
      @inbounds for (ti, i) in enumerate(nonslack)
         g = real(drow[i])
         b = imag(drow[i])
         r_re = 2 * (ti - 1) + 1
         r_im = r_re + 1
         push!(Iidx, r_re); push!(Jidx, ti); push!(Vals, -g)
         push!(Iidx, r_re); push!(Jidx, nn + ti); push!(Vals, b)
         push!(Iidx, r_im); push!(Jidx, ti); push!(Vals, -b)
         push!(Iidx, r_im); push!(Jidx, nn + ti); push!(Vals, -g)
      end
   end

   # PV coupling term on LHS for PV equation buses i:
   #   (-j) Q^(n-1)_i * conj(W^(0)_i) moved to LHS
   #   contributes to the network equation rows at bus i, column c_q
   @inbounds for i in pv
      ti = pos_ns[i]
      ti == 0 && continue
      kpv = pv_pos[i]
      @assert kpv > 0
      c_q = 2 * nn + kpv

      r_re = 2 * (ti - 1) + 1
      r_im = r_re + 1

      # RHS term (-j) Q^(n-1) W0* moved to the LHS: (+j)Q*W0* = -Q*imag(W0*) + j*(Q*real(W0*))
      # (sign consistent with the known lower-order Q terms on the RHS)
      push!(Iidx, r_re)
      push!(Jidx, c_q)
      push!(Vals, -imag(W0c[i]))
      push!(Iidx, r_im)
      push!(Jidx, c_q)
      push!(Vals, real(W0c[i]))
   end

   # PV magnitude constraints rows:
   #   Re(conj(V0_i)*V_i^(n)) = eps_i^(n)
   # With V_i^(n)=Vr + jVi:
   #   Re(conj(V0)*(Vr+jVi)) = Re(V0)*Vr + Im(V0)*Vi
   @inbounds for (k, bus) in enumerate(pv)
      ti = pos_ns[bus]
      ti == 0 && continue
      r = 2 * nn + k
      push!(Iidx, r)
      push!(Jidx, ti)
      push!(Vals, real(V0[bus]))
      push!(Iidx, r)
      push!(Jidx, nn + ti)
      push!(Vals, imag(V0[bus]))
   end

   Asys = sparse(Iidx, Jidx, Vals, nunk, nunk)
   Fsys = lu(Asys)  # sparse LU

   # -----------------------
   # Debug helpers
   # -----------------------
   function dbg_print_header()
      println("\n[apslf_pf_pv_direct_sparse] DEBUG")
      println(
         "  nbus=$(nbus), slack=$(slack), npv=$(npv), order=$(order), use_pade=$(use_pade), germ=$(germ)",
      )
      println("  Asys size = $(size(Asys)), nnz(Asys)=$(nnz(Asys))")
      println("  columns: n | ||rhs||₂ | max|Vn| | ||Vn||₂ | max|Q_{n-1}| | ||Q_{n-1}||₂")
   end

   function dbg_print_row(n::Int, rhs::Vector{Float64})
      Vn = Vcoeff[:, n+1]
      maxV = maximum(abs.(Vn))
      nrmV = norm(Vn)

      maxQ = 0.0
      nrmQ = 0.0
      if npv > 0
         Qn1 = Qcoeff[:, n]  # Q_{n-1} stored at column n
         maxQ = maximum(abs.(Qn1))
         nrmQ = norm(Qn1)
      end

      @printf("  %3d | %10.3e | %10.3e | %10.3e | %10.3e | %10.3e\n", n, norm(rhs), maxV, nrmV, maxQ, nrmQ)
   end

   debug && dbg_print_header()

   # -----------------------
   # Convolution helpers (same semantics as dense kernel)
   # -----------------------
   function conv_Q_Wc_known_shifted(n::Int, kpv::Int, bus::Int)
      acc = 0.0 + 0.0im
      n <= 1 && return acc
      @inbounds for m = 0:(n-2)
         Qm = Qcoeff[kpv, m+1]                      # Q^(m)
         Wc = conj(Wcoeff[bus, (n-1-m)+1])          # conj(W^(n-1-m))
         acc += Qm * Wc
      end
      return acc
   end

   function conv_V_Vc_mid(n::Int, bus::Int)
      acc = 0.0 + 0.0im
      n <= 1 && return acc
      @inbounds for m = 1:(n-1)
         acc += Vcoeff[bus, m+1] * conj(Vcoeff[bus, (n-m)+1])
      end
      return acc
   end

   # -----------------------
   # Main recursion
   # -----------------------
   S = ComplexF64.(Pspec, Qspec)
   Sstar = conj.(S)

   rhs = zeros(Float64, nunk)
   x = zeros(Float64, nunk)
   Wc_prev = zeros(ComplexF64, nbus)  # conj(W^(n-1)) buffer

   for n = 1:order
      _maybe_check_timeout(timeout_check, :apslf_direct_sparse_order)
      @inbounds for i = 1:nbus
         Wc_prev[i] = conj(Wcoeff[i, n])  # Wcoeff col n is W^(n-1)
      end

      fill!(rhs, 0.0)

      # Network equation RHS rows (per non-slack bus i)
      @inbounds for i in nonslack
         ti = pos_ns[i]
         r_re = 2 * (ti - 1) + 1
         r_im = r_re + 1

         if _is_pq(bustype[i])
            f = Sstar[i] * Wc_prev[i]
            germ == :deviation && (f -= drow[i] * Vcoeff[i, n])   # - d_i V_i^(n-1)
            rhs[r_re] = real(f)
            rhs[r_im] = imag(f)

         elseif _is_pv(bustype[i])
            kpv = pv_pos[i]
            @assert kpv > 0
            acc = conv_Q_Wc_known_shifted(n, kpv, i)
            f = Pspec[i] * Wc_prev[i] - 1.0im * acc
            germ == :deviation && (f -= drow[i] * Vcoeff[i, n])   # - d_i V_i^(n-1)
            rhs[r_re] = real(f)
            rhs[r_im] = imag(f)

         else
            error("apslf_pf_pv_direct_sparse: unsupported bustype at bus $i: $(bustype[i])")
         end
      end

      # PV magnitude constraint RHS eps_i^(n)
      @inbounds for (k, bus) in enumerate(pv)
         r = 2 * nn + k
         δn1 = (n == 1) ? 1.0 : 0.0
         eps = 0.5 * δn1 * (Vm[bus]^2 - abs(V0[bus])^2) - 0.5 * real(conv_V_Vc_mid(n, bus))
         rhs[r] = eps
      end

      # Solve sparse LU once per order
      ldiv!(x, Fsys, rhs)

      # Write back V^(n) (complex) for non-slack
      @inbounds for (t, bus) in enumerate(nonslack)
         Vcoeff[bus, n+1] = complex(x[t], x[nn+t])
      end
      Vcoeff[slack, n+1] = 0.0 + 0.0im

      # Write back Q^(n-1) (SHIFTED) into column n
      @inbounds for k = 1:npv
         Qcoeff[k, n] = x[2*nn+k]
      end

      # Update inverse series W = 1/V
      @inbounds for i = 1:nbus
         if i == slack
            Wcoeff[i, n+1] = 0.0 + 0.0im
            continue
         end
         acc = 0.0 + 0.0im
         for m = 1:n
            acc += Vcoeff[i, m+1] * Wcoeff[i, (n-m)+1]
         end
         Wcoeff[i, n+1] = -acc / Vcoeff[i, 1]
      end

      if debug && (n <= debug_maxn) && (n % debug_every == 0)
         dbg_print_row(n, rhs)
      end
   end

   # -----------------------
   # Evaluate at s=1 via Padé or series sum
   # -----------------------
   V = zeros(ComplexF64, nbus)
   V[slack] = Vslack
   Qpv = zeros(Float64, npv)

   eval_opts = isnothing(evaluation_options) ? APSLFEvaluationOptions(mode = use_pade ? :pade : :taylor) : evaluation_options
   @inbounds for i = 1:nbus
      i == slack && continue
      cV = @view Vcoeff[i, :]
      V[i] = evaluate_series(cV, eval_opts; collect_logs = false).voltage
   end
   if npv > 0
      @inbounds for k = 1:npv
         cQ = ComplexF64.(Qcoeff[k, :], 0.0)
         Qpv[k] = real(evaluate_series(cQ, eval_opts; collect_logs = false).voltage)
      end
   end

   # Optional: recompute PV Q from injections for consistency
   if qpv_from_inj && npv > 0
      Sinj = calc_injections(Yuse, V)
      @inbounds for (k, bus) in enumerate(pv)
         Qpv[k] = imag(Sinj[bus])
      end
   end

   # Self-check (prints only on violation)
   if self_check
      Sinj = calc_injections(Yuse, V)

      maxP = 0.0
      @inbounds for i = 1:nbus
         i == slack && continue
         maxP = max(maxP, abs(real(Sinj[i]) - Pspec[i]))
      end

      maxQpq = 0.0
      @inbounds for i = 1:nbus
         i == slack && continue
         _is_pq(bustype[i]) || continue
         maxQpq = max(maxQpq, abs(imag(Sinj[i]) - Qspec[i]))
      end

      maxVpv = 0.0
      @inbounds for bus in pv
         maxVpv = max(maxVpv, abs(abs(V[bus]) - Vm[bus]))
      end

      maxQpvDiff = 0.0
      if npv > 0
         @inbounds for (k, bus) in enumerate(pv)
            maxQpvDiff = max(maxQpvDiff, abs(Qpv[k] - imag(Sinj[bus])))
         end
      end

      if (maxP > mis_tol_p) || (maxQpq > mis_tol_q) || (maxVpv > v_tol) || (maxQpvDiff > qpv_tol)
         println("\n[apslf_pf_pv_direct_sparse] SELF-CHECK (only printed on violation)")
         @printf("  max |P mismatch| (p.u.)         = %.6e   (tol %.1e)\n", maxP, mis_tol_p)
         @printf("  max |Q mismatch| on PQ (p.u.)   = %.6e   (tol %.1e)\n", maxQpq, mis_tol_q)
         @printf("  max | |V|-Vsp | on PV (p.u.)    = %.6e   (tol %.1e)\n", maxVpv, v_tol)
         @printf("  max |Qpv_series - imag(Sinj)|   = %.6e   (tol %.1e)\n", maxQpvDiff, qpv_tol)
      end
   end

   return V, Qpv, Vcoeff, Wcoeff, Qcoeff
end


"""
    apslf_pf_pv_direct(Y, bustype, Pspec, Qspec, Vm; kwargs...) -> (V, Qpv, Vcoeff, Wcoeff, Qcoeff)

Direct APSLF PV (GridCal-style augmented real system), forced flat germ.

Purpose:
Solve AC power flow with PV buses directly in the APSLF recursion by solving,
at each series order `n`, a single augmented real linear system whose unknowns
include the PV reactive power coefficient `Q^(n-1)` (shifted storage).

Key series:
- Voltages: `V(s) = Σ_{n=0..order} V^(n) s^n`
- Inverse:  `W(s) = 1/V(s) = Σ_{n=0..order} W^(n) s^n`
- PV Q:     `Q_pv(s) = Σ_{n=0..order-1} Q^(n) s^n`

SHIFTED Q storage:
`Qcoeff :: npv × order` where column `n` stores `Q^(n-1)`.

Augmented real system (constant for all `n>0`):
Unknown vector:
`x = [ Vr(nonslack); Vi(nonslack); Q^(n-1)(pv) ]`
Size: `2*(nbus-1) + npv`.

Evaluation at `s=1`:
- Via Padé or direct series sum for `V` and `Qpv`.

Optional:
- NR polishing (finite-difference Jacobian) for robustness/benchmark parity.
- Consistency override `qpv_from_inj=true` recomputes PV `Q` from final injections.
"""
function apslf_pf_pv_direct(
   Y::AbstractMatrix{ComplexF64},
   bustype::Vector{Symbol},
   Pspec::Vector{Float64},
   Qspec::Vector{Float64},
   Vm::Vector{Float64};
   slack::Int = 1,
   Vslack::ComplexF64 = 1.0 + 0.0im,
   order::Int = 24,
   use_pade::Bool = true,
   germ::Symbol = :deviation,
   evaluation_options::Union{Nothing,APSLFEvaluationOptions} = nothing,
   debug::Bool = false,
   debug_every::Int = 1,
   debug_maxn::Int = order,
   self_check::Bool = true,
   mis_tol_p::Float64 = 1e-6,
   mis_tol_q::Float64 = 1e-6,
   v_tol::Float64 = 1e-6,
   qpv_tol::Float64 = 1e-6,
   qpv_from_inj::Bool = true,
   timeout_check = nothing,
)
   nbus = length(bustype)
   @assert size(Y, 1) == nbus && size(Y, 2) == nbus
   @assert length(Pspec) == nbus && length(Qspec) == nbus && length(Vm) == nbus
   @assert 1 <= slack <= nbus
   @assert order >= 1
   Yuse = Y
   pv = findall(i -> _is_pv(bustype[i]), 1:nbus)
   npv = length(pv)
   # Fast maps (avoid Dict/findfirst in hot paths)
   nonslack = [i for i = 1:nbus if i != slack]
   nn = length(nonslack)

   pos_ns = zeros(Int, nbus)         # bus -> position in nonslack (0 if slack)
   @inbounds for (t, bus) in enumerate(nonslack)
      pos_ns[bus] = t
   end

   pv_pos = zeros(Int, nbus)         # bus -> position in pv list (0 if not pv)
   @inbounds for (k, bus) in enumerate(pv)
      pv_pos[bus] = k
   end

   # -----------------------
   # Allocate coefficient arrays (PF ordering)
   #
   # Convention: column 1 stores order 0 coefficient:
   #   Vcoeff[:,1] = V^(0), Vcoeff[:,n+1] = V^(n)
   #   Wcoeff similarly.
   #   Qcoeff[:,n] stores Q^(n-1) (SHIFTED).
   # -----------------------
   Vcoeff = zeros(ComplexF64, nbus, order + 1)  # col 1 => V⁽0⁾
   Wcoeff = zeros(ComplexF64, nbus, order + 1)  # col 1 => W⁽0⁾
   Qcoeff = zeros(Float64, npv, order)        # col 1 => Q⁽0⁾  (SHIFTED: stores Q⁽n-1⁾ at step n)

   # -----------------------
   # Germ: canonical flat 1∠0 (germ=:flat) or no-load solution (germ=:noload)
   #
   # V^(0) = 1 at all buses (flat) or the no-load voltage profile (noload),
   # slack set to Vslack. W^(0) = 1/V^(0). For PV: Q^(0) = 0.
   # -----------------------
   V0 = fill(1.0 + 0.0im, nbus)
   V0[nonslack] .= apslf_germ(Yuse, nonslack, slack, Vslack, germ)
   V0[slack] = Vslack
   # :deviation embedding: constant matrix Y0 = Y - diag(d), d = Y·1 ramped with s on the RHS
   drow = germ == :deviation ? apslf_row_sums(Yuse) : zeros(ComplexF64, nbus)
   Vcoeff[:, 1] .= V0
   Wcoeff[:, 1] .= 1.0 ./ V0
   Wcoeff[slack, 1] = 1.0 / Vslack

   if npv > 0
      Qcoeff[:, 1] .= 0.0  # Q⁽0⁾ = 0
   end

   # -----------------------
   # Build constant augmented real system Asys for all n>0
   #
   # Unknown vector:
   #   x = [ Vr(nonslack); Vi(nonslack); Q^(n-1)(pv) ]
   #
   # Network equations (for each nonslack bus i) in real/imag split:
   #   (Y V^(n))_i + PV-coupling-term = rhs_i^(n)
   #
   # PV coupling term corresponds to:
   #   (-j) Q^(n-1) * conj(W^(0))   moved to LHS
   #
   # PV magnitude constraint row for each PV bus i:
   #   Re( conj(V0_i) * V_i^(n) ) = eps_i^(n)
   # With V0_i = 1 (flat), this is simply Re(V_i^(n)) = eps_i^(n).
   # (But code uses general expression via V0[bus].)
   # -----------------------
   nunk = 2 * nn + npv

   Asys = zeros(Float64, nunk, nunk)

   # conj(W^(0)) used in PV coupling term
   W0c = conj.(Wcoeff[:, 1])

   # Network equation rows for each nonslack bus i
   for (row_i, i) in enumerate(nonslack)
      r_re = 2 * (row_i - 1) + 1
      r_im = r_re + 1

      # Contribution of Y * V^(n) for nonslack voltage unknowns.
      # If V_k^(n) = Vr + jVi, then:
      #   Y_ik V_k = (g + jb)(Vr + jVi)
      #           = (g*Vr - b*Vi) + j(b*Vr + g*Vi)
      @inbounds for (col_t, kbus) in enumerate(nonslack)
         c_vr = col_t
         c_vi = nn + col_t

         yik = Yuse[i, kbus]
         (germ == :deviation && kbus == i) && (yik -= drow[i])   # Y0 = Y - diag(d)
         g = real(yik)
         b = imag(yik)

         Asys[r_re, c_vr] += g
         Asys[r_re, c_vi] += -b

         Asys[r_im, c_vr] += b
         Asys[r_im, c_vi] += g
      end

      # PV coupling: the RHS term (-j) Q^(n-1) * conj(W^(0)) moved to the LHS
      # becomes (+j) Q W0*, with (+j)QW0* = -Q*imag(W0*) + j*(Q*real(W0*)).
      # The sign must match the known lower-order Q terms on the RHS
      # (-j Σ Q^(m) conj(W^(n-1-m))); with the opposite sign the recursion is
      # inconsistent from order 2 on and the PV active power is not met.
      if _is_pv(bustype[i])
         pvk = pv_pos[i]
         @assert pvk > 0
         c_q = 2 * nn + pvk
         Asys[r_re, c_q] += -imag(W0c[i])
         Asys[r_im, c_q] += real(W0c[i])
      end
   end

   # PV magnitude constraints rows:
   # Re( conj(V0_i) * V_i^(n) ) = eps_i^(n)
   # With V_i^(n) = Vr + jVi, and conj(V0)=Re(V0)-jIm(V0):
   #   Re( conj(V0)*(Vr+jVi) ) = Re(V0)*Vr + Im(V0)*Vi
   for (k, bus) in enumerate(pv)
      r = 2 * nn + k
      t = pos_ns[bus]
      @assert t > 0
      Asys[r, t] = real(V0[bus])     # Vr coefficient
      Asys[r, nn+t] = imag(V0[bus])     # Vi coefficient
   end

   Fsys = lu(Asys)

   # -----------------------
   # Debug helpers (optional)
   # -----------------------
   function dbg_print_header()
      println("\n[apslf_pf_pv_direct] DEBUG")
      println(
         "  nbus=$(nbus), slack=$(slack), npv=$(npv), order=$(order), use_pade=$(use_pade), germ=$(germ)",
      )
      println("  Asys size = $(size(Asys))")
      println("  columns: n | ||rhs||₂ | max|Vn| | ||Vn||₂ | max|Q_{n-1}| | ||Q_{n-1}||₂")
   end

   function dbg_print_row(n::Int, rhs::Vector{Float64})
      Vn = Vcoeff[:, n+1]
      maxV = maximum(abs.(Vn))
      nrmV = norm(Vn)

      maxQ = 0.0
      nrmQ = 0.0
      if npv > 0
         Qn1 = Qcoeff[:, n]  # Q_{n-1} stored at column n
         maxQ = maximum(abs.(Qn1))
         nrmQ = norm(Qn1)
      end

      @printf("  %3d | %10.3e | %10.3e | %10.3e | %10.3e | %10.3e\n", n, norm(rhs), maxV, nrmV, maxQ, nrmQ)
   end

   if debug
      dbg_print_header()
   end

   # -----------------------
   # Convolution helpers (known-only, SHIFTED Q storage)
   #
   # At recursion step n:
   # - Unknown: Q^(n-1) for PV buses (solved via augmented system)
   # - Known:   Q^(0..n-2) already stored in Qcoeff[:,1..n-1]
   #
   # conv_Q_Wc_known_shifted(n):
   #   Σ_{m=0..n-2} Q^(m) * conj(W^(n-1-m))
   #
   # conv_V_Vc_mid(n):
   #   Σ_{m=1..n-1} V^(m) * conj(V^(n-m))
   # -----------------------
   function conv_Q_Wc_known_shifted(n::Int, kpv::Int, bus::Int)
      acc = 0.0 + 0.0im
      if n <= 1
         return acc
      end
      @inbounds for m = 0:(n-2)
         Qm = Qcoeff[kpv, m+1]                      # Q⁽m⁾
         Wc = conj(Wcoeff[bus, (n-1-m)+1])        # W*⁽n-1-m⁾
         acc += Qm * Wc
      end
      return acc
   end

   function conv_V_Vc_mid(n::Int, bus::Int)
      acc = 0.0 + 0.0im
      if n <= 1
         return acc
      end
      @inbounds for m = 1:(n-1)
         acc += Vcoeff[bus, m+1] * conj(Vcoeff[bus, (n-m)+1])
      end
      return acc
   end

   # -----------------------
   # Main recursion
   # -----------------------
   S = ComplexF64.(Pspec, Qspec)
   Sstar = conj.(S)
   rhs = zeros(Float64, nunk)
   x = zeros(Float64, nunk)
   Wc_prev = zeros(ComplexF64, nbus)   # conj(W^(n-1)) for all buses (buffer)
   for n = 1:order
      _maybe_check_timeout(timeout_check, :apslf_direct_order)
      # W*^(n-1) (since Wcoeff col n is W^(n-1))
      @inbounds for i = 1:nbus
         Wc_prev[i] = conj(Wcoeff[i, n])
      end

      fill!(rhs, 0.0)


      # Network equations for each nonslack bus i
      for (row_i, i) in enumerate(nonslack)
         r_re = 2 * (row_i - 1) + 1
         r_im = r_re + 1

         if _is_pq(bustype[i])
            # PQ coefficient equation:
            #   (Y V^(n))_i = conj(S_i) * conj(W_i^(n-1))
            f = Sstar[i] * Wc_prev[i]
            germ == :deviation && (f -= drow[i] * Vcoeff[i, n])   # - d_i V_i^(n-1)
            rhs[r_re] = real(f)
            rhs[r_im] = imag(f)

         elseif _is_pv(bustype[i])
            # PV coefficient equation arranged as:
            #   (Y V^(n))_i + (-j) Q^(n-1)_i * conj(W_i^(0))
            #       = P_i * conj(W_i^(n-1)) - j * Σ_{m=0..n-2} Q^(m)_i * conj(W_i^(n-1-m))
            #
            # Here: RHS is the known part
            kpv = pv_pos[i]
            @assert kpv > 0
            acc = conv_Q_Wc_known_shifted(n, kpv, i)

            f = Pspec[i] * Wc_prev[i] - 1.0im * acc
            germ == :deviation && (f -= drow[i] * Vcoeff[i, n])   # - d_i V_i^(n-1)
            rhs[r_re] = real(f)
            rhs[r_im] = imag(f)

         else
            error("Unsupported bustype at bus $i: $(bustype[i])")
         end
      end

      # PV magnitude constraint RHS ε_i^(n):
      #   Re(conj(V0_i) V_i^(n)) = 0.5 δ_{n1}(Vm_i^2 - |V0_i|^2) - 0.5 Re( Σ_{m=1..n-1} V_i^(m) conj(V_i^(n-m)) )
      for (k, bus) in enumerate(pv)
         r = 2 * nn + k
         δn1 = (n == 1) ? 1.0 : 0.0
         eps = 0.5 * δn1 * (Vm[bus]^2 - abs(V0[bus])^2) - 0.5 * real(conv_V_Vc_mid(n, bus))
         rhs[r] = eps
      end

      # Solve for unknowns at order n
      #x = Fsys \ rhs
      ldiv!(x, Fsys, rhs)

      # Write back V^(n) for nonslack buses (complex form)
      @inbounds for (t, bus) in enumerate(nonslack)
         Vr = x[t]
         Vi = x[nn+t]
         Vcoeff[bus, n+1] = complex(Vr, Vi)
      end
      Vcoeff[slack, n+1] = 0.0 + 0.0im

      # Write back Q^(n-1) for PV buses into shifted storage column n
      @inbounds for k = 1:npv
         Qcoeff[k, n] = x[2*nn+k]
      end

      # Update inverse series W = 1/V with general germ V^(0)
      # w_n = -(1/v0) * Σ_{m=1..n} v_m w_{n-m}
      @inbounds for i = 1:nbus
         if i == slack
            Wcoeff[i, n+1] = 0.0 + 0.0im
            continue
         end
         acc = 0.0 + 0.0im
         for m = 1:n
            acc += Vcoeff[i, m+1] * Wcoeff[i, (n-m)+1]
         end
         Wcoeff[i, n+1] = -acc / Vcoeff[i, 1]
      end

      if debug && (n <= debug_maxn) && (n % debug_every == 0)
         dbg_print_row(n, rhs)
      end
   end

   # -----------------------
   # Evaluate at s=1 via Padé or series sum
   # -----------------------
   V = zeros(ComplexF64, nbus)
   V[slack] = Vslack

   Qpv = zeros(Float64, npv)

   eval_opts = isnothing(evaluation_options) ? APSLFEvaluationOptions(mode = use_pade ? :pade : :taylor) : evaluation_options
   @inbounds for i = 1:nbus
      i == slack && continue
      cV = @view Vcoeff[i, :]
      V[i] = evaluate_series(cV, eval_opts; collect_logs = false).voltage
   end
   if npv > 0
      @inbounds for k = 1:npv
         cQ = ComplexF64.(Qcoeff[k, :], 0.0)
         Qpv[k] = real(evaluate_series(cQ, eval_opts; collect_logs = false).voltage)
      end
   end

   # Optionally recompute PV Q from injections for consistency:
   #   Q_i = imag( S_inj(V)_i )
   if qpv_from_inj && npv > 0
      Sinj = calc_injections(Yuse, V)
      @inbounds for (k, bus) in enumerate(pv)
         Qpv[k] = imag(Sinj[bus])
      end
   end

   # -----------------------
   # Self-check (prints ONLY if something is off)
   #
   # Checks at the final evaluated V (s=1):
   # - P mismatch on all non-slack buses
   # - Q mismatch on PQ buses
   # - |V| mismatch on PV buses
   # - consistency between Qpv output and imag(Sinj) (optional)
   # -----------------------
   if self_check
      Sinj = calc_injections(Yuse, V)

      maxP = 0.0
      @inbounds for i = 1:nbus
         i == slack && continue
         maxP = max(maxP, abs(real(Sinj[i]) - Pspec[i]))
      end

      maxQpq = 0.0
      @inbounds for i = 1:nbus
         i == slack && continue
         _is_pq(bustype[i]) || continue
         maxQpq = max(maxQpq, abs(imag(Sinj[i]) - Qspec[i]))
      end

      maxVpv = 0.0
      @inbounds for bus in pv
         maxVpv = max(maxVpv, abs(abs(V[bus]) - Vm[bus]))
      end

      maxQpvDiff = 0.0
      if npv > 0
         @inbounds for (k, bus) in enumerate(pv)
            maxQpvDiff = max(maxQpvDiff, abs(Qpv[k] - imag(Sinj[bus])))
         end
      end

      if (maxP > mis_tol_p) || (maxQpq > mis_tol_q) || (maxVpv > v_tol) || (maxQpvDiff > qpv_tol)
         println("\n[apslf_pf_pv_direct] SELF-CHECK (only printed on violation)")
         @printf("  max |P mismatch| (p.u.)         = %.6e   (tol %.1e)\n", maxP, mis_tol_p)
         @printf("  max |Q mismatch| on PQ (p.u.)   = %.6e   (tol %.1e)\n", maxQpq, mis_tol_q)
         @printf("  max | |V|-Vsp | on PV (p.u.)    = %.6e   (tol %.1e)\n", maxVpv, v_tol)
         @printf("  max |Qpv_series - imag(Sinj)|   = %.6e   (tol %.1e)\n", maxQpvDiff, qpv_tol)
      end
   end

   return V, Qpv, Vcoeff, Wcoeff, Qcoeff
end

"""
    solve_pf_apslf_with_pv_q_limits(Y, bustype, Pspec, Qspec, Vm, Qmin, Qmax; kwargs...)
        -> NamedTuple

Outer-loop PV handling + Q-limits around an inner APSLF solver.

Inner modes:
- `inner = :pq`          # :pq | :direct_pv | :direct_pv_sparse
  Solve using PQ-only APSLF (`apslf_pq`) and enforce PV `|V|` by iteratively
  adjusting `Q` at PV buses (secant loop). Then enforce Q-limits (PV→PQ switching)
  in an outer loop.
- `inner = :direct_pv`:
  Use the direct PV kernel (`apslf_pf_pv_direct`) inside the same Q-limit outer loop.

Q-limit logic:
After a solve, for any PV bus `i`:
- `Qinj_i = imag(Sinj_i)`
- If `Qinj_i > Qmax_i + qtol`: switch PV→PQ, set `Q[i]=Qmax[i]`
- If `Qinj_i < Qmin_i - qtol`: switch PV→PQ, set `Q[i]=Qmin[i]`

Degenerate PV preprocessing:
- PV buses with `|Qmax - Qmin| <= qdeg_tol` are treated as fixed-Q and demoted to PQ
  before solving. This tolerance is intentionally independent from `qtol`.

Set `enforce_q_limits=false` to disable PV→PQ switching while still reporting
computed PV reactive injections in `Q`.

Convergence summary:
- For `inner=:pq`: converged if PV `|V|` errors `< vtol` and no switching occurs; when `nr_polish_Y` is active, post-polish P/Q mismatches on that Y-bus must also satisfy mismatch tolerances.
- For `inner=:direct_pv`: converged if mismatches (`maxP/maxQpq`) are below tolerances
  and no switching occurs.

APSLF germ semantics (`germ` keyword):
- `germ = :flat` (default): canonical analytic germ `V(s=0)=1∠0`. Exact at
  order 0 only for a pure series network with `Vslack = 1`; with line shunts,
  transformer taps or phase shifters the result needs the NR polish.
- `germ = :noload`: the germ is the linear no-load solution of the
  full Y-bus (theory Section 6.5, variant 2). Exact at order 0 for shunts,
  off-nominal ratios, phase-shifting transformers and `Vslack ≠ 1`; the pure
  APSLF result is then a load-flow solution without NR polish.
- `germ = :deviation` (default): flat germ `Vslack·1` with the deviation
  embedding `Y(s) = Y0 + s (Y - Y0)`, `Y0 = Y - diag(Y·1)` (theory Section
  6.5, variant 1). Also exact, but with a different path in `s`; it keeps the
  germ at nominal voltage and is the robust choice for large networks whose
  no-load state is far from the operating point (PEGASE cases converge with
  `:deviation` and diverge with `:noload`).
- The germ is not a Newton-Raphson start value. Legacy `flatstart`/`V0_germ`
  kwargs are deprecated and ignored with a warning.

Optional NR polish:
After each inner solve, an NR refinement may be applied to the final `V`
starting from the APSLF solution itself. Pass `nr_polish_Y` to use a different
admittance matrix for this NR refinement (for example, APSLF without bus shunts
followed by NR polish with the full Y-bus). When this alternate matrix is used,
the post-polish injections and convergence checks also use it so the returned
`converged` flag is consistent with the polished state.
"""
function solve_pf_apslf_with_pv_q_limits(
   Y::AbstractMatrix{ComplexF64},
   bustype::Vector{Symbol},
   Pspec::Vector{Float64},
   Qspec::Vector{Float64},
   Vm::Vector{Float64},
   Qmin::Vector{Float64},
   Qmax::Vector{Float64};
   slack::Int = 1,
   Vslack::Union{Nothing,ComplexF64} = nothing,
   order::Int = 24,
   use_pade::Bool = false,
   germ::Symbol = :deviation,    # :deviation | :noload | :flat
   evaluation_options::Union{Nothing,APSLFEvaluationOptions} = nothing,
   inner::Symbol = :pq,          # :pq | :direct_pv
   max_outer::Int = 30,
   max_pv_iter::Int = 12,
   vtol::Float64 = 1e-6,
   qtol::Float64 = 1e-6,
   qdeg_tol::Float64 = 1e-12,
   enforce_q_limits::Bool = true,
   q_limit_switch_min_outer::Int = 1,
   q_limit_switch_require_stable::Bool = false,
   q_limit_switch_stability_tol::Float64 = 1e-4,
   pv_step0::Float64 = 0.05,
   pv_secant_damping::Float64 = 1.0,
   nr_polish::Bool = true,
   nr_polish_Y::Union{Nothing,AbstractMatrix{ComplexF64}} = nothing,
   nr_q_polish::Bool = false,
   nr_max_iter::Int = 15,
   nr_tol::Float64 = 1e-10,
   nr_damping::Float64 = 1.0,
   nr_polish_autodamping::Bool = true,
   nr_polish_min_step::Float64 = 0.03125,
   nr_polish_max_damping_trials::Int = 6,
   nr_polish_fallback_to_apslf::Bool = true,
   nr_polish_allow_worse_result::Bool = false,
   nr_polish_max_voltage_pu::Float64 = 2.0,
   mis_tol_p::Float64 = 1e-8,
   mis_tol_q::Float64 = 1e-8,
   verbose::Int = 0,
   return_coeffs::Bool = false,
   use_sparse::Union{Bool,Symbol} = :auto,
   sparse_nbus_min::Int = 110,
   timeout_s::Real = 0.0,
   kwargs...,
)
   _handle_deprecated_apslf_germ_kwargs(kwargs; context = "solve_pf_apslf_with_pv_q_limits")

   solve_start_ns = time_ns()
   timeout_s_float = Float64(timeout_s)
   timeout_ns = timeout_s_float > 0.0 ? UInt64(round(timeout_s_float * 1e9)) : UInt64(0)
   function check_solver_timeout!(where::Symbol)
      timeout_ns == 0 && return nothing
      elapsed_ns = time_ns() - solve_start_ns
      elapsed_ns <= timeout_ns && return nothing
      throw(APSLFTimeoutError(where, elapsed_ns / 1e9, timeout_s_float))
   end
   check_solver_timeout!(:start)

   nbus = length(bustype)
   @assert size(Y, 1) == nbus && size(Y, 2) == nbus
   @assert length(Pspec) == nbus && length(Qspec) == nbus && length(Vm) == nbus
   @assert length(Qmin) == nbus && length(Qmax) == nbus
   q_limit_switch_min_outer >= 1 || throw(ArgumentError("q_limit_switch_min_outer must be at least 1."))
   q_limit_switch_stability_tol >= 0.0 || throw(ArgumentError("q_limit_switch_stability_tol must be non-negative."))
   nr_polish_min_step > 0.0 || throw(ArgumentError("nr_polish_min_step must be positive."))
   nr_polish_max_damping_trials >= 1 || throw(ArgumentError("nr_polish_max_damping_trials must be at least 1."))
   nr_polish_max_voltage_pu > 0.0 || throw(ArgumentError("nr_polish_max_voltage_pu must be positive."))
   if nr_polish_Y !== nothing
      @assert size(nr_polish_Y, 1) == nbus && size(nr_polish_Y, 2) == nbus
   end
   @assert inner in (:pq, :direct_pv, :direct_pv_sparse) "inner must be :pq | :direct_pv | :direct_pv_sparse"
   germ in (:flat, :noload, :deviation) || throw(ArgumentError("germ must be :noload, :deviation or :flat, got :$(germ)"))
   qdeg_tol >= 0.0 || throw(ArgumentError("qdeg_tol must be non-negative."))
   _warn_degenerate_q_limits(bustype, Qmin, Qmax; atol = qdeg_tol)
   Yuse = maybe_sparse_Y(Y; nbus = nbus, use_sparse = use_sparse, sparse_nbus_min = sparse_nbus_min)
   Ynr =
      nr_polish_Y === nothing ? Yuse :
      maybe_sparse_Y(nr_polish_Y; nbus = nbus, use_sparse = use_sparse, sparse_nbus_min = sparse_nbus_min)
   use_nr_polish_y_for_post_metrics = nr_polish && nr_max_iter > 0 && nr_polish_Y !== nothing
   Ypost = use_nr_polish_y_for_post_metrics ? Ynr : Yuse
   # Debug output for sparse policy
   if verbose >= 2
      y_type = issparse(Yuse) ? "sparse" : "dense"
      println("APSLF: Using $y_type Y matrix (nbus=$nbus, threshold=$sparse_nbus_min)")
   end


   bt = _canonicalize_bustype_lower(bustype)
   Q = copy(Qspec)
   demoted_pv = _demote_degenerate_pv_buses!(bt, Q, Qmin, Qmax; atol = qdeg_tol)
   if !isempty(demoted_pv)
      @warn "Demoting PV buses with Qmin ≈ Qmax to PQ before solving." buses = demoted_pv
   end

   make_S(P::Vector{Float64}, Qv::Vector{Float64}) = ComplexF64.(P, Qv)
   Vsl = (Vslack === nothing) ? ComplexF64(Vm[slack], 0.0) : Vslack

   all_switch_log =
      Vector{NamedTuple{(:outer, :bus, :qinj, :qmin, :qmax, :side),Tuple{Int,Int,Float64,Float64,Float64,Symbol}}}()

   # Keep the last coefficient matrix we computed (only meaningful if return_coeffs=true)
   Vcoeff_last = return_coeffs ? zeros(ComplexF64, nbus, order + 1) : nothing
   # APSLF PQ workspace (reuses Yred LU across many inner=:pq solves)
   ws_pq =
      (inner == :pq) ?
      build_apslf_pq_workspace(
         Yuse;
         slack = slack,
         order = order,
         use_sparse = false, # already applied policy above, TODO: check comment.
         sparse_nbus_min = sparse_nbus_min,
         germ = germ,
      ) : nothing

   # -----------------------
   # Optional NR polishing step (outer) - rectangular analytic Jacobian
   # Cached across outer iterations (rebuild only if bus-type layout for polish changes)
   # -----------------------
   nr_cache_ref = Ref{Union{Nothing,NRRectCache}}(nothing)
   nr_cache_sig = Ref{UInt64}(0)
   nr_polish_rejected = Ref(false)
   nr_polish_reject_reason = Ref(:none)
   nr_polish_score_before = Ref(NaN)
   nr_polish_score_after = Ref(NaN)
   nr_polish_rejection_warned = Ref(false)
   nr_polish_success = Ref(false)
   nr_polish_improved = Ref(false)
   nr_polish_damped = Ref(false)
   nr_polish_fallback_to_apslf_used = Ref(false)
   nr_polish_failed_no_acceptable_step = Ref(false)
   q_limit_switch_deferred = Ref(false)
   q_limit_switch_deferred_count = Ref(0)

   function voltage_max_abs(V::Vector{ComplexF64})
      _all_finite_complex(V) || return Inf
      return maximum(abs, V)
   end

   function nr_mismatch_score(cache::NRRectCache, V::Vector{ComplexF64}, Qwork::Vector{Float64})
      _all_finite_complex(V) || return Inf
      F = mismatch_rectangular!(cache, Ynr, V, Pspec, Qwork, Vm)
      all(isfinite, F) || return Inf
      return isempty(F) ? 0.0 : maximum(abs, F)
   end

   function valid_nr_polish_state(Vtrial::Vector{ComplexF64})
      _all_finite_complex(Vtrial) || return false
      vmax = voltage_max_abs(Vtrial)
      return isfinite(vmax) && vmax <= nr_polish_max_voltage_pu
   end

   function damping_factors()
      factors = Float64[]
      step = 1.0
      for _ = 1:nr_polish_max_damping_trials
         push!(factors, step)
         step <= nr_polish_min_step && break
         step = max(step / 2.0, nr_polish_min_step)
      end
      return factors
   end

   function nr_polish!(V::Vector{ComplexF64}, Qwork::Vector{Float64})
      nr_polish || return V
      nr_max_iter <= 0 && return V
      nr_polish_rejected[] = false
      nr_polish_reject_reason[] = :none
      nr_polish_success[] = false
      nr_polish_improved[] = false
      nr_polish_damped[] = false
      nr_polish_fallback_to_apslf_used[] = false
      nr_polish_failed_no_acceptable_step[] = false
      v_apslf = copy(V)
      # If nr_q_polish: treat PV as PQ during polish (i.e., ΔQ instead of ΔV constraints)
      # This is only meaningful for inner=:pq (outer PV enforcement); for :direct_pv it typically shouldn't be used.
      bt_polish = bt
      if nr_q_polish && inner == :pq
         bt_polish = copy(bt)
         @inbounds for i in axes(bt_polish, 1)
            bt_polish[i] == :pv && (bt_polish[i] = :pq)
         end
      end

      # Build/reuse cache (rebuild only if bt_polish layout changed)
      sig = _bt_sig(bt_polish)
      cache = nr_cache_ref[]
      if cache === nothing || sig != nr_cache_sig[]
         cache = build_nr_rect_cache(bt_polish, slack)
         nr_cache_ref[] = cache
         nr_cache_sig[] = sig
      end

      nr_polish_score_before[] = nr_mismatch_score(cache, V, Qwork)
      @debug "NR polish initial mismatch" status = :nr_polish_initial_mismatch mismatch = nr_polish_score_before[]

      best_V = copy(V)
      best_score = nr_polish_score_before[]
      original_score = nr_polish_score_before[]
      margin = max(sqrt(eps(Float64)), nr_tol)
      factors = nr_polish_autodamping ? damping_factors() : [1.0]

      # Rectangular NR steps
      for iter = 1:nr_max_iter
         check_solver_timeout!(:nr_polish_iteration)
         current_score = nr_mismatch_score(cache, V, Qwork)
         isfinite(current_score) || break
         accepted = false
         accepted_dx_norm = Inf
         accepted_residual_norm = current_score

         for factor in factors
            trial_V = copy(V)
            actual_damping = nr_damping * factor
            dxn, rn = try
               nr_refine_step_rect!(
                  cache,
                  Ynr,
                  trial_V,
                  Pspec,
                  Qwork,
                  Vm;
                  Vslack = Vsl,
                  damping = actual_damping,
                  use_sparse = (Ynr isa SparseMatrixCSC{ComplexF64}),
               )
            catch err
               _is_recoverable_linear_solve_error(err) || rethrow(err)
               (Inf, Inf)
            end
            trial_score = valid_nr_polish_state(trial_V) ? nr_mismatch_score(cache, trial_V, Qwork) : Inf
            acceptable = trial_score <= current_score * (1.0 + margin)
            if acceptable || nr_polish_allow_worse_result
               V .= trial_V
               accepted = true
               accepted_dx_norm = dxn * abs(actual_damping)
               accepted_residual_norm = rn
               nr_polish_damped[] |= nr_polish_autodamping && factor < 1.0
               @debug "NR polish accepted damping factor" status = :nr_polish_damped iteration = iter damping_factor =
                  actual_damping mismatch_before = current_score mismatch_after = trial_score
               if trial_score <= best_score
                  best_score = trial_score
                  best_V .= trial_V
               end
               break
            elseif nr_polish_autodamping
               @debug "NR polish rejected damped trial" status = :nr_polish_rejected_step iteration = iter damping_factor =
                  actual_damping mismatch_before = current_score mismatch_after = trial_score
            end
         end

         if !accepted
            nr_polish_failed_no_acceptable_step[] = true
            nr_polish_rejected[] = true
            nr_polish_reject_reason[] = :nr_polish_failed_no_acceptable_step
            V .= best_V
            if !nr_polish_rejection_warned[]
               @warn "Rejecting NR polish step because no damping factor improved the APSLF candidate" status =
                  :nr_polish_failed_no_acceptable_step iteration = iter mismatch_before = current_score best_mismatch =
                  best_score
               nr_polish_rejection_warned[] = true
            else
               @debug "Rejecting NR polish step because no damping factor improved the APSLF candidate" status =
                  :nr_polish_failed_no_acceptable_step iteration = iter mismatch_before = current_score best_mismatch =
                  best_score
            end
            break
         end
         max(accepted_dx_norm, accepted_residual_norm) < nr_tol && break
      end

      nr_polish_score_after[] = nr_mismatch_score(cache, V, Qwork)
      final_vmax = voltage_max_abs(V)
      final_finite = _all_finite_complex(V) && isfinite(final_vmax) && isfinite(nr_polish_score_after[])
      final_voltage_ok = final_vmax <= nr_polish_max_voltage_pu
      final_valid = final_finite && final_voltage_ok
      final_allowed = nr_polish_allow_worse_result || (nr_polish_score_after[] <= original_score * (1.0 + margin))
      if !(final_valid && final_allowed)
         if nr_polish_fallback_to_apslf
            V .= v_apslf
            nr_polish_score_after[] = original_score
            nr_polish_fallback_to_apslf_used[] = true
            nr_polish_reject_reason[] =
               !final_finite ? :nonfinite_or_infinite_nr_polish :
               (!final_voltage_ok ? :nr_polish_voltage_explosion : :nr_polish_residual_worsened)
         else
            V .= best_V
            nr_polish_score_after[] = best_score
            nr_polish_reject_reason[] = :nr_polish_restored_best_valid_state
         end
         nr_polish_rejected[] = true
         if !nr_polish_rejection_warned[]
            @warn "Rejecting NR polish because it worsened the APSLF candidate" status = nr_polish_reject_reason[] score_before =
               original_score score_after = nr_polish_score_after[] reason = nr_polish_reject_reason[]
            nr_polish_rejection_warned[] = true
         end
      end

      nr_polish_improved[] = nr_polish_score_after[] < original_score
      nr_polish_success[] =
         final_valid &&
         (nr_polish_improved[] || nr_polish_score_after[] <= max(nr_tol, eps(Float64))) &&
         !nr_polish_fallback_to_apslf_used[]
      @debug "NR polish final mismatch" status =
         nr_polish_success[] ? :nr_polish_success :
         (nr_polish_improved[] ? :nr_polish_improved : :nr_polish_fallback_to_apslf) mismatch = nr_polish_score_after[] improved =
         nr_polish_improved[] damped = nr_polish_damped[] fallback = nr_polish_fallback_to_apslf_used[]

      return V
   end

   # -----------------------
   # Inner solve wrapper:
   # Returns always (Vvec, Sinj_loc, Vcoeff_full_or_nothing)
   # -----------------------
   function inner_solve_given_Q(Qwork::Vector{Float64}; want_coeffs::Bool = false)
      if inner == :pq
         S = make_S(Pspec, Qwork)

         Vvec, Vcoeff_red, _ = apslf_pq_solve!(
            ws_pq,
            S;
            Vslack = Vsl,
            use_pade = use_pade,
            germ = germ,
            evaluation_options = evaluation_options,
            timeout_check = check_solver_timeout!,
         )

         _assert_finite_vec(Vvec, "inner voltage state V")
         want_coeffs && _assert_finite_coeff_matrix(Vcoeff_red, "inner voltage coefficients")

         Sinj_loc = calc_injections(Yuse, Vvec)

         if want_coeffs
            Vcoeff_full = zeros(ComplexF64, nbus, order + 1)
            Vcoeff_full[slack, 1] = Vsl
            nonslack = [i for i = 1:nbus if i != slack]
            @inbounds for (idx, bus) in enumerate(nonslack)
               Vcoeff_full[bus, :] .= Vcoeff_red[idx, :]
            end
            return Vvec, Sinj_loc, Vcoeff_full
         else
            return Vvec, Sinj_loc, nothing
         end

      elseif inner == :direct_pv
         Vvec, _, Vcoeff_full, _, _ = apslf_pf_pv_direct(
            Yuse,
            bt,
            Pspec,
            Qwork,
            Vm;
            slack = slack,
            Vslack = Vsl,
            order = order,
            use_pade = use_pade,
            germ = germ,
            evaluation_options = evaluation_options,
            debug = false,
            self_check = false,
            qpv_from_inj = true,
            timeout_check = check_solver_timeout!,
         )
         _assert_finite_vec(Vvec, "inner voltage state V")
         want_coeffs && _assert_finite_coeff_matrix(Vcoeff_full, "inner voltage coefficients")
         Sinj_loc = calc_injections(Yuse, Vvec)
         return want_coeffs ? (Vvec, Sinj_loc, Vcoeff_full) : (Vvec, Sinj_loc, nothing)

      else
         # :direct_pv_sparse
         Vvec, _, Vcoeff_full, _, _ = apslf_pf_pv_direct_sparse(
            Yuse,
            bt,
            Pspec,
            Qwork,
            Vm;
            slack = slack,
            Vslack = Vsl,
            order = order,
            use_pade = use_pade,
            germ = germ,
            evaluation_options = evaluation_options,
            debug = false,
            self_check = false,
            qpv_from_inj = true,
            timeout_check = check_solver_timeout!,
         )
         _assert_finite_vec(Vvec, "inner voltage state V")
         want_coeffs && _assert_finite_coeff_matrix(Vcoeff_full, "inner voltage coefficients")
         Sinj_loc = calc_injections(Yuse, Vvec)
         return want_coeffs ? (Vvec, Sinj_loc, Vcoeff_full) : (Vvec, Sinj_loc, nothing)
      end
   end
   # -----------------------
   # Mismatch metrics (used for convergence checks)
   # -----------------------
   function mismatch_metrics(V::Vector{ComplexF64}, Qwork::Vector{Float64})
      # reuse outer buffers (Sinj/Iinj) to avoid allocations
      calc_injections!(Sinj, Iinj, Ypost, V)

      maxP = 0.0
      maxQpq = 0.0
      @inbounds for i = 1:nbus
         i == slack && continue
         maxP = max(maxP, abs(real(Sinj[i]) - Pspec[i]))
         if bt[i] == :pq
            maxQpq = max(maxQpq, abs(imag(Sinj[i]) - Qwork[i]))
         end
      end
      return maxP, maxQpq
   end
   # -----------------------
   # Verbose final report (only prints if verbose>=1)
   # -----------------------
   function print_final_report(converged::Bool, outer_iters::Int, V::Vector{ComplexF64}, Qwork::Vector{Float64})
      verbose >= 1 || return

      maxP, maxQpq = mismatch_metrics(V, Qwork)

      max_v_err = 0.0
      if inner == :pq
         @inbounds for i = 1:nbus
            bt[i] == :pv || continue
            max_v_err = max(max_v_err, abs(abs(V[i]) - Vm[i]))
         end
      end

      npv_now = count(==(:pv), bt)
      npq_now = count(==(:pq), bt)

      println("\n--- APSLF Final Report ---")
      @printf("converged           = %s\n", string(converged))
      @printf("inner               = %s\n", string(inner))
      @printf("outer_iters         = %d\n", outer_iters)
      @printf("nbus                = %d\n", nbus)
      @printf("slack               = %d\n", slack)
      @printf("final npv / npq     = %d / %d\n", npv_now, npq_now)
      @printf("maxP mismatch (pu)  = %.6e  (tol %.1e)\n", maxP, mis_tol_p)
      @printf("maxQpq mismatch(pu) = %.6e  (tol %.1e)\n", maxQpq, mis_tol_q)
      if inner == :pq
         @printf("max PV |V|-Vsp (pu) = %.6e  (vtol %.1e)\n", max_v_err, vtol)
      end
      @printf("PV→PQ switches      = %d\n", length(all_switch_log))

      if !isempty(all_switch_log)
         if verbose >= 2
            println("\nSwitch events (PV→PQ):")
            println("  outer  bus  side   Qinj(pu)        Qmin(pu)        Qmax(pu)")
            for ev in all_switch_log
               @printf(
                  "  %5d  %3d  %-4s  %+12.6e  %+12.6e  %+12.6e\n",
                  ev.outer,
                  ev.bus,
                  string(ev.side),
                  ev.qinj,
                  ev.qmin,
                  ev.qmax
               )
            end
         else
            println("Switch events suppressed at verbose=1; use verbose=2 to print every PV→PQ event.")
         end
      end

      println("--- End Final Report ---\n")
   end

   function make_result(; V, Sinj, bt, Q, converged, outer_iters, reason = nothing, Vcoeff = nothing)
      base = (
         V = V,
         Sinj = Sinj,
         bustype = bt,
         Q = Q,
         switch_log = copy(all_switch_log),
         converged = converged,
         outer_iters = outer_iters,
         apslf_germ = germ,
         nr_polish_enabled = nr_polish,
         nr_polish_start = nr_polish ? :apslf_solution : :none,
         nr_polish_rejected = nr_polish_rejected[],
         nr_polish_reject_reason = nr_polish_reject_reason[],
         nr_polish_score_before = nr_polish_score_before[],
         nr_polish_score_after = nr_polish_score_after[],
         nr_polish_success = nr_polish_success[],
         nr_polish_improved = nr_polish_improved[],
         nr_polish_damped = nr_polish_damped[],
         nr_polish_fallback_to_apslf = nr_polish_fallback_to_apslf_used[],
         nr_polish_failed_no_acceptable_step = nr_polish_failed_no_acceptable_step[],
         q_limit_switch_deferred = q_limit_switch_deferred[],
         q_limit_switch_deferred_count = q_limit_switch_deferred_count[],
      )
      with_reason = (reason === nothing) ? base : merge(base, (reason = reason,))
      return (Vcoeff === nothing) ? with_reason : merge(with_reason, (Vcoeff = Vcoeff,))
   end

   # -----------------------
   # Main outer loop
   # -----------------------
   V = zeros(ComplexF64, nbus)
   Sinj = zeros(ComplexF64, nbus)
   Iinj = zeros(ComplexF64, nbus)

   for outer = 1:max_outer
      check_solver_timeout!(:outer_iteration)
      pv_buses = findall(i -> bt[i] == :pv, 1:nbus)

      # Step A: enforce |V| at PV buses by adjusting Q (secant per PV) for inner=:pq only
      if inner == :pq
         # Secant step damping (kwarg). The inner APSLF solve is exact for the
         # given Q, so the undamped secant (1.0) converges superlinearly; a
         # damped value (e.g. 0.5) only slows convergence and can leave the PV
         # |V| error above vtol within max_pv_iter.
         pv_secant_damping > 0.0 || throw(ArgumentError("pv_secant_damping must be positive."))
         q_probe_margin = 0.1
         overflow_f_limit = 1e6
         for i in pv_buses
            check_solver_timeout!(:pv_secant_bus)
            bt[i] != :pv && continue

            Qi0 = Q[i]
            Qi1 = Qi0 + pv_step0

            Qtmp = copy(Q)
            Qtmp[i] = Qi0
            V0, _, _ = try
               inner_solve_given_Q(Qtmp; want_coeffs = false)
            catch err
               _is_nonfinite_inner_error(err) || rethrow(err)
               # Some hard cases can trigger non-finite intermediate states for a
               # particular PV probe value. Skip secant enforcement on this bus in
               # this outer step and keep current Q.
               continue
            end
            f0 = abs(V0[i]) - Vm[i]
            (isfinite(f0) && abs(f0) <= overflow_f_limit) || continue
            abs(f0) < vtol && continue

            # If the default step is too aggressive, back off a few times.
            got_probe = false
            V1 = similar(V0)
            f1 = 0.0
            for _ = 1:4
               check_solver_timeout!(:pv_secant_probe)
               Qtmp[i] = Qi1
               ok = true
               Vcand = V0
               try
                  Vcand, _, _ = inner_solve_given_Q(Qtmp; want_coeffs = false)
               catch err
                  _is_nonfinite_inner_error(err) || rethrow(err)
                  ok = false
               end
               if ok
                  V1 = Vcand
                  f1 = abs(V1[i]) - Vm[i]
                  (isfinite(f1) && abs(f1) <= overflow_f_limit) || begin
                     got_probe = false
                     break
                  end
                  got_probe = true
                  break
               end
               Qi1 = Qi0 + 0.5 * (Qi1 - Qi0)
            end
            got_probe || continue

            Qi = Qi1
            for _ = 1:max_pv_iter
               check_solver_timeout!(:pv_secant_iteration)
               abs(f1 - f0) < 1e-12 && break

               Qi2_raw = Qi1 - f1 * (Qi1 - Qi0) / (f1 - f0)
               Qi2 = Qi1 + pv_secant_damping * (Qi2_raw - Qi1)
               Qi2 = clamp(Qi2, Qmin[i] - q_probe_margin, Qmax[i] + q_probe_margin)

               Qtmp[i] = Qi2
               V2, _, _ = try
                  inner_solve_given_Q(Qtmp; want_coeffs = false)
               catch err
                  _is_nonfinite_inner_error(err) || rethrow(err)
                  break
               end
               f2 = abs(V2[i]) - Vm[i]
               (isfinite(f2) && abs(f2) <= overflow_f_limit) || break

               Qi = Qi2
               abs(f2) < vtol && break

               Qi0, f0 = Qi1, f1
               Qi1, f1 = Qi2, f2
            end

            Q[i] = Qi
         end
      end

      # Step B: inner solve with current Q
      Vcoeff_opt = nothing
      try
         V, Sinj, Vcoeff_opt = inner_solve_given_Q(Q; want_coeffs = return_coeffs)
      catch err
         if _is_nonfinite_inner_error(err)
            # Fallback for hard cases: when the PQ inner solve becomes non-finite,
            # try one direct-PV inner evaluation before declaring failure.
            if inner == :pq
               fallback_inner = (nbus >= sparse_nbus_min) ? :direct_pv_sparse : :direct_pv
               prev_inner = inner
               local Vfb, Sfb, Cfb
               try
                  inner = fallback_inner
                  Vfb, Sfb, Cfb = inner_solve_given_Q(Q; want_coeffs = return_coeffs)
                  V .= Vfb
                  Sinj .= Sfb
                  Vcoeff_opt = Cfb
               catch fb_err
                  _is_nonfinite_inner_error(fb_err) || rethrow(fb_err)
                  inner = prev_inner
                  print_final_report(false, outer, V, Q)
                  return make_result(
                     V = V,
                     Sinj = Sinj,
                     bt = bt,
                     Q = Q,
                     converged = false,
                     outer_iters = outer,
                     reason = :nonfinite_inner_state,
                     Vcoeff = return_coeffs ? Vcoeff_last : nothing,
                  )
               end
               inner = prev_inner
               # Continue with fallback state.
            else
               print_final_report(false, outer, V, Q)
               return make_result(
                  V = V,
                  Sinj = Sinj,
                  bt = bt,
                  Q = Q,
                  converged = false,
                  outer_iters = outer,
                  reason = :nonfinite_inner_state,
                  Vcoeff = return_coeffs ? Vcoeff_last : nothing,
               )
            end
         end
         rethrow(err)
      end
      if return_coeffs
         Vcoeff_last .= Vcoeff_opt
      end

      # Step B2: NR polish (outer)
      nr_polish!(V, Q)
      if !_all_finite_complex(V)
         print_final_report(false, outer, V, Q)
         return make_result(
            V = V,
            Sinj = Sinj,
            bt = bt,
            Q = Q,
            converged = false,
            outer_iters = outer,
            reason = :nonfinite_voltage_after_nr,
            Vcoeff = return_coeffs ? Vcoeff_last : nothing,
         )
      end
      calc_injections!(Sinj, Iinj, Ypost, V)

      maxP, maxQpq = mismatch_metrics(V, Q)

      max_v_err = 0.0
      if inner == :pq
         @inbounds for i in axes(bt, 1)
            bt[i] == :pv || continue
            max_v_err = max(max_v_err, abs(abs(V[i]) - Vm[i]))
         end
      end

      q_switch_stability_score = inner == :pq ? max(maxP, maxQpq, max_v_err) : max(maxP, maxQpq)
      q_switch_stable_enough = !q_limit_switch_require_stable || q_switch_stability_score <= q_limit_switch_stability_tol
      q_switch_outer_allowed = outer >= q_limit_switch_min_outer && q_switch_stable_enough

      # Step C: Q-limit enforcement (PV→PQ switching), optional
      switched = false
      q_limit_violation_pending = false
      for i in axes(bt, 1)
         if bt[i] == :pv
            Qi_calc = imag(Sinj[i])
            if enforce_q_limits
               if Qi_calc > Qmax[i] + qtol
                  if q_switch_outer_allowed
                     push!(
                        all_switch_log,
                        (outer = outer, bus = i, qinj = Qi_calc, qmin = Qmin[i], qmax = Qmax[i], side = :max),
                     )
                     Q[i] = Qmax[i]
                     bt[i] = :pq
                     switched = true
                  else
                     q_limit_violation_pending = true
                     q_limit_switch_deferred[] = true
                     q_limit_switch_deferred_count[] += 1
                     Q[i] = Qi_calc
                     @debug "Deferring PV→PQ Q-limit switch" status = :q_limit_switch_deferred outer = outer bus = i side =
                        :max qinj = Qi_calc qmax = Qmax[i] stability_score = q_switch_stability_score
                  end
               elseif Qi_calc < Qmin[i] - qtol
                  if q_switch_outer_allowed
                     push!(
                        all_switch_log,
                        (outer = outer, bus = i, qinj = Qi_calc, qmin = Qmin[i], qmax = Qmax[i], side = :min),
                     )
                     Q[i] = Qmin[i]
                     bt[i] = :pq
                     switched = true
                  else
                     q_limit_violation_pending = true
                     q_limit_switch_deferred[] = true
                     q_limit_switch_deferred_count[] += 1
                     Q[i] = Qi_calc
                     @debug "Deferring PV→PQ Q-limit switch" status = :q_limit_switch_deferred outer = outer bus = i side =
                        :min qinj = Qi_calc qmin = Qmin[i] stability_score = q_switch_stability_score
                  end
               else
                  # keep Q consistent with injections
                  Q[i] = Qi_calc
               end
            else
               # Keep PV buses in PV mode and track computed reactive output.
               Q[i] = Qi_calc
            end
         end
      end

      # Step D: convergence check
      if verbose >= 2
         npv_now = count(==(:pv), bt)
         npq_now = count(==(:pq), bt)
         if inner == :pq
            @printf(
               "[outer %2d] inner=%s  switched=%s  npv=%d npq=%d  max_v_err=%.3e  maxP=%.3e  maxQpq=%.3e\n",
               outer,
               string(inner),
               string(switched),
               npv_now,
               npq_now,
               max_v_err,
               maxP,
               maxQpq
            )
         else
            @printf(
               "[outer %2d] inner=%s  switched=%s  npv=%d npq=%d  maxP=%.3e  maxQpq=%.3e\n",
               outer,
               string(inner),
               string(switched),
               npv_now,
               npq_now,
               maxP,
               maxQpq
            )
         end
      end

      # Return on success
      if inner == :pq
         pq_mismatch_ok = !use_nr_polish_y_for_post_metrics || ((maxP < mis_tol_p) && (maxQpq < mis_tol_q))
         if max_v_err < vtol && !switched && !q_limit_violation_pending && pq_mismatch_ok
            print_final_report(true, outer, V, Q)
            return make_result(
               V = V,
               Sinj = Sinj,
               bt = bt,
               Q = Q,
               converged = true,
               outer_iters = outer,
               Vcoeff = return_coeffs ? Vcoeff_last : nothing,
            )
         end
      else
         if !switched && !q_limit_violation_pending && (maxP < mis_tol_p) && (maxQpq < mis_tol_q)
            print_final_report(true, outer, V, Q)
            return make_result(
               V = V,
               Sinj = Sinj,
               bt = bt,
               Q = Q,
               converged = true,
               outer_iters = outer,
               Vcoeff = return_coeffs ? Vcoeff_last : nothing,
            )
         end
      end
   end

   # Non-converged
   print_final_report(false, max_outer, V, Q)
   return make_result(
      V = V,
      Sinj = Sinj,
      bt = bt,
      Q = Q,
      converged = false,
      outer_iters = max_outer,
      Vcoeff = return_coeffs ? Vcoeff_last : nothing,
   )
end

"""
    solve_pf_apslf(spec; mode::Symbol=:direct, kwargs...) -> NamedTuple

Unified public entrypoint for APSLF (Holomorphic Embedding Load Flow) power flow solutions.

# Arguments
- `spec`: Power flow specification object containing system data with fields:
  - `Y::AbstractMatrix{ComplexF64}`: Nodal admittance matrix (nbus × nbus)
  - `bustype::Vector{Symbol}`: Bus types (:PQ, :PV, :Slack)
  - `Pspec::Vector{Float64}`: Specified active power injections (p.u.)
  - `Qspec::Vector{Float64}`: Specified reactive power injections (p.u.)
  - `Vm::Vector{Float64}`: Specified voltage magnitudes for PV buses (p.u.)
  - `Qmin::Vector{Float64}`: Minimum reactive power limits (p.u.)
  - `Qmax::Vector{Float64}`: Maximum reactive power limits (p.u.)
  - `slack::Int`: Slack bus index

# Keyword Arguments
- `mode::Symbol = :direct`: Solution method
  - `:direct` → Uses direct PV formulation (`apslf_pf_pv_direct`)
  - `:outer` → Uses PQ-only APSLF with outer PV loop (`apslf_pq` + secant)
- Additional kwargs passed to `solve_pf_apslf_with_pv_q_limits`

# Returns
`NamedTuple` containing:
- `V::Vector{ComplexF64}`: Final voltage solution (length nbus)
- `Sinj::Vector{ComplexF64}`: Complex power injections at solution
- `bustype::Vector{Symbol}`: Final bus types (may change due to Q-limit switching)
- `Q::Vector{Float64}`: Final reactive power schedule
- `converged::Bool`: Convergence flag
- `outer_iters::Int`: Number of outer iterations required
- `Vcoeff::Matrix{ComplexF64}`: Voltage series coefficients (if `return_coeffs=true`)

# Solution Methods

## Direct Mode (`:direct`)
- Uses `apslf_pf_pv_direct` kernel with augmented real system
- Solves PV buses directly in APSLF recursion via voltage magnitude constraints
- Generally more robust and faster for systems with many PV buses
- Handles PV buses natively without outer iteration loops

## Outer Mode (`:outer`)
- Uses `apslf_pq` for PQ-only APSLF + secant method for PV enforcement
- Iteratively adjusts Q at PV buses to meet voltage magnitude specifications
- Falls back method compatible with traditional APSLF implementations
- May require more iterations but provides step-by-step PV adjustment visibility

# Q-Limit Handling
Both modes include automatic PV↔PQ bus type switching:
- If PV bus reactive injection exceeds `Qmax + qtol`: switch to PQ at `Qmax`
- If PV bus reactive injection falls below `Qmin - qtol`: switch to PQ at `Qmin`
- Continues until no more switching occurs and power flow converges

This behavior can be disabled with `enforce_q_limits=false` (keeps PV buses as PV).
Switching can be delayed with `q_limit_switch_min_outer > 1`, and can be
gated on a stable mismatch/voltage-error score with
`q_limit_switch_require_stable=true` plus `q_limit_switch_stability_tol`.

# Slack Bus Voltage
Constructed automatically as `Vslack = Vm[slack] + 0im` (zero angle reference).
Can be overridden via kwargs to `solve_pf_apslf_with_pv_q_limits`.

# Series Evaluation
- Default uses Padé `[L/M]` approximants for improved convergence radius
- Series order and evaluation method configurable via kwargs
- Uses the canonical APSLF germ `V(s=0)=1∠0`

# Performance Notes
- Automatically selects sparse/dense matrix algorithms based on system size
- Reuses LU factorizations across iterations for efficiency
- Optional Newton-Raphson polishing for higher accuracy

# See Also
- `solve_pf_apslf_with_pv_q_limits`: Lower-level interface with full parameter control
- `apslf_pf_pv_direct`: Direct PV formulation kernel
- `apslf_pq`: PQ-only APSLF solver
"""
function solve_pf_apslf(spec; mode::Symbol = :direct, kwargs...)
   @assert mode == :direct || mode == :outer "mode must be :direct or :outer"

   Y = getfield(spec, :Y)
   nbus = size(Y, 1)
   kw = (; kwargs...)
   _handle_deprecated_apslf_germ_kwargs(kw; context = "solve_pf_apslf", strict = false)

   # default inner choice
   inner_auto = (mode == :direct) ? :direct_pv : :pq

   # auto-upgrade to sparse kernel when Y is sparse or system is large
   if mode == :direct
      if issparse(Y) || nbus >= 110
         inner_auto = :direct_pv_sparse
      end
   end
   inner_requested = haskey(kw, :inner) ? kw.inner : inner_auto
   kw_no_inner = haskey(kw, :inner) ? Base.structdiff(kw, NamedTuple{(:inner,)}((kw.inner,))) : kw
   kw_forward =
      haskey(kw_no_inner, :flatstart) ? Base.structdiff(kw_no_inner, NamedTuple{(:flatstart,)}((kw_no_inner.flatstart,))) :
      kw_no_inner
   verbose_level = haskey(kw_forward, :verbose) ? Int(kw_forward.verbose) : 0
   dense_fallback_nbus_max = haskey(kw_forward, :dense_fallback_nbus_max) ? Int(kw_forward.dense_fallback_nbus_max) : 2000
   kw_forward =
      haskey(kw_forward, :dense_fallback_nbus_max) ?
      Base.structdiff(kw_forward, NamedTuple{(:dense_fallback_nbus_max,)}((kw_forward.dense_fallback_nbus_max,))) : kw_forward
   outer_fallback_nbus_max = haskey(kw_forward, :outer_fallback_nbus_max) ? Int(kw_forward.outer_fallback_nbus_max) : 2000
   kw_forward =
      haskey(kw_forward, :outer_fallback_nbus_max) ?
      Base.structdiff(kw_forward, NamedTuple{(:outer_fallback_nbus_max,)}((kw_forward.outer_fallback_nbus_max,))) : kw_forward
   outer_fallback_pv_bus_max =
      haskey(kw_forward, :outer_fallback_pv_bus_max) ? Int(kw_forward.outer_fallback_pv_bus_max) : 250
   kw_forward =
      haskey(kw_forward, :outer_fallback_pv_bus_max) ?
      Base.structdiff(kw_forward, NamedTuple{(:outer_fallback_pv_bus_max,)}((kw_forward.outer_fallback_pv_bus_max,))) :
      kw_forward
   kw_forward =
      haskey(kw_forward, :V0_germ) ? Base.structdiff(kw_forward, NamedTuple{(:V0_germ,)}((kw_forward.V0_germ,))) : kw_forward

   Vslack = ComplexF64(getfield(spec, :Vm)[getfield(spec, :slack)], 0.0)
   res = solve_pf_apslf_with_pv_q_limits(
      getfield(spec, :Y),
      getfield(spec, :bustype),
      getfield(spec, :Pspec),
      getfield(spec, :Qspec),
      getfield(spec, :Vm),
      getfield(spec, :Qmin),
      getfield(spec, :Qmax);
      slack = getfield(spec, :slack),
      Vslack = Vslack,
      inner = inner_requested,
      kw_forward...,
   )

   # Robustness fallback path for direct mode:
   # 1) If sparse-direct does not converge, optionally retry with dense-direct kernel.
   # 2) If direct kernels still do not converge, retry with outer PQ loop.
   #
   # Dense-direct retry is skipped for large problems because it builds/factors
   # dense augmented systems that are not viable for large sparse networks.
   if mode == :direct && !res.converged
      if inner_requested == :direct_pv_sparse && nbus <= dense_fallback_nbus_max
         res_dense = solve_pf_apslf_with_pv_q_limits(
            getfield(spec, :Y),
            getfield(spec, :bustype),
            getfield(spec, :Pspec),
            getfield(spec, :Qspec),
            getfield(spec, :Vm),
            getfield(spec, :Qmin),
            getfield(spec, :Qmax);
            slack = getfield(spec, :slack),
            Vslack = Vslack,
            inner = :direct_pv,
            kw_forward...,
         )
         res_dense.converged && return merge(res_dense, (effective_mode = :direct,))
      end

      npv_for_outer = count(bt -> _is_pv(bt), getfield(spec, :bustype))
      if nbus <= outer_fallback_nbus_max && npv_for_outer <= outer_fallback_pv_bus_max
         res_outer = solve_pf_apslf_with_pv_q_limits(
            getfield(spec, :Y),
            getfield(spec, :bustype),
            getfield(spec, :Pspec),
            getfield(spec, :Qspec),
            getfield(spec, :Vm),
            getfield(spec, :Qmin),
            getfield(spec, :Qmax);
            slack = getfield(spec, :slack),
            Vslack = Vslack,
            inner = :pq,
            kw_forward...,
         )
         res_outer.converged && return merge(res_outer, (effective_mode = :outer,))
      elseif verbose_level >= 1
         @printf(
            "APSLF: skipping outer fallback because nbus=%d (limit %d), npv=%d (limit %d)\n",
            nbus,
            outer_fallback_nbus_max,
            npv_for_outer,
            outer_fallback_pv_bus_max
         )
      end
   end

   eff_mode = mode == :outer ? :outer : (inner_requested == :pq ? :outer : :direct)
   return merge(res, (effective_mode = eff_mode,))
end
