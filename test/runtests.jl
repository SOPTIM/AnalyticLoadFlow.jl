# =============================================================================
# File: test/runtests.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Test-suite entry point that includes the maintained APSLF test sets.
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
using Logging
using AnalyticLoadFlow

const _SUITE_FILES = (
   ("Core Solver Functions", "test_solver_core.jl"),
   ("Performance and Benchmarks", "test_performance.jl"),
   ("Robustness and Edge Cases", "test_robustness.jl"),
   ("Minimal Y-Bus Example", "test_minimal_example.jl"),
   ("Germ and Inner Solver", "test_germ_fix.jl"),
)

function _render_progress(done::Int, total::Int, label::AbstractString = "")
   width = 28
   filled = total == 0 ? width : round(Int, width * done / total)
   filled = clamp(filled, 0, width)
   empty = width - filled
   pct = total == 0 ? 100 : round(Int, 100 * done / total)
   line = "[" * repeat('█', filled) * repeat('░', empty) * "] " * lpad(string(pct), 3) * "%"
   if !isempty(label)
      line *= "  " * label
   end
   print('\r', rpad(line, 96))
   flush(stdout)
end

# Silence the per-testset summaries while the progress bar runs. On Julia ≥ 1.13
# `Test.TESTSET_PRINT_ENABLE` is a `ScopedValue` and must be set with `with`;
# on older versions it is a `Ref`.
function _with_testset_print_disabled(f)
   flag = Test.TESTSET_PRINT_ENABLE
   if flag isa Base.ScopedValues.ScopedValue
      return Base.ScopedValues.with(f, flag => false)
   end
   prev = flag[]
   flag[] = false
   try
      return f()
   finally
      flag[] = prev
   end
end

test_results = _with_testset_print_disabled() do
   @testset "AnalyticLoadFlow.jl Complete Test Suite" begin
      total = length(_SUITE_FILES)
      _render_progress(0, total, "starting")
      for idx in eachindex(_SUITE_FILES)
         label, file = _SUITE_FILES[idx]
         @testset "$label" begin
            with_logger(Logging.NullLogger()) do
               redirect_stdout(devnull) do
                  redirect_stderr(devnull) do
                     include(file)
                  end
               end
            end
         end
         _render_progress(idx, total, "$label done")
      end
   end
end

print("\n")
if Test.get_test_counts(test_results).fails == 0 && Test.get_test_counts(test_results).errors == 0
   println("Result: PASS")
else
   println("Result: FAIL")
end

nothing

