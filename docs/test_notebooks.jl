# =============================================================================
# File: docs/test_notebooks.jl
# Date: 2026-10-05
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Run every workshop notebook source docs/lit/*.jl as a script, each
#          in a fresh Julia process like a new Colab session, and fail on an
#          error or a warning. The sources carry `@assert ... #src` lines
#          that Literate drops from the notebooks and pages; here they run
#          and check the results against the statements of the text.
#          Run from the repository root: `julia --project=. docs/test_notebooks.jl`
#          (the large-network notebook downloads case2869pegase.m once).
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

const ROOT = normpath(joinpath(@__DIR__, ".."))
const LIT_DIR = joinpath(@__DIR__, "lit")

# Runs one source in its own process from the repository root (so the
# large-network notebook caches its download in data/_downloaded/).
# Returns (ok, seconds, output); a warning or error log entry fails the run.
function run_notebook(path::AbstractString)
   out = IOBuffer()
   cmd = Cmd(`$(Base.julia_cmd()) --project=$(ROOT) --startup-file=no $(path)`; dir = ROOT)
   t = @elapsed proc = run(pipeline(ignorestatus(cmd); stdout = out, stderr = out))
   text = String(take!(out))
   logged_problem = occursin("┌ Warning", text) || occursin("┌ Error", text)
   return (ok = success(proc) && !logged_problem, seconds = t, output = text)
end

function main()
   sources = filter(endswith(".jl"), sort(readdir(LIT_DIR)))
   failed = String[]
   for source in sources
      r = run_notebook(joinpath(LIT_DIR, source))
      println(rpad(source, 34), r.ok ? "ok" : "FAILED", "  ($(round(r.seconds, digits = 1)) s)")
      if !r.ok
         push!(failed, source)
         println("----- output of $(source) -----")
         println(r.output)
         println("----- end of $(source) -----")
      end
   end
   isempty(failed) || error("notebook sources failed: $(join(failed, ", "))")
   println("all $(length(sources)) notebook sources ran without errors or warnings")
end

main()
