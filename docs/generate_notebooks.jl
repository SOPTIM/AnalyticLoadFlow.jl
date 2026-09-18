# =============================================================================
# File: docs/generate_notebooks.jl
# Date: 2026-09-17
# Author: Udo Schmitz
# Organization: SOPTIM AG
# Purpose: Regenerate the committed Literate.jl outputs from docs/lit/*.jl: a
#          Documenter page under docs/src/generated/ and a Colab-ready notebook
#          under notebooks/. Run manually after editing a source:
#          `julia --project=docs docs/generate_notebooks.jl`. Both outputs are
#          committed; docs/make.jl does not call this script.
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

using Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()

using Literate
using JSON3

const LIT_DIR = joinpath(@__DIR__, "lit")
const GENERATED_DIR = joinpath(@__DIR__, "src", "generated")
const NOTEBOOK_DIR = normpath(joinpath(@__DIR__, "..", "notebooks"))

"""
    patch_notebook_metadata!(path)

Rewrite the notebook's `metadata.kernelspec` to the Julia 1.13 kernel
(`julia-1.13`, the version the package targets) and pin `language_info` to the
same version, so the committed file does not depend on the Julia that ran the
generator.
"""
function patch_notebook_metadata!(path::AbstractString)
   nb = copy(JSON3.read(read(path, String)))
   metadata = get!(nb, :metadata, Dict{Symbol,Any}())
   metadata[:kernelspec] = Dict(:name => "julia-1.13", :display_name => "Julia 1.13.0", :language => "julia")
   metadata[:language_info] = Dict(:name => "julia", :version => "1.13.0", :file_extension => ".jl", :mimetype => "application/julia")
   open(path, "w") do io
      JSON3.pretty(io, nb)
      println(io)
   end
   ks = JSON3.read(read(path, String)).metadata.kernelspec
   @assert ks.name == "julia-1.13" && ks.language == "julia" && ks.display_name == "Julia 1.13.0"
   return nothing
end

"""
    check_outputs(mdpath, nbpath)

The Colab install cell (`#nb` lines) must be present in the notebook and absent
from the Documenter page, where the docs project already provides the package.
"""
function check_outputs(mdpath::AbstractString, nbpath::AbstractString)
   nb = JSON3.read(read(nbpath, String))
   install_cell = any(occursin("Pkg.add(", join(cell.source)) for cell in nb.cells)
   @assert install_cell "Colab install cell missing from $(nbpath)"
   @assert !occursin("Pkg.add(", read(mdpath, String)) "install cell leaked into $(mdpath)"
   return nothing
end

function generate()
   mkpath(GENERATED_DIR)
   mkpath(NOTEBOOK_DIR)
   for source in filter(endswith(".jl"), sort(readdir(LIT_DIR)))
      srcpath = joinpath(LIT_DIR, source)
      stem = first(splitext(source))
      Literate.markdown(srcpath, GENERATED_DIR; credit = false)
      Literate.notebook(srcpath, NOTEBOOK_DIR; execute = false, credit = false)
      nbpath = joinpath(NOTEBOOK_DIR, stem * ".ipynb")
      mdpath = joinpath(GENERATED_DIR, stem * ".md")
      patch_notebook_metadata!(nbpath)
      check_outputs(mdpath, nbpath)
      println("generated: $(relpath(mdpath)) and $(relpath(nbpath))")
   end
end

Base.invokelatest(generate)
