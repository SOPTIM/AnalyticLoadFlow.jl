# =============================================================================
# File: docs/make.jl
# Date: 2026-07-08
# Author: Udo Schmitz
# Organization: SOPTIM AG
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

# Purpose: Builds the public Documenter.jl documentation from checked-in Markdown sources.

# file: docs/make.jl

using Documenter
using AnalyticLoadFlow

let
   docs_dir = @__DIR__
   docs_src = joinpath(docs_dir, "src")
   assets_dir = joinpath(docs_src, "assets")
   repo_root = normpath(joinpath(docs_dir, ".."))

   mkpath(docs_src)
   mkpath(assets_dir)

   # Copy README into docs/src so it becomes the landing page.
   readme_in = joinpath(repo_root, "README.md")
   readme_out = joinpath(docs_src, "README.md")
   index_out = joinpath(docs_src, "index.md")
   changelog_in = joinpath(repo_root, "CHANGELOG.md")
   changelog_out = joinpath(docs_src, "CHANGELOG.md")

   if isfile(readme_in)
      cp(readme_in, readme_out; force = true)
      cp(readme_in, index_out; force = true)   # <- makes index.html possible
   else
      @warn "Missing README.md at repo root: $readme_in"
   end

   if isfile(changelog_in)
      cp(changelog_in, changelog_out; force = true)
   else
      @warn "Missing CHANGELOG.md at repo root: $changelog_in"
   end
   # API page is kept as a checked-in source file. Do not rewrite it during
   # the docs build: an interrupted build can leave docs/src/api.md empty.
end

architecture_page = joinpath(@__DIR__, "src", "architecture_overview.md")
isfile(architecture_page) || error("Missing documentation page: $architecture_page")

pages = [
   "Home" => "index.md",
   "Architecture" => "architecture_overview.md",
   "Theory" => "theorie-eng.md",
   "API" => "api.md",
   "Example" => "minimal_ybus_demo.md",
   "Workshop notebook" => "generated/workshop_apslf.md",
   "Changelog" => "CHANGELOG.md",
]

makedocs(
   modules = [AnalyticLoadFlow],
   sitename = "AnalyticLoadFlow.jl",
   authors = "Udo Schmitz and contributors",
   format = Documenter.HTML(
      prettyurls = get(ENV, "CI", "false") == "true",
      canonical = nothing,
      repolink = nothing,
      mathengine = Documenter.MathJax3(),
   ),
   remotes = nothing,
   checkdocs = :none,
   doctest = :none,
   warnonly = true,
   pages = pages,
)
