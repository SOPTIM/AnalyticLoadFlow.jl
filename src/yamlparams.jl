# =============================================================================
# File: src/yamlparams.jl
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

# Purpose:
#   Provides lightweight parameter-loading helpers for YAML-like configuration inputs used by examples or experiments.

# Small, dependency-free YAML subset used by example programs.
# The repository examples intentionally avoid adding a full YAML dependency for
# simple benchmark configuration files. Supported syntax is limited to nested
# dictionaries (2-space indentation) and scalar/list values.

"""
    parse_yaml_scalar(raw::AbstractString)

Parse a small YAML scalar subset into Julia values.

Supported values: booleans, `null`, integers, floats, symbols (`:name`), quoted
or unquoted strings, and one-line lists of scalar values like `[1, 2, auto]`.
"""
function parse_yaml_scalar(raw::AbstractString)
   s = strip(raw)
   isempty(s) && return ""

   if startswith(s, "[") && endswith(s, "]")
      inner = strip(s[2:(end - 1)])
      isempty(inner) && return Any[]
      return Any[parse_yaml_scalar(tok) for tok in split(inner, ",")]
   end

   if (startswith(s, "\"") && endswith(s, "\"")) || (startswith(s, "'") && endswith(s, "'"))
      return s[2:(end - 1)]
   end

   endswith(s, ",") && (s = strip(s[1:(end - 1)]))
   lo = lowercase(s)
   lo in ("true", "yes", "on") && return true
   lo in ("false", "no", "off") && return false
   lo in ("null", "~") && return nothing
   startswith(s, ":") && length(s) > 1 && return Symbol(s[2:end])

   i = tryparse(Int, s)
   i !== nothing && return i
   f = tryparse(Float64, s)
   f !== nothing && return f
   return s
end

"""
    load_yaml_dict(path::AbstractString)::Dict{String,Any}

Load the small YAML subset used by APSLF example configurations.

This is intentionally not a general-purpose YAML parser. It supports nested
maps using two-space indentation and scalar/list values parsed by
[`parse_yaml_scalar`](@ref).
"""
function load_yaml_dict(path::AbstractString)::Dict{String,Any}
   root = Dict{String,Any}()
   dict_stack = [root]

   for raw in readlines(path)
      stripped = strip(raw)
      startswith(stripped, "#") && continue
      line = replace(raw, r"\s+#.*$" => "")
      isempty(strip(line)) && continue

      indent = length(line) - length(lstrip(line))
      (indent % 2 == 0) || error("Invalid YAML indentation in $path: '$raw'")
      level = (indent ÷ 2) + 1
      content = strip(line)

      while length(dict_stack) > level
         pop!(dict_stack)
      end

      m = match(r"^([^:]+):(.*)$", content)
      m === nothing && error("Unsupported YAML line in $path: '$raw'")
      key = strip(m.captures[1])
      value_raw = strip(m.captures[2])
      parent = dict_stack[end]

      if isempty(value_raw)
         child = Dict{String,Any}()
         parent[key] = child
         push!(dict_stack, child)
      else
         parent[key] = parse_yaml_scalar(value_raw)
      end
   end

   return root
end

"""
    merge_yaml_dict!(dst::Dict{String,Any}, src::Dict{String,Any})

Recursively merge `src` into `dst`, preserving nested dictionaries where both
sides are dictionaries. Returns `dst`.
"""
function merge_yaml_dict!(dst::Dict{String,Any}, src::Dict{String,Any})
   for (k, v) in src
      if haskey(dst, k) && dst[k] isa Dict{String,Any} && v isa Dict{String,Any}
         merge_yaml_dict!(dst[k], v)
      else
         dst[k] = v
      end
   end
   return dst
end

"""
    load_yaml_layers(paths)::Dict{String,Any}

Load and merge multiple YAML files in order.

Later files override earlier values. Non-existing paths are ignored so callers
can pass optional default locations.
"""
function load_yaml_layers(paths)::Dict{String,Any}
   cfg = Dict{String,Any}()
   for path in paths
      p = String(path)
      if isfile(p)
         merge_yaml_dict!(cfg, load_yaml_dict(p))
      end
   end
   return cfg
end

"""
    as_bool(x)::Bool

Normalize a bool-like YAML value to `Bool`.
"""
function as_bool(x)::Bool
   x isa Bool && return x
   x isa Integer && return x != 0
   x isa AbstractString && return lowercase(strip(x)) in ("1", "true", "yes", "y", "on")
   error("Expected boolean-compatible value, got $(typeof(x))")
end

"""
    as_int_vector(x)::Vector{Int}

Normalize either a scalar integer-like value or a YAML list into `Vector{Int}`.
"""
function as_int_vector(x)::Vector{Int}
   if x isa AbstractVector
      return Int[Int(v) for v in x]
   end
   return Int[Int(x)]
end
