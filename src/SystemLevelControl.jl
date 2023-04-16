# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
module SystemLevelControl
# --

export AbstractGeneralizedPlant, StateFeedbackPlant, GeneralizedPlant, Plant
export SLS_𝓗₂
export sparsity_dim_reduction
export generateTree

# Sub-modules (in separate files)
using Distributed, SharedArrays
using LinearAlgebra, Statistics, Distributions
using JuMP, Ipopt
using SparseArrays

include("types/GeneralizedPlant.jl") 
include("types/StateFeedbackPlant.jl") 

include("types/operations.jl") 
include("synthesis.jl")
include("reduction.jl")
include("utils.jl")

# --
end
