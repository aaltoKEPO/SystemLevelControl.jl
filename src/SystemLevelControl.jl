# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
module SystemLevelControl
# --

export # Plant and controller types
       AbstractFeedbackStructure, 
       StateFeedback, 
       OutputFeedback,
       AbstractGeneralizedPlant, 
       GeneralizedPlant, 
       DualGeneralizedPlant, 
       GeneralizedSubPlant, 
       Plant,
       # Synthesis methods
       SLS_𝓗₂,
       # Dimensionality reduction methods
       sparsity_dim_reduction,
       # Utils and auxiliary functions
       generateTree

# Sub-modules (in separate files)
using Distributed, SharedArrays
using LinearAlgebra, Statistics, Distributions
using JuMP, Ipopt
using SparseArrays

include("types/FeedbackStructures.jl") 
include("types/GeneralizedPlant.jl") 
include("types/conversions.jl") 

include("types/operations.jl") 
include("synthesis.jl")
include("reduction.jl")
include("utils.jl")

# --
end
