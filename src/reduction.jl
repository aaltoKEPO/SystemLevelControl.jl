# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines the functions and auxiliary data-types for the 
##  dimensionality reduction of generalized plant models 

# FUNCTIONS _____________________________________________________________

function sparsity_dim_reduction(P::AbstractGeneralizedPlant, cⱼ::AbstractVector, 𝓢::AbstractVector)
    # Defines the reduced model
    if P isa AbstractGeneralizedPlant{T,StateFeedback} where T
        # sₓ = unique(findnz((𝓢[1][end]*(P.A.≠0))[:,cⱼ])[1])
        # sᵤ = unique(findnz(𝓢[2][end][:,cⱼ])[1])
        sₓ,sᵤ = (unique(findnz((𝓢ⱼ[end]*(P.A.≠0))[:,cⱼ])[1]) for 𝓢ⱼ in 𝓢);   
        P̃ = view(P, (sₓ, [sₓ;P.Nx.+sᵤ]), (sₓ, cⱼ, sᵤ));
    else
        sₓ,sᵤ,sᵧ = (unique(findnz((𝓢ⱼ[end])[:,cⱼ])[1]) for 𝓢ⱼ in 𝓢);   
        P̃ = view(P, (sₓ, [sₓ;P.Nx.+sᵤ], sᵧ), (sₓ, cⱼ, sᵤ));
    end

    # Appropriate identity matrix from original state-space
    iiₓ = indexin(sₓ,cⱼ) .≠ nothing; 
    Ĩ = [I(P̃.Nx)[:,iiₓ] spzeros(P̃.Nx, P̃.Nw-sum(iiₓ))] #(cⱼ[end] ≤ P.Nx)*I;
    # --

    return P̃, Ĩ, iiₓ, sₓ, sᵤ
end

# _______________________________________________________________________
