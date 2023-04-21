# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines the functions and auxiliary data-types for the 
##  dimensionality reduction of generalized plant models 

# FUNCTIONS _____________________________________________________________

function sparsity_dim_reduction(P::AbstractGeneralizedPlant, cⱼ::AbstractVector, 𝓢::AbstractVector)
    # Defines the reduced model
    if P isa AbstractGeneralizedPlant{T,StateFeedback} where T
        sₓ,sᵤ = (unique(findnz((𝓢ⱼ[end]*(P.A.≠0))[:,cⱼ])[1]) for 𝓢ⱼ in 𝓢);   
        P̃ = Plant(P.A[sₓ,sₓ],              P.B₁[sₓ,cⱼ],            P.B₂[sₓ,sᵤ],
                  P.C₁[[sₓ;P.Nx.+sᵤ], sₓ], P.D₁₁[[sₓ;P.Nx.+sᵤ],cⱼ], P.D₁₂[[sₓ;P.Nx.+sᵤ],sᵤ])
    else
        sₓ,sᵤ,sᵧ = (unique(findnz((𝓢ⱼ[end])[:,cⱼ])[1]) for 𝓢ⱼ in 𝓢);   
        P̃ = Plant(P.A[sₓ,sₓ],              P.B₁[sₓ,cⱼ],            P.B₂[sₓ,sᵤ],
                  P.C₁[[sₓ;P.Nx.+sᵤ], sₓ], P.D₁₁[ [sₓ;P.Nx.+sᵤ],cⱼ], P.D₁₂[[sₓ;P.Nx.+sᵤ],sᵤ],
                  P.C₂[sᵧ,sₓ],              P.B₁[sᵧ,cⱼ],            P.B₂[sᵧ,sᵤ])
    end

    # Appropriate identity matrix from original state-space
    iiₓ = indexin(sₓ,cⱼ) .≠ nothing; 
    Ĩ = [I(P̃.Nx)[:,iiₓ] spzeros(P̃.Nx, P̃.Nw-sum(iiₓ))] #(cⱼ[end] ≤ P.Nx)*I;
    # --

    return P̃, Ĩ, iiₓ, sₓ, sᵤ
end

# _______________________________________________________________________
