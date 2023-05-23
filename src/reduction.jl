# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines the functions and auxiliary data-types for the 
##  dimensionality reduction of generalized plant models 

# FUNCTIONS _____________________________________________________________

function sparsity_dim_reduction(P::AbstractGeneralizedPlant, cⱼ::AbstractVector, Sₓ::AbstractVector, Sᵤ::AbstractVector)
    # Obtains the (sorted) indexes of the (d+1)-Localized neighborhood of x(cⱼ)
    sₓ = unique(findnz(((P.A .≠ 0) * Sₓ[end])[:,cⱼ])[1]);   
    sᵤ = unique(findnz(Sᵤ[end][:,cⱼ])[1]);   
    
    # Creates the reduced-order model
    P̃ = view(P, (sₓ, [sₓ; P.Nx.+sᵤ], :), (sₓ, cⱼ, sᵤ));

    # Appropriate identity matrix from original state-space
    iᵣ = [indexin(sₓ, cⱼ); indexin(1:P̃.Ny, cⱼ .- P.Nx)] .≠ nothing;
    Ĩ = I(P̃.Nx+P̃.Ny)[1:P̃.Nx,iᵣ];
    
    # --
    return P̃, Ĩ, iᵣ, sₓ, sᵤ
end

# _______________________________________________________________________
