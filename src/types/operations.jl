# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines operations/functions that can be applied to AbstractGeneralizedPlant's 
##  and their specialized subtypes

# PLANT OPERATIONS __________________________________________________________
function Base.:(==)(P1::AbstractGeneralizedPlant, P2::AbstractGeneralizedPlant)
    if P1 isa DualGeneralizedPlant || P2 isa DualGeneralizedPlant   # This avoids elementwise comparison in adjoint arrays
        return all(norm(getfield(P1,f)-getfield(P2,f)) <= eps() for f in fieldnames(GeneralizedPlant)[1:9])
    else 
        return all(getfield(P1,f)==getfield(P2,f) for f in fieldnames(GeneralizedPlant)[1:9])
    end
end

Base.:size(P::AbstractGeneralizedPlant) = (P.Nx+P.Nz+P.Ny, P.Nx+P.Nu+P.Nw);
Base.:size(P::AbstractGeneralizedPlant, i::Int) = (P.Nx+P.Nz+P.Ny, P.Nx+P.Nu+P.Nw)[i];
Base.:ndims(P::AbstractGeneralizedPlant) = 2;

# Overloads the iterate() interface to unpack the system matrices
Base.:iterate(P::AbstractGeneralizedPlant) = return (P.A, Val(:B₁));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:B₁}) = return (P.B₁, Val(:B₂));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:B₂}) = return (P.B₂, Val(:C₁));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:C₁}) = return (P.C₁, Val(:D₁₁));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:D₁₁}) = return (P.D₁₁, Val(:D₁₂));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:D₁₂}) = return (P.D₁₂, Val(:C₂));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:C₂}) = return (P.C₂, Val(:D₂₁));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:D₂₁}) = return (P.D₂₁, Val(:D₂₂));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:D₂₂}) = return (P.D₂₂, Val(:done));
Base.:iterate(P::AbstractGeneralizedPlant, ::Val{:done}) = return nothing;

# System/Block algebra
Base.:adjoint(P::GeneralizedPlant{T,Ts}) where {T,Ts} = DualGeneralizedPlant{T,Ts}(P)
Base.:adjoint(P::DualGeneralizedPlant{T,Ts}) where {T,Ts} = P.parent

# ___________________________________________________________________________
