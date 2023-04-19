# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines operations/functions that can be applied to AbstractGeneralizedPlant's 
##  and their specialized subtypes

# PLANT OPERATIONS __________________________________________________________
Base.:(==)(P1::GeneralizedPlant, P2::GeneralizedPlant) = all(getfield(P1,f)==getfield(P2,f) for f in fieldnames(GeneralizedPlant))
Base.:size(P::AbstractGeneralizedPlant) = (P.Nx+P.Nz+P.Ny, P.Nu+P.Nw);
Base.:size(P::AbstractGeneralizedPlant, i::Int) = (P.Nx+P.Nz+P.Ny, P.Nu+P.Nw)[i];
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
Base.:adjoint(P::AbstractGeneralizedPlant{T,OutputFeedback}) where {T} = Plant(P.A', P.C₁', P.C₂', P.B₁', P.D₁₁', P.D₂₁', P.B₂', P.D₁₂', P.D₂₂')
Base.:adjoint(P::AbstractGeneralizedPlant{T,StateFeedback}) where {T} = Plant(P.A', P.C₁', P.C₂', P.B₁', P.D₁₁', spzeros(size(P.B₁))', P.B₂', P.D₁₂', spzeros(size(P.B₂))')

# ___________________________________________________________________________
