# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
using SystemLevelControl, Test 
using LinearAlgebra, SparseArrays
# --

## Sparse large-scale systems ___________________________________________
Nx,Nu,Nw,Ny = (100000, 52000, 51000, 50000)
A = sprandn(Nx, Nx, 1/Nx);
B₁ = sprandn(Nx, Nw, 1/Nw);
B₂ = sprandn(Nx, Nu, 1/Nu);

C₁ = sprandn(Nx+Nu, Nx, 1/Nx);
D₁₁ = sprandn(Nx+Nu, Nw, 1/Nw);
D₁₂ = sprandn(Nx+Nu, Nu, 1/Nu);

C₂ = sprandn(Ny, Nx, 1/Nx);
D₂₁ = sprandn(Ny, Nw, 1/Nw);
D₂₂ = sprandn(Ny, Nu, 1/Nu);

P_large = Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
P_large_adjoint = Plant(A', C₁', C₂', B₁', D₁₁', D₂₁', B₂', D₁₂', D₂₂')

P_large_dual = DualGeneralizedPlant{Float64,OutputFeedback}(P_large);

@test P_large === (P_large')'
@test P_large_adjoint == P_large'
@test P_large_dual == P_large'
@test P_large_dual == P_large_adjoint
