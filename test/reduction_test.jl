# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
using SystemLevelControl, Test 
using LinearAlgebra, SparseArrays
# --

# Test sparsity-based dimensionality reduction
Nx = 59
A = I + spdiagm(1 => 0.2ones(Nx-1)) - spdiagm(-1 => 0.2ones(Nx-1));
B₁ = I(Nx);
B₂ = spdiagm(0 => ones(Nx))[:,vec((1:2).+6(0:9)')];

P_full = Plant(A, B₁, B₂)

Sₓ = [           (A .≠ 0)^9 .≠ 0];
Sᵤ = [(B₂' .≠ 0)*(A .≠ 0)^9 .≠ 0];

cⱼ,Ĩ,iiₓ,sₓ,sᵤ = (1:20, I(30)[:,1:20], [ones(Int,20);zeros(Int,10)], 1:30, 1:10);
P_redu = Plant(A[sₓ,sₓ], B₁[sₓ,cⱼ], B₂[sₓ,sᵤ])

@test (P_redu,Ĩ,iiₓ,sₓ,sᵤ) == sparsity_dim_reduction(P_full, cⱼ, [Sₓ,Sᵤ])