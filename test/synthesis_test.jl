# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
using SystemLevelControl, Test 
using LinearAlgebra, SparseArrays, MatrixEquations
# --

## STATE-FEEDBACK SLS ___________________________________________________
## AUXILIARY VARIABLES (59-Chain System)
Nx, Nu = (59, 20);
d,T,α = (9, 29, 1.5);
    
A = I + spdiagm(1 => 0.2ones(Nx-1)) - spdiagm(-1 => 0.2ones(Nx-1));
B₁ = I(Nx);
B₂ = spdiagm(0 => ones(Nx))[:,vec((1:2).+6(0:9)')];

P = Plant(A, B₁, B₂);
𝓢ₓ = [           (A .≠ 0)^min(d,  floor(α*(t-1))) .≠ 0 for t = 1:T];
𝓢ᵤ = [(B₂' .≠ 0)*(A .≠ 0)^min(d+1,floor(α*(t-1))) .≠ 0 for t = 1:T];

Φₓ,Φᵤ = SLS_𝓗₂(P, [𝓢ₓ,𝓢ᵤ]);

## 1st Test: The closed-loop 𝓗₂-norm of the SLS solutions is approximately
#   that of the (LQR) centralized solution 
S,_,_,_ = ared(Matrix(A), Matrix(B₂), 1.0I(Nu), 1.0I(Nx), 0I(Nx)[:,1:Nu]);

H2_CLQR = tr(S)/2π;
H2_LLQR = norm([P.C₁*Φ[1]+P.D₁₂*Φ[2] for Φ in zip(Φₓ,Φᵤ)], :𝓗₂);

@test (H2_LLQR / H2_CLQR) < 1.18 

## 2nd Test: The system is controllable, thus it should be (d,T)-localized
# Simulates the closed-loop system 
w(t) = (t==1)*I(59)[:,30]
x = spzeros(Nx,2T); 
u = spzeros(Nu,2T)
β = similar(x); 
for t = 1:(2T-1)
    β[:,t+1] = sum([Φₓ[τ+1]*(x[:,t+1-τ] - β[:,t+1-τ]) for τ = 1:min(t,T-1)]);
    u[:,t]   = sum([Φᵤ[ τ ]*(x[:,t+1-τ] - β[:,t+1-τ]) for τ = 1:min(t,T)  ]);
    
    x[:,t+1] = A*x[:,t] + B₁*w(t) + B₂*u[:,t];
end

@test norm(x[:,T+2])^2 <= 1e-10     # Tests T-localization
@test norm(x[[1:29-d; 31+d:Nx],:])^2 <= 1e-10   # Tests d-Localization
@test all([norm(x[[29-t; 31+t],t+1])^2 <= 1e-10 for t in 1:T-1])    # Tests α-Localization

## 3rd Test: Check if dimensionality reduction is working as intended,
#   i.e., if the parallel and complete optimizations are equivalent
Φₓ_2,Φᵤ_2 = SLS_𝓗₂(P, [𝓢ₓ,𝓢ᵤ], 𝓘=[(1+2j):min(2+2j, Nx) for j in 0:(Nx÷2)]);
Φₓ_T,Φᵤ_T = SLS_𝓗₂(P, [𝓢ₓ,𝓢ᵤ], 𝓘=[1:Nx]);

for t = 1:T
    @test norm(Φₓ[t] - Φₓ_2[t])^2 < 1.5e-4
    @test norm(Φₓ[t] - Φₓ_T[t])^2 < 1.5e-4
    
    @test norm(Φᵤ[t] - Φᵤ_2[t])^2 < 1.5e-4
    @test norm(Φᵤ[t] - Φᵤ_T[t])^2 < 1.5e-4
end

# -----------------------------------------------------------------------
