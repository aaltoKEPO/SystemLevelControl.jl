# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
using SystemLevelControl, Test 
using LinearAlgebra, SparseArrays, MatrixEquations
# --

## AUXILIARY VARIABLES (59-Chain System)
# System and controller definitions
Nx, Nu, Ny = (59, 20, 20);
d,T,α = (9, 29, 1.5);
    
A = I + spdiagm(1 => 0.2ones(Nx-1)) - spdiagm(-1 => 0.2ones(Nx-1));
B₁ = I(Nx);
B₂ = spdiagm(0 => ones(Nx))[:,vec((1:2).+6(0:9)')];

C₁ = [I(Nx); 0I(Nx)[1:Nu,:]];
D₁₂ = [0I(Nx)[:,1:Nu]; I(Nu)];

C₂ = spdiagm(0 => ones(Nx))[vec((1:2).+6(0:9)'),:];
D₂₁ = 1.2*C₂*I(Nx)

# Sparisty constraints
A_sp, B_sp, C_sp = (A .≠ 0, B₂ .≠ 0, C₂ .≠ 0)

𝓢ₓₓ = [     A_sp^min(d,  floor(α*(t-1)))       .≠ 0 for t = 1:T];
𝓢ᵤₓ = [B_sp'A_sp^min(d+1,floor(α*(t-1)))       .≠ 0 for t = 1:T];
𝓢ₓᵧ = [     A_sp^min(d,  floor(α*(t-1)))*C_sp' .≠ 0 for t = 1:T];
𝓢ᵤᵧ = [B_sp'A_sp^min(d+1,floor(α*(t-1)))*C_sp' .≠ 0 for t = 1:T];


## STATE-FEEDBACK SLS ___________________________________________________
P_chain_SF = Plant(A, B₁, B₂);

Φₓ,Φᵤ = SLS(P_chain_SF, [𝓢ₓₓ,𝓢ᵤₓ]);

## 1st Test: The closed-loop 𝓗₂-norm of the SLS solutions is approximately
#   that of the (LQR) centralized solution 
Sx,_,_,_ = ared(Matrix(A), Matrix(B₂), D₁₂'D₁₂, C₁'C₁, C₁'D₁₂);

H2_CLQR = tr(B₁'Sx*B₁);
H2_LLQR = 2π * norm([C₁*Φ[1]+D₁₂*Φ[2] for Φ in zip(Φₓ,Φᵤ)], :𝓗₂);

@test (H2_LLQR / H2_CLQR) < 1.18 

## 2nd Test: The system is controllable, thus it should be (d,T)-localized
# Simulates the closed-loop system 
x = spzeros(Nx,2T); 
u = spzeros(Nu,2T)
β = similar(x); 

w(t) = (t==1)*I(Nx)[:,1+(Nx÷2)]

for t = 1:(2T-1)
    β[:,t+1] = sum([Φₓ[τ+1]*(x[:,t+1-τ] - β[:,t+1-τ]) for τ = 1:min(t,T-1)]);
    u[:,t]   = sum([Φᵤ[ τ ]*(x[:,t+1-τ] - β[:,t+1-τ]) for τ = 1:min(t,T)  ]);
    
    x[:,t+1] = A*x[:,t] + B₁*w(t) + B₂*u[:,t];
end

@test norm(x[:,T+2])^2 <= eps()     # Tests T-localization
@test norm(x[[1:29-d; 31+d:Nx],:])^2 <= eps()   # Tests d-Localization
@test all([norm(x[[29-t; 31+t],t+1])^2 <= eps() for t in 1:T-1])    # Tests α-Localization

## 3rd Test: Check if dimensionality reduction is working as intended,
#   i.e., if the parallel and complete optimizations are equivalent
Φₓ_2,Φᵤ_2 = SLS(P_chain_SF, [𝓢ₓₓ,𝓢ᵤₓ], J=[(1+2j):min(2+2j, Nx) for j in 0:(Nx÷2)]);
Φₓ_T,Φᵤ_T = SLS(P_chain_SF, [𝓢ₓₓ,𝓢ᵤₓ], J=[1:Nx]);

@test all([norm(Φₓ[t] - Φₓ_2[t])^2 < 1e-8  for t in 1:T])
@test all([norm(Φᵤ[t] - Φᵤ_2[t])^2 < 1e-8  for t in 1:T])

@test all([norm(Φₓ[t] - Φₓ_T[t])^2 < 1e-8  for t in 1:T])
@test all([norm(Φᵤ[t] - Φᵤ_T[t])^2 < 1e-8  for t in 1:T])

# -----------------------------------------------------------------------

# ## STATE-FEEDBACK SLS w/ DISTURBANCES ___________________________________
# B₁d = [0.5I(Nx÷2+1) zeros(Nx÷2+1,Nx÷2); zeros(Nx÷2,Nx÷2+1) 1.5I(Nx÷2)];
# P_chain_dSF = Plant(A, B₁d, B₂, C₁, 0, D₁₂);

# Φₓ,Φᵤ = SLS(P_chain_dSF, [𝓢ₓₓ,𝓢ᵤₓ]);

# ## 1st Test: The closed-loop 𝓗₂-norm of the SLS solutions is approximately
# #   that of the (LQR) centralized solution 
# Sx,_,_,_ = ared(Matrix(A), Matrix(B₂), 1.0I(Nu), 1.0I(Nx), 0I(Nx)[:,1:Nu]);
# Sy,_,_,_ = ared(Matrix(A)', Matrix(C₂)', 0I(Ny), B₁d*B₁d', 0I(Nx)[:,1:Ny]);
# K = (I + B₂'Sx*B₂)\B₂'Sx*A;

# H2_CLQR = tr(B₁'S*B₁) + tr(D₁₂'D₁₂*K*Sy*K');
# H2_LLQR = 2π * norm([C₁*Φ[1]+D₁₂*Φ[2] for Φ in zip(Φₓ,Φᵤ)], :𝓗₂);

# @test (H2_LLQR / H2_CLQR) < 1.18

# -----------------------------------------------------------------------

## OUTPUT-FEEDBACK SLS __________________________________________________
# B₁ = [0.5I(Nx÷2+1) zeros(Nx÷2+1,Nx÷2); zeros(Nx÷2,Nx÷2+1) 1.5I(Nx÷2)];  # Sligthly more fun disturbance matrix

# P_chain_OF = Plant(A, B₁, B₂, 
#                    C₁, 0, D₁₂, 
#                    C₂, D₂₁, 0);

# Φₓₓ,Φᵤₓ,Φₓᵧ,Φᵤᵧ = SLS(P_chain_OF, [𝓢ₓₓ,𝓢ᵤₓ,𝓢ₓᵧ,𝓢ᵤᵧ]);

# ## 1st Test: The closed-loop 𝓗₂-norm of the SLS solutions is approximately
# #   that of the (LQR) centralized solution 
# Sx,_,_,_ = ared(Matrix(A), Matrix(B₂), D₁₂'D₁₂, C₁'C₁, C₁'D₁₂);
# Sy,_,_,_ = ared(Matrix(A)', Matrix(C₂)', D₂₁'D₂₁, B₁*B₁', B₁'D₂₁);
# K = (I + B₂'Sx*B₂)\B₂'Sx*A;

# H2_CLQG = tr(B₁'Sx*B₁) + tr(D₁₂'D₁₂*K*Sy*K');
# H2_LLQG = 2π * norm([(C₁*Φ[1]+D₁₂*Φ[2])*B₁ + (C₁*Φ[3]+D₁₂*Φ[4])*D₂₁ for Φ in zip(Φₓₓ,Φᵤₓ,Φₓᵧ,Φᵤᵧ)], :𝓗₂);

# @test (H2_LLQG / H2_CLQG) < 1.18

# ## 2nd Test: The system is controllable and observable, thus it should be (d,T)-localized
# # Simulates the closed-loop system 
# x = spzeros(Nx,2T); 
# u = spzeros(Nu,2T)
# y = spzeros(Ny,2T); 
# β = similar(x); 

# w(t) = (t==1)*I(Nx)[:,1+(Nx÷2)]

# for t = 1:(2T-1)
#     β[:,t+1] = sum([Φₓ[τ+1]*(x[:,t+1-τ] - β[:,t+1-τ]) for τ = 1:min(t,T-1)]);
#     u[:,t]   = sum([Φᵤ[ τ ]*(x[:,t+1-τ] - β[:,t+1-τ]) for τ = 1:min(t,T)  ]);
    
#     x[:,t+1] = A*x[:,t] + B₁*w(t) + B₂*u[:,t];
#     y[:,t] = C₂*x[:,t] + D₂₁*w(t);
# end

# @test norm(x[:,T+2])^2 <= eps()     # Tests T-localization
# @test norm(x[[1:29-d; 31+d:Nx],:])^2 <= eps()   # Tests d-Localization
# @test all([norm(x[[29-t; 31+t],t+1])^2 <= eps() for t in 1:T-1])    # Tests α-Localization

# ## 3rd Test: Check if dimensionality reduction is working as intended,
# #   i.e., if the parallel and complete optimizations are equivalent
# Φₓₓ_2,Φᵤₓ_2,Φₓᵧ_2,Φᵤᵧ_2 = SLS(P_chain_OF, [𝓢ₓₓ,𝓢ᵤₓ,𝓢ₓᵧ,𝓢ᵤᵧ], J=[(1+2j):min(2+2j, Nx) for j in 0:(Nx÷2)]);
# Φₓₓ_T,Φᵤₓ_T,Φₓᵧ_T,Φᵤᵧ_T = SLS(P_chain_OF, [𝓢ₓₓ,𝓢ᵤₓ,𝓢ₓᵧ,𝓢ᵤᵧ], J=[1:Nx]);

# @test all([norm(Φₓₓ[t] - Φₓₓ_2[t])^2 < 1e-8  for t in 1:T])
# @test all([norm(Φᵤₓ[t] - Φᵤₓ_2[t])^2 < 1e-8  for t in 1:T])
# @test all([norm(Φₓᵧ[t] - Φₓᵧ_2[t])^2 < 1e-8  for t in 1:T])
# @test all([norm(Φᵤᵧ[t] - Φᵤᵧ_2[t])^2 < 1e-8  for t in 1:T])

# @test all([norm(Φₓₓ[t] - Φₓₓ_T[t])^2 < 1e-8  for t in 1:T])
# @test all([norm(Φᵤₓ[t] - Φᵤₓ_T[t])^2 < 1e-8  for t in 1:T])
# @test all([norm(Φₓᵧ[t] - Φₓᵧ_T[t])^2 < 1e-8  for t in 1:T])
# @test all([norm(Φᵤᵧ[t] - Φᵤᵧ_T[t])^2 < 1e-8  for t in 1:T])

# -----------------------------------------------------------------------
