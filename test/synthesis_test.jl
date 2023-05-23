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
D₂₁ = 1.2*I(Ny)

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
H2_LLQR = 2π * norm([C₁*Φ[1]*B₁+D₁₂*Φ[2]*B₁ for Φ in zip(Φₓ,Φᵤ)], :𝓗₂);

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

## STATE-FEEDBACK SLS w/ DISTURBANCES ___________________________________
B₁d = [0.5I(Nx÷2+1) zeros(Nx÷2+1,Nx÷2); zeros(Nx÷2,Nx÷2+1) 1.5I(Nx÷2)];
P_chain_dSF = Plant(A, B₁d, B₂, C₁, 0, D₁₂);

Φₓ,Φᵤ = SLS(P_chain_dSF, [𝓢ₓₓ,𝓢ᵤₓ]);

## 1st Test: The closed-loop 𝓗₂-norm of the SLS solutions is approximately
#   that of the (LQR) centralized solution 
Sx,_,_,_ = ared(Matrix(A), Matrix(B₂), D₁₂'D₁₂, C₁'C₁, C₁'D₁₂);

H2_CLQR = tr(B₁d'Sx*B₁d);
H2_LLQR = 2π * norm([C₁*Φ[1]*B₁d + D₁₂*Φ[2]*B₁d for Φ in zip(Φₓ,Φᵤ)], :𝓗₂);

@test (H2_LLQR / H2_CLQR) < 1.18

# -----------------------------------------------------------------------

# OUTPUT-FEEDBACK SLS __________________________________________________
BB₁ = [B₁ zeros(Nx,Ny)]
DD₂₁ = [zeros(Ny,Nx) D₂₁]

P_chain_OF = Plant(A, BB₁, B₂, 
                   C₁, 0, D₁₂, 
                   C₂, DD₂₁, 0);

Φₓₓ,Φᵤₓ,Φₓᵧ,Φᵤᵧ = SLS(P_chain_OF, [𝓢ₓₓ,𝓢ᵤₓ,𝓢ₓᵧ,𝓢ᵤᵧ]);

## 1st Test: The closed-loop 𝓗₂-norm of the SLS solutions is approximately
#   that of the (LQR) centralized solution 
Sx,_,_,_ = ared(Matrix(A), Matrix(B₂), D₁₂'D₁₂, C₁'C₁, C₁'D₁₂);
Sy,_,_,_ = ared(Matrix(A)', Matrix(C₂)', DD₂₁*DD₂₁', BB₁*BB₁', BB₁*DD₂₁');

Rb = (I + B₂'Sx*B₂);
F₂ = -Rb\B₂'Sx*A;
F₀ = -Rb\B₂'Sx*BB₁;
L₀ = (F₂*Sy*C₂' + F₀*DD₂₁')/(I + C₂*Sy*C₂');

H2_CLQG = tr(F₀'D₁₂'D₁₂*F₀ + (BB₁+B₂*F₀)'Sx*(BB₁+B₂*F₀)) + tr(Rb*((L₀*DD₂₁ - F₀)*(L₀*DD₂₁ - F₀)' + (L₀*C₂ - F₂)*Sy*(L₀*C₂ - F₂)'));
H2_LLQG = 2π * norm([(C₁*Φ[1]+D₁₂*Φ[2])*BB₁ + (C₁*Φ[3]+D₁₂*Φ[4])*DD₂₁ for Φ in zip(Φ̃ₓₓ,Φ̃ᵤₓ,Φ̃ₓᵧ,Φ̃ᵤᵧ)], :𝓗₂);

@test (H2_LLQG / H2_CLQG) < 2

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


# !!! The code below solves the Output-SLS problem directly without
# !!!  using partial separability. It is used only for testing purposes
using JuMP, SCS 

BB₁ = [B₁ zeros(Nx,Ny)]
DD₂₁ = [zeros(Ny,Nx) D₂₁]

problem = Model(SCS.Optimizer);
Φₓₓ = [@variable(problem, [1:Nx,1:Nx]) for _ in 1:T];
Φᵤₓ = [@variable(problem, [1:Nu,1:Nx]) for _ in 1:T];
Φₓᵧ = [@variable(problem, [1:Nx,1:Ny]) for _ in 1:T];
Φᵤᵧ = [@variable(problem, [1:Nu,1:Ny]) for _ in 1:T];

T_zw = [@expression(problem, [C₁ D₁₂]*[Φ[1] Φ[3]; Φ[2] Φ[4]]*[BB₁; DD₂₁]) for Φ in zip(Φₓₓ,Φᵤₓ,Φₓᵧ,Φᵤᵧ)];

@objective(problem,      Min,       norm(T_zw, :𝓗₂));

@constraint(problem,                Φₓₓ[1]   .== I);
@constraint(problem, [t = 1:(T-1)], Φₓₓ[t+1] .== A*Φₓₓ[t] + B₂*Φᵤₓ[t]);
@constraint(problem,                    0    .== A*Φₓₓ[T] + B₂*Φᵤₓ[T]);

@constraint(problem,                Φₓᵧ[1]   .== 0);
@constraint(problem, [t = 1:(T-1)], Φₓᵧ[t+1] .== A*Φₓᵧ[t] + B₂*Φᵤᵧ[t]);
@constraint(problem,                    0    .== A*Φₓᵧ[T] + B₂*Φᵤᵧ[T]);

@constraint(problem, [t = 1:(T-1)], Φₓₓ[t+1] .== Φₓₓ[t]*A + Φₓᵧ[t]*C₂);
@constraint(problem,                    0    .== Φₓₓ[T]*A + Φₓᵧ[T]*C₂);

@constraint(problem,                Φᵤₓ[1]   .== 0);
@constraint(problem, [t = 1:(T-1)], Φᵤₓ[t+1] .== Φᵤₓ[t]*A + Φᵤᵧ[t]*C₂);
@constraint(problem,                    0    .== Φᵤₓ[T]*A + Φᵤᵧ[T]*C₂);

for t in 1:T
    fix.(Φₓₓ[t][𝓢ₓₓ[t] .≠ 1], 0.0, force=true);
    fix.(Φᵤₓ[t][𝓢ᵤₓ[t] .≠ 1], 0.0, force=true);
    fix.(Φₓᵧ[t][𝓢ₓᵧ[t] .≠ 1], 0.0, force=true);
    fix.(Φᵤᵧ[t][𝓢ᵤᵧ[t] .≠ 1], 0.0, force=true);
end

optimize!(problem)

Φ̃ₓₓ = [value.(Φₓₓ[t]) .* 𝓢ₓₓ[t] for t in 1:T];
Φ̃ᵤₓ = [value.(Φᵤₓ[t]) .* 𝓢ᵤₓ[t] for t in 1:T];
Φ̃ₓᵧ = [value.(Φₓᵧ[t]) .* 𝓢ₓᵧ[t] for t in 1:T];
Φ̃ᵤᵧ = [value.(Φᵤᵧ[t]) .* 𝓢ᵤᵧ[t] for t in 1:T];
# -- --