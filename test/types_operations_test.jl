# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
using SystemLevelControl, Test 
using LinearAlgebra, SparseArrays
# --

## AUXILIARY VARIALBES __________________________________________________
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

## ADJOINT PLANTS _______________________________________________________
## Sparse large-scale systems ___________________________________________
P_large = Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
P_large_adjoint = Plant(A', C₁', C₂', B₁', D₁₁', D₂₁', B₂', D₁₂', D₂₂')

P_large_dual = @inferred DualGeneralizedPlant{Float64,OutputFeedback}(P_large);

@test P_large_dual === P_large'
@test P_large === (P_large')'
@test P_large_adjoint == P_large'

@test P_large' isa DualGeneralizedPlant{Float64,OutputFeedback}

## Sparse large-scale state-feedback systems ____________________________
P_SF_large = Plant(A, B₁, B₂, C₁, D₁₁, D₁₂)
P_SF_large_adjoint = Plant(A', C₁', I(Nx), B₁', D₁₁', spzeros(Nw,Nx), B₂', D₁₂', spzeros(Nu,Nx))

P_SF_large_dual = @inferred DualGeneralizedPlant{Float64,StateFeedback}(P_SF_large);

@test P_SF_large_dual == P_SF_large'
@test P_SF_large === (P_SF_large')'
@test P_SF_large_adjoint == P_SF_large'

@test P_SF_large' isa DualGeneralizedPlant{Float64,StateFeedback}

## SUBPLANTS / VIEWS / SLICING __________________________________________
## Sparse large-scale systems ___________________________________________
II = (1:2000, [1:500; 900:1200], :);
JJ = (1:2000, 1, [1; 3; 6]);

P_large = Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
P_large_slice = Plant(A[II[1],JJ[1]], B₁[II[1],JJ[2]], B₂[II[1],JJ[3]], 
                      C₁[II[2],JJ[1]], D₁₁[II[2],JJ[2]], D₁₂[II[2],JJ[3]], 
                      C₂[II[3],JJ[1]], D₂₁[II[3],JJ[2]], D₂₂[II[3],JJ[3]])

P_large_subsystem = @inferred GeneralizedSubPlant{Float64,OutputFeedback}(P_large, II, JJ)

@test P_large_subsystem === view(P_large, II, JJ)
@test P_large_slice == view(P_large, II, JJ)

@test P_large_slice == P_large[II,JJ]

@test P_large_subsystem !== copy(P_large_subsystem)
@test P_large_subsystem == copy(P_large_subsystem)

@test view(P_large, II, JJ) isa GeneralizedSubPlant{Float64,OutputFeedback}
@test P_large[II,JJ] isa GeneralizedPlant{Float64,OutputFeedback}

## Sparse large-scale state-feedback systems ____________________________
# Explicit reference signal model
II = (1:2000, [1:500; 900:1200]);
JJ = (1:2000, 1, [1; 3; 6]);

P_SF_large = Plant(A, B₁, B₂, C₁, D₁₁, D₁₂)
P_SF_large_slice = Plant(A[II[1],JJ[1]], B₁[II[1],JJ[2]], B₂[II[1],JJ[3]], 
                      C₁[II[2],JJ[1]], D₁₁[II[2],JJ[2]], D₁₂[II[2],JJ[3]])

P_SF_large_subsystem = @inferred GeneralizedSubPlant{Float64,StateFeedback}(P_SF_large, II, JJ)

@test P_SF_large_subsystem === view(P_SF_large, II, JJ)
@test P_SF_large_slice == view(P_SF_large, II, JJ)

@test P_SF_large_slice == P_SF_large[II,JJ]

@test P_SF_large_subsystem !== copy(P_SF_large_subsystem)
@test P_SF_large_subsystem == copy(P_SF_large_subsystem)

@test view(P_SF_large, II, JJ) isa GeneralizedSubPlant{Float64,StateFeedback}
@test P_SF_large[II,JJ] isa GeneralizedPlant{Float64,StateFeedback}

# Standard reference signal model 
II = (1:2000, [1:2000; Nx.+[1;3;6]]);    # The second argument must be explictly given
JJ = (1:2000, 1, [1; 3; 6]);

P_SF_large = Plant(A, B₁, B₂)
P_SF_large_slice = Plant(A[II[1],JJ[1]], B₁[II[1],JJ[2]], B₂[II[1],JJ[3]])

P_SF_large_subsystem = @inferred GeneralizedSubPlant{Float64,StateFeedback}(P_SF_large, II, JJ)

@test P_SF_large_subsystem === view(P_SF_large, II, JJ)
@test P_SF_large_slice == view(P_SF_large, II, JJ)

@test P_SF_large_slice == P_SF_large[II,JJ]

@test P_SF_large_subsystem !== copy(P_SF_large_subsystem)
@test P_SF_large_subsystem == copy(P_SF_large_subsystem)

@test view(P_SF_large, II, JJ) isa GeneralizedSubPlant{Float64,StateFeedback}
@test P_SF_large[II,JJ] isa GeneralizedPlant{Float64,StateFeedback}

## NESTED OPERATIONS __________________________________________
## Sparse large-scale systems _________________________________
II = (1:2000, [1:500; 900:1200], :);
JJ = (1:2000, 1, [1; 3; 6]);

P_large = Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
P_large_adjoint_subsystem = Plant(A'[II[1],JJ[1]], C₁'[II[1],JJ[2]], C₂'[II[1],JJ[3]], 
                                  B₁'[II[2],JJ[1]], D₁₁'[II[2],JJ[2]], D₂₁'[II[2],JJ[3]], 
                                  B₂'[II[3],JJ[1]], D₁₂'[II[3],JJ[2]], D₂₂'[II[3],JJ[3]])

@test P_large_adjoint_subsystem == view(P_large', II, JJ)
@test P_large_adjoint_subsystem == P_large'[II,JJ]

P_large_subsystem_adjoint = Plant(A[II[1],JJ[1]]', C₁[II[2],JJ[1]]', C₂[II[3],JJ[1]]', 
                                  B₁[II[1],JJ[2]]', D₁₁[II[2],JJ[2]]', D₂₁[II[3],JJ[2]]', 
                                  B₂[II[1],JJ[3]]', D₁₂[II[2],JJ[3]]', D₂₂[II[3],JJ[3]]')

@test P_large_subsystem_adjoint == view(P_large, II, JJ)'
@test P_large_subsystem_adjoint == P_large[II,JJ]'
