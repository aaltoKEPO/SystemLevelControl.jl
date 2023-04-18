# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
using SystemLevelControl, Test 
using LinearAlgebra, SparseArrays
# --

## Scalar systems _______________________________________________________
a,b₁,b₂,c₁,d₁₁,d₁₂,c₂,d₂₁,d₂₂ = [sparse([i]') for i in 1.0:9.0]
P_scalar = @inferred GeneralizedPlant{Float64,OutputFeedback}(a, b₁, b₂, c₁, d₁₁, d₁₂, c₂, d₂₁, d₂₂)

@test typeof(P_scalar) <: GeneralizedPlant{Float64,OutputFeedback}

# Test user-friendly constructor and promotion/conversion rules
@test P_scalar === Plant(a,b₁,b₂,c₁,d₁₁,d₁₂,c₂,d₂₁,d₂₂)
@test P_scalar !== Plant(1, 2, 3, 4, 5, 6, 7, 8, 9)
@test P_scalar == Plant(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0)
@test P_scalar == Plant(1, 2, 3, 4, 5, 6, 7, 8, 9)

@test typeof(Plant(1, 2, 3, 4, 5, 6, 7, 8, 9)) <: GeneralizedPlant{Int,OutputFeedback}
@test typeof(Plant(1, 2, 3, 4, 5.0, 6, 7, 8, 9)) <: GeneralizedPlant{Float64,OutputFeedback}
@test typeof(Plant(1, [2], 3, 4, 5.0, 6, 7, 8, 9)) <: GeneralizedPlant{Float64,OutputFeedback}

@test P_scalar == Plant([1 2 3; 4 5 6; 7 8 9], [1, 1, 1, 1, 1])

# Test if the field->values are all correct
@test all(getfield.([P_scalar], [:Nx, :Nw, :Nu, :Nz, :Ny]) .== 1) 

for (i,f) in enumerate([:A, :B₁, :B₂, :C₁, :D₁₁, :D₁₂, :C₂, :D₂₁, :D₂₂])
    @test getfield(P_scalar, f) == sparse([i]')
end

## Vector + scalar systems ______________________________________________
A = sparse([1.0 2.0; 3.0 4.0]);
B₁ = sparse(reshape([5.0; 6.0],:,1));
B₂ = sparse(reshape([7.0; 8.0],:,1));
C₁ = C₂ = sparse([10.0 11.0]);

P_vector = @inferred GeneralizedPlant{Float64,OutputFeedback}(A, B₁, B₂, C₁, d₁₁, d₁₂, C₂, d₂₁, d₂₂)

@test typeof(P_vector) <: GeneralizedPlant{Float64,OutputFeedback}

# Test user-friendly constructor and promotion/conversion rules
@test P_vector === Plant(A, B₁, B₂, C₁, d₁₁, d₁₂, C₂, d₂₁, d₂₂)
@test P_vector == Plant(A, B₁, B₂, [10 11], 5, 6, [10 11], 8, 9)

@test typeof(Plant(A, B₁, B₂, [9 10], 11.0, 12, [13 14], 15, 16)) <: GeneralizedPlant{Float64,OutputFeedback}
@test typeof(Plant(A, B₁, B₂, [9.0 10], 11, 12, [13 14], 15, 16)) <: GeneralizedPlant{Float64,OutputFeedback}

@test P_vector == Plant([A B₁ B₂; [10 11] 5 6; [10 11] 8 9], [2, 1, 1, 1, 1])

# Test if the field->values are all correct, and if the matrices are being passed by reference
@test P_vector.Nx == 2
@test all(getfield.([P_vector], [:Nw, :Nu, :Nz, :Ny]) .== 1) 

@test P_vector.A === A
@test P_vector.B₁ === B₁
@test P_vector.B₂ === B₂
@test P_vector.C₁ === C₁
@test P_vector.C₂ === C₂

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

P_large = @inferred GeneralizedPlant{Float64,OutputFeedback}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)

@test typeof(P_large) <: GeneralizedPlant{Float64,OutputFeedback}

# Test user-friendly constructor and promotion/conversion rules
@test P_large === Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
@test P_large == Plant([A B₁ B₂; C₁ D₁₁ D₁₂; C₂ D₂₁ D₂₂], [Nx, Nx+Nu, Ny, Nw, Nu])

P_large_D0 = Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, 0D₂₂)
@test P_large_D0 == Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, 0)

## State Feedback plants ________________________________________________
C₂ = SparseMatrixCSC{Float64,Int}(I, Nx, Nx)
D₂₁ = SparseMatrixCSC{Float64,Int}(I, 0, Nw)
D₂₂ = SparseMatrixCSC{Float64,Int}(I, 0, Nu)

P_SF = @inferred GeneralizedPlant{Float64,StateFeedback}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)

@test typeof(P_SF) <: GeneralizedPlant{Float64,StateFeedback}

# Test user-friendly constructor and promotion/conversion rules
@test P_SF == Plant(A, B₁, B₂, C₁, D₁₁, D₁₂)
@test P_SF == Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, I, 0, 0)
@test P_SF == Plant([A B₁ B₂; C₁ D₁₁ D₁₂], [Nx, Nx+Nu, Nw, Nu])

@test P_SF != Plant(A, B₁, B₂, C₁, D₁₁, D₁₂, I(Nx), spzeros(Nx,Nw), spzeros(Nx,Nu))

# Test if output matrices are correct 
@test P_SF.C₂ == I
@test isempty(P_SF.D₂₁)
@test isempty(P_SF.D₂₂)

## Special constructors _________________________________________________
# State-feedback plant with LQR-style unitary weights
Q = I(Nx); R = I(Nu);
C₁ = [Q; spzeros(Nu,Nx)];
D₁₁ = spzeros(Nx+Nu,Nw);
D₁₂ = [spzeros(Nx,Nu); R];

P_SF_LQR = @inferred GeneralizedPlant{Float64,StateFeedback}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)

@test P_SF_LQR == Plant(A, B₁, B₂);
@test P_SF_LQR == Plant(A, B₁, B₂, C₁, 0, D₁₂);

## Error handling and warnings __________________________________________
@test_throws ErrorException Plant([1 2], 3, 4)  # A is not square
@test_throws ErrorException Plant([1 2; 3 4], 5, [6; 7])  # Dimension mismatch w/ B₁ 
@test_throws ErrorException Plant(A, B₁, B₂, Q, D₁₁, R)  # Dimension mismatch w/ C₁ and D₁₂ 

