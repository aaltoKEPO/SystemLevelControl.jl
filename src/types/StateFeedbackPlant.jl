# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
abstract type AbstractStateFeedbackPlant{T} <: AbstractGeneralizedPlant{T} end
# --

# PLANT DATA TYPE DECLARATIONS _______________________________________________
struct StateFeedbackPlant{T} <: AbstractStateFeedbackPlant{T}
    A::SparseMatrixCSC{T,Int64}
    B₁::SparseMatrixCSC{T,Int64}
    B₂::SparseMatrixCSC{T,Int64}
    C₁::SparseMatrixCSC{T,Int64}
    D₁₁::SparseMatrixCSC{T,Int64}
    D₁₂::SparseMatrixCSC{T,Int64}
    C₂::UniformScaling{Bool}
    D₂₁::SparseMatrixCSC{T,Int64}
    D₂₂::SparseMatrixCSC{T,Int64}
    Nx::Integer; Nz::Integer; Ny::Integer; Nw::Integer; Nu::Integer;

    function StateFeedbackPlant{T}(A::SparseMatrixCSC{T,Int64},B₁::SparseMatrixCSC{T,Int64},B₂::SparseMatrixCSC{T,Int64},
                                   C₁::SparseMatrixCSC{T,Int64},D₁₁::SparseMatrixCSC{T,Int64},D₁₂::SparseMatrixCSC{T,Int64}) where {T}
        validate_GeneralizedPlant(A, B₁, B₂, C₁, D₁₁, D₁₂,  A, B₁, B₂)
        D₂₁ = SparseMatrixCSC{T,Int64}(I, 0, size(B₁,2))
        D₂₂ = SparseMatrixCSC{T,Int64}(I, 0, size(B₂,2))
        new{T}(A, B₁, B₂, C₁, D₁₁, D₁₂, I, D₂₁, D₂₂, size(A,1), size(C₁,1), size(A,1), size(D₂₁,2), size(D₂₂,2))
    end

    function StateFeedbackPlant(A::AbstractMatrix,B₁::AbstractMatrix,B₂::AbstractMatrix,
                                C₁::AbstractMatrix,D₁₁::AbstractMatrix,D₁₂::AbstractMatrix)
        validate_GeneralizedPlant(A, B₁, B₂, C₁, D₁₁, D₁₂,  A, B₁, B₂)
        D₂₁ = SparseMatrixCSC{Float64,Int64}(I, 0, size(B₁,2))
        D₂₂ = SparseMatrixCSC{Float64,Int64}(I, 0, size(B₂,2))
        new{Float64}(A, B₁, B₂, C₁, D₁₁, D₁₂, I, D₂₁, D₂₂, size(A,1), size(C₁,1), size(A,1), size(D₂₁,2), size(D₂₂,2))
    end
end

# Custom constructors
function StateFeedbackPlant(Σ::Union{SparseMatrixCSC{T,Int64}, AbstractMatrix{T}}, DIMS::Vector{Int64}) where {T}
    nx, nz, nu, nw = DIMS;
    if (nx+nz) ≠ size(Σ,1) || (nx+nw+nu) ≠ size(Σ,2)
        error("Dimensions do not match! Expected: ($(nx+nz),$(nx+nu+nw)), got $(size(Σ))")
    end

    A  = sparse(Σ[1:nx,       1:nx]);   B₁ = sparse(Σ[1:nx,       nx.+(1:nw)]);   B₂ = sparse(Σ[1:nx,       (nx+nw).+(1:nu)]);
    C₁ = sparse(Σ[nx.+(1:nz), 1:nx]);  D₁₁ = sparse(Σ[nx.+(1:nz), nx.+(1:nw)]);  D₁₂ = sparse(Σ[nx.+(1:nz), (nx+nw).+(1:nu)]);
    StateFeedbackPlant(A, B₁, B₂, C₁, D₁₁, D₁₂)
end
