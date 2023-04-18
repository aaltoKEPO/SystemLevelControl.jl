# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
abstract type AbstractGeneralizedPlant{T<:Number,Ts<:AbstractFeedbackStructure} end
# --

# AUXILIARY TYPE ALIASES
const NumberOrAbstractArray = Union{Number, AbstractArray}
const EyeOrNumberOrAbstractArray = Union{UniformScaling, NumberOrAbstractArray}

# PLANT DATA TYPE DECLARATIONS _______________________________________________
struct GeneralizedPlant{T,Ts} <: AbstractGeneralizedPlant{T,Ts}
    # Fields
    A::SparseMatrixCSC{T,Int64}
    B₁::SparseMatrixCSC{T,Int64}
    B₂::SparseMatrixCSC{T,Int64}
    C₁::SparseMatrixCSC{T,Int64}
    D₁₁::SparseMatrixCSC{T,Int64}
    D₁₂::SparseMatrixCSC{T,Int64}
    C₂::SparseMatrixCSC{T,Int64}
    D₂₁::SparseMatrixCSC{T,Int64}
    D₂₂::SparseMatrixCSC{T,Int64}
    Nx::Integer; Nz::Integer; Ny::Integer; Nw::Integer; Nu::Integer;

    # Explicit constructor
    function GeneralizedPlant{T,Ts}(A::SparseMatrixCSC{T,Int64},B₁::SparseMatrixCSC{T,Int64},B₂::SparseMatrixCSC{T,Int64},
                                 C₁::SparseMatrixCSC{T,Int64},D₁₁::SparseMatrixCSC{T,Int64},D₁₂::SparseMatrixCSC{T,Int64},
                                 C₂::SparseMatrixCSC{T,Int64},D₂₁::SparseMatrixCSC{T,Int64},D₂₂::SparseMatrixCSC{T,Int64}
    ) where {T<:Number,Ts<:AbstractFeedbackStructure}
        
        validate_GeneralizedPlant(Ts, A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
        new{T,Ts}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂, size(A,1), size(C₁,1), size(C₂,1), size(B₁,2), size(B₂,2))
    end
end

# Custom constructors
function GeneralizedPlant(A::NumberOrAbstractArray,B₁::NumberOrAbstractArray,B₂::NumberOrAbstractArray,
                          C₁::NumberOrAbstractArray,D₁₁::NumberOrAbstractArray,D₁₂::NumberOrAbstractArray,
                          C₂::EyeOrNumberOrAbstractArray,D₂₁::NumberOrAbstractArray,D₂₂::NumberOrAbstractArray
)
    # Firstly, obtains the most generic eltype from all the matrices and verify feedback structure
    T = promote_type(eltype(A), eltype(B₁), eltype(B₂), eltype(C₁), eltype(D₁₁), eltype(D₁₂), eltype(C₂), eltype(D₂₁), eltype(D₂₂))
    Ts = (C₂ == I && (isempty(D₂₁) || D₂₁==0)) ? StateFeedback : OutputFeedback;

    # Converts all numbers and vectors to sparse matrices 
    A = to_sparse_matrix(T,A)
    B₁ = to_sparse_matrix(T, B₁)
    B₂ = to_sparse_matrix(T, B₂)

    C₁ = to_sparse_matrix(T, C₁)
    D₁₁ = fix_feedthrough(to_sparse_matrix(T, D₁₁), size(C₁,1), size(B₁,2))
    D₁₂ = to_sparse_matrix(T, D₁₂)

    if Ts <: OutputFeedback 
        C₂ = (C₂ isa UniformScaling) ? SparseMatrixCSC{T,Int64}(C₂, size(A)) : to_sparse_matrix(T, C₂)
        D₂₁ = to_sparse_matrix(T, D₂₁)
        D₂₂ = fix_feedthrough(to_sparse_matrix(T, D₂₂), size(C₂,1), size(B₂,2))
    else 
        C₂ = SparseMatrixCSC{T,Int64}(I, size(A))
        D₂₁ = SparseMatrixCSC{T,Int64}(I, 0, size(B₁,2))
        D₂₂ = SparseMatrixCSC{T,Int64}(I, 0, size(B₂,2))
    end

    # Returns the GeneralizedPlant object
    return GeneralizedPlant{T,Ts}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
end

function GeneralizedPlant(A::NumberOrAbstractArray,B₁::NumberOrAbstractArray,B₂::NumberOrAbstractArray,C₁::NumberOrAbstractArray,D₁₁::NumberOrAbstractArray,D₁₂::NumberOrAbstractArray)
    return GeneralizedPlant(A,B₁,B₂,C₁,D₁₁,D₁₂,I,Bool[],Bool[])
end

function GeneralizedPlant(A::NumberOrAbstractArray,B₁::NumberOrAbstractArray,B₂::NumberOrAbstractArray)
    CD₁ = SparseMatrixCSC{Bool,Int64}(I, size(A,1)+size(B₂,2), size(A,1)+size(B₂,2))
    C₁  = CD₁[:,1:size(A,1)]
    D₁₂ = CD₁[:,1+size(A,1):end]
    return GeneralizedPlant(A,B₁,B₂,C₁,0,D₁₂,I,Bool[],Bool[])
end

function GeneralizedPlant(Σ::AbstractArray{T}, DIMS::AbstractArray{Int64}) where {T<:Number}
    if length(DIMS) == 5
        Ts = OutputFeedback 
        nx, nz, ny, nw, nu = DIMS;
    else
        Ts = StateFeedback 
        nx, nz, nw, nu = DIMS;
        ny = 0
    end
 
    if (nx+nz+ny) ≠ size(Σ,1) || (nx+nw+nu) ≠ size(Σ,2)
        error("Dimensions mismatch! Expected: ($(nx+nz+ny),$(nx+nu+nw)), got $(size(Σ))")
    end

    Σ = sparse(Σ);
    A  = Σ[1:nx,                1:nx];   B₁ = Σ[1:nx,                (nx+1):(nx+nw)];   B₂ = Σ[1:nx,                 (nx+nw+1):(nx+nw+nu)];
    C₁ = Σ[(nx+1):(nx+nz),      1:nx];  D₁₁ = Σ[(nx+1):(nx+nz),      (nx+1):(nx+nw)];  D₁₂ = Σ[(nx+1):(nx+nz),       (nx+nw+1):(nx+nw+nu)];
    C₂ = Σ[(nx+nz+1):(nx+nz+ny),1:nx];  D₂₁ = Σ[(nx+nz+1):(nx+nz+ny),(nx+1):(nx+nw)];  D₂₂ = Σ[(nx+nz+1):(nx+nz+ny), (nx+nw+1):(nx+nw+nu)];
    
    C₂ = isempty(C₂) ? SparseMatrixCSC{T,Int64}(I, size(A)) : C₂;

    GeneralizedPlant{T,Ts}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
end

# User-friendly constructor interface
Plant(args...; kwargs...) = GeneralizedPlant(args...; kwargs...)


# VALIDATIONS AND AUXILIARY FUNCTIONS __________________________________________
function validate_GeneralizedPlant(Ts, A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
    C₂,D₂₁,D₂₂ = (Ts <: StateFeedback) ? (A,B₁,B₂) : (C₂,D₂₁,D₂₂);

    nx, nw, nu, nz, ny = (size(A,1), size(B₁,2), size(B₂,2), size(C₁,1), size(C₂,1))

    if size(A, 2) != nx && 
        error("A must be nonempty and square, but has dimensions ($(size(A,1))×$(size(A,2))).")
    elseif size(B₁, 1) != nx || size(B₂, 1) != nx
        error("The number of rows of A (=$(size(A,1))) does not match either B₁ (=$(size(B₁,1))) or B₂ (=$(size(B₂,1))).")
    elseif size(C₁, 2) != nx || size(C₂, 2) != nx
        error("The number of columns of A (=$(size(A,2))) does not match either C₁ (=$(size(C₁,2))) or C₂ (=$(size(C₂,2))).")
    elseif size(D₁₁, 1) != nz || size(D₁₂, 1) != nz
        error("The number of rows of C₁ (=$(size(C₁,1))) does not match either D₁₁ (=$(size(D₁₁,1))) or D₁₂ (=$(size(D₁₂,1))).")
    elseif size(D₁₁, 2) != nw || size(D₂₁, 2) != nw
        error("The number of columns of B₁ (=$(size(B₁,2))) does not match either D₁₁ (=$(size(D₁₁,2))) or D₂₁ (=$(size(D₂₁,2))).")
    elseif size(D₂₁, 1) != ny || size(D₂₂, 1) != ny
        error("The number of rows of C₂ (=$(size(C₂,1))) does not match either D₂₁ (=$(size(D₂₁,1))) or D₂₂ (=$(size(D₂₂,1))).")
    elseif size(D₁₂, 2) != nu || size(D₂₂, 2) != nu
        error("The number of columns of B₂ (=$(size(B₂,2))) does not match either D₁₂ (=$(size(D₁₂,2))) or D₂₂ (=$(size(D₂₂,2))).")
    end
end

function Base.show(io::IO, P::AbstractGeneralizedPlant)
    # io_Buffer = IOBuffer();
    # fNames = [:A :B₁ :B₂; :C₁ :D₁₁ :D₁₂; :C₂ :D₂₁ :D₂₂];
    # for i = 1:3, j = 1:3
    # end
    if isa(P, AbstractGeneralizedPlant{T,StateFeedback} where T)
        Σ = [P.A P.B₁ P.B₂; P.C₁ P.D₁₁ P.D₁₂];
        print(io, "$(size(Σ,1))×$(size(Σ,2)) $(typeof(P)) w/ $(size(P.A,1)) states and $(size(P.B₂,2)) controls.")
    else
        Σ = [P.A P.B₁ P.B₂; P.C₁ P.D₁₁ P.D₁₂; P.C₂ P.D₂₁ P.D₂₂];
        print(io, "$(size(Σ,1))×$(size(Σ,2)) $(typeof(P)) w/ $(size(P.A,1)) states, $(size(P.B₂,2)) controls, $(size(P.C₂,1)) outputs.")
    end
end


# ______________________________________________________________________________
