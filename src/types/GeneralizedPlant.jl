# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
abstract type AbstractGeneralizedPlant{T<:Number} end
# --

# PLANT DATA TYPE DECLARATIONS _______________________________________________
struct GeneralizedPlant{T} <: AbstractGeneralizedPlant{T}
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

    function GeneralizedPlant{T}(A::SparseMatrixCSC{T,Int64},B₁::SparseMatrixCSC{T,Int64},B₂::SparseMatrixCSC{T,Int64},
                                 C₁::SparseMatrixCSC{T,Int64},D₁₁::SparseMatrixCSC{T,Int64},D₁₂::SparseMatrixCSC{T,Int64},
                                 C₂::SparseMatrixCSC{T,Int64},D₂₁::SparseMatrixCSC{T,Int64},D₂₂::SparseMatrixCSC{T,Int64}) where {T}
        
        validate_GeneralizedPlant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
        new{T}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂, size(A,1), size(C₁,1), size(C₂,1), size(B₁,2), size(B₂,2))
    end

    function GeneralizedPlant(A::AbstractArray,B₁::AbstractArray,B₂::AbstractArray,
                              C₁::AbstractArray,D₁₁::AbstractArray,D₁₂::AbstractArray,
                              C₂::AbstractArray,D₂₁::AbstractArray,D₂₂::AbstractArray)
        
        
        validate_GeneralizedPlant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
        new{Float64}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂, size(A,1), size(C₁,1), size(C₂,1), size(B₁,2), size(B₂,2))
    end
end

# Custom constructors
function GeneralizedPlant(Σ::AbstractArray{T}, DIMS::AbstractArray{Int64}) where {T<:Number}
    if length(DIMS) == 4
        StateFeedbackPlant(Σ, DIMS)
    elseif length(DIMS) == 5
        nx, nz, ny, nu, nw = DIMS;
        if (nx+nz+ny) ≠ size(Σ,1) || (nx+nw+nu) ≠ size(Σ,2)
            error("Dimensions do not match! Expected: ($(nx+nz+ny),$(nx+nu+nw)), got $(size(Σ))")
        end

        Σ = sparse(Σ);
        A = Σ[1:nx,           1:nx];   B₁ = Σ[1:nx,           nx.+(1:nw)];   B₂ = Σ[1:nx,            (nx+nw).+(1:nu)];
        C₁ = Σ[nx.+(1:nz),     1:nx];  D₁₁ = Σ[nx.+(1:nz),     nx.+(1:nw)];  D₁₂ = Σ[nx.+(1:nz),      (nx+nw).+(1:nu)];
        C₂ = Σ[(nx+nz).+(1:ny),1:nx];  D₂₁ = Σ[(nx+nz).+(1:ny),nx.+(1:nw)];  D₂₂ = Σ[(nx+nz).+(1:ny), (nx+nw).+(1:nu)];
        
        GeneralizedPlant{T}(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
    end
end

GeneralizedPlant(A::AbstractArray,B₁::AbstractArray,B₂::AbstractArray,C₁::AbstractArray,D₁₁::AbstractArray,D₁₂::AbstractArray) = StateFeedbackPlant(A, B₁, B₂, C₁, D₁₁, D₁₂)
GeneralizedPlant(A::AbstractArray,B₁::AbstractArray,B₂::AbstractArray) = StateFeedbackPlant(A, B₁, B₂)

# User-friendly constructor 
Plant(args...; kwargs...) = GeneralizedPlant(args...; kwargs...)


# VALIDATIONS AND AUXILIARY FUNCTIONS __________________________________________
function validate_GeneralizedPlant(A, B₁, B₂, C₁, D₁₁, D₁₂, C₂, D₂₁, D₂₂)
    nx, nw, nu, nz, ny = (size(A,1), size(B₁,2), size(B₂,2), size(C₁,1), size(C₂,1))

    if size(A, 2) != nx && nx != 0
        error("A has dimensions $(size(A)), but must be square.")
    elseif size(B₁, 1) != nx || size(B₂, 1) != nx
        error("The number of rows of A ($(size(A,1))) does not match either B₁ ($(size(B₁,1))) or B₂ ($(size(B₂,1))).")
    elseif size(C₁, 2) != nx || size(C₂, 2) != nx
        error("The number of columns of A ($(size(A,2))) does not match either C₁ ($(size(C₁,2))) or C₂ ($(size(C₂,2))).")
    elseif size(D₁₁, 1) != nz || size(D₁₂, 1) != nz
        error("The number of rows of C₁ ($(size(C₁,1))) does not match either D₁₁ ($(size(D₁₁,1))) or D₁₂ ($(size(D₁₂,1))).")
    elseif size(D₁₁, 2) != nw || size(D₂₁, 2) != nw
        error("The number of columns of B₁ ($(size(B₁,2))) does not match either D₁₁ ($(size(D₁₁,2))) or D₂₁ ($(size(D₂₁,2))).")
    elseif size(D₂₁, 1) != ny || size(D₂₂, 1) != ny
        error("The number of rows of C₂ ($(size(C₂,1))) does not match either D₂₁ ($(size(D₂₁,1))) or D₂₂ ($(size(D₂₂,1))).")
    elseif size(D₁₂, 2) != nu || size(D₂₂, 2) != nu
        error("The number of columns of B₂ ($(size(B₂,2))) does not match either D₁₂ ($(size(D₁₂,2))) or D₂₂ ($(size(D₂₂,2))).")
    end
end

Base.show(io::IO, P::AbstractGeneralizedPlant{T}) where {T} = begin 
    # io_Buffer = IOBuffer();
    # fNames = [:A :B₁ :B₂; :C₁ :D₁₁ :D₁₂; :C₂ :D₂₁ :D₂₂];
    # for i = 1:3, j = 1:3
    # end
    if isa(P, AbstractStateFeedbackPlant)
        Σ = [P.A P.B₁ P.B₂; P.C₁ P.D₁₁ P.D₁₂];
        print(io, "$(size(Σ,1))×$(size(Σ,2)) $(typeof(P)) w/ $(size(P.A,1)) states and $(size(P.B₂,2)) controls.")
    else
        Σ = [P.A P.B₁ P.B₂; P.C₁ P.D₁₁ P.D₁₂; P.C₂ P.D₂₁ P.D₂₂];
        print(io, "$(size(Σ,1))×$(size(Σ,2)) $(typeof(P)) w/ $(size(P.A,1)) states, $(size(P.B₂,2)) controls, $(size(P.C₂,1)) outputs.")
    end
end

# USEFUL ALIASES _______________________________________________________________


# ______________________________________________________________________________
