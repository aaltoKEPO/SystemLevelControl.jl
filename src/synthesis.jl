# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines the functions and auxiliary data-types for the 
##  synthesis of SLS controllers 

# FUNCTIONS _____________________________________________________________

function SLS_𝓗₂(P::AbstractGeneralizedPlant, 𝓢::AbstractVector; 𝓘=nothing)
    # -- --
    if P isa GeneralizedPlant{T,StateFeedback} where {T}
        # Auxiliary variables
        𝓘 = (𝓘 === nothing) ? [[i] for i in 1:P.Nx] : 𝓘;
        𝓒 = Iterators.partition( 𝓘, ceil(Int, length(𝓘)/nworkers()) );
        
        # Unpack the internal function arguments
        𝓢ₓ,𝓢ᵤ = 𝓢;
        T = length(𝓢ₓ)
        L⁺(Φ,j) = 0;

        let P=P, T=T, 𝓢ₓ=𝓢ₓ, 𝓢ᵤ=𝓢ᵤ, L⁺=L⁺
            Φ = @distributed (+) for Cⱼ in collect(𝓒)
                _SLS_𝓗₂(Cⱼ, P, T, 𝓢ₓ, 𝓢ᵤ, L⁺)
            end
            return eachcol(Φ)
        end

    end
# --
end 

function _SLS_𝓗₂(Cⱼ, P::AbstractGeneralizedPlant, T::Integer, 𝓢ₓ::AbstractVector, 𝓢ᵤ::AbstractVector, L⁺::Function)
    # Optimization loop _________________________________________________
    Φ̃ = [[spzeros(P.Nx,P.Nx) for _ in 1:T] [spzeros(P.Nu,P.Nx) for _ in 1:T]];      # SLS Mappings
    for cⱼ in Cⱼ
        (P̃,Ĩ,iiₓ,sₓ,sᵤ) = sparsity_dim_reduction(P, cⱼ, [𝓢ₓ,𝓢ᵤ]);
        Ã,B̃₁,B̃₂, C̃₁,D̃₁₁,D̃₁₂, C̃₂,D̃₂₁,D̃₂₂  = P̃;
        B̃₁ = isempty(B̃₁) ? B̃₁ : B̃₁[iiₓ,:];
        D̃₂₁ = isempty(D̃₂₁) ? D̃₂₁ : D̃₂₁[iiₓ,:];

        problem = Model(Ipopt.Optimizer); set_silent(problem)
        Φ̃ₓ = [@variable(problem, [1:P̃.Nx,1:P̃.Nw]) for _ in 1:T];
        Φ̃ᵤ = [@variable(problem, [1:P̃.Nu,1:P̃.Nw]) for _ in 1:T];
        
        @objective( problem,    Min,    norm([C̃₁,D̃₁₂]*[Φ̃ₓ,Φ̃ᵤ]*[B̃₁,D̃₂₁] + D̃₁₁, :𝓗₂) + L⁺([Φ̃ₓ,Φ̃ᵤ],cⱼ) ); # <~ L^+ is not parallelized
        @constraint(problem,                Φ̃ₓ[1]   .== Ĩ);
        @constraint(problem, [t = 1:(T-1)], Φ̃ₓ[t+1] .== Ã*Φ̃ₓ[t] + B̃₂*Φ̃ᵤ[t]);
        @constraint(problem,                   0    .== Ã*Φ̃ₓ[T] + B̃₂*Φ̃ᵤ[T]);

        for t in 1:T
            fix.(Φ̃ₓ[t][𝓢ₓ[t][sₓ,cⱼ] .≠ 1], 0.0, force=true);
            fix.(Φ̃ᵤ[t][𝓢ᵤ[t][sᵤ,cⱼ] .≠ 1], 0.0, force=true);
        end
        
        optimize!(problem)

        # TODO: Verify dimensions
        Φₓ = [sparse(vec(repeat(sₓ,P̃.Nw,1)'), vec(repeat(cⱼ,1,P̃.Nx)'), vec(value.(Φ̃ₓ[t]).*𝓢ₓ[t][sₓ,cⱼ]), P.Nx, P.Nx) for t in 1:T];
        Φᵤ = [sparse(vec(repeat(sᵤ,P̃.Nw,1)'), vec(repeat(cⱼ,1,P̃.Nu)'), vec(value.(Φ̃ᵤ[t]).*𝓢ᵤ[t][sᵤ,cⱼ]), P.Nu, P.Nx) for t in 1:T];
        Φ̃ += [Φₓ Φᵤ];
    end
    # ___________________________________________________________________
    return Φ̃
# --
end 


# OPERATOR OVERLOADS (AUXILIARY) ________________________________________

# This is working but SHOULD be done better. Specifically, it's causing some memory leak
LinearAlgebra.:*(A::Vector{SparseMatrixCSC{T,Int64}}, B::Vector{Vector{Matrix{VariableRef}}}) where {T<:Real} = hcat(A...) * vcat.(B...)
LinearAlgebra.:*(A::Vector{Matrix{AffExpr}}, B::Vector{SparseMatrixCSC{T,Int64}}) where {T<:Real} = A * vcat(B...)
LinearAlgebra.:*(A::SparseMatrixCSC{T,Int64}, B::Vector{Matrix{VariableRef}}) where {T<:Real} = [A*b for b in B]
LinearAlgebra.:*(B::Vector{Matrix{AffExpr}}, A::SparseMatrixCSC{T,Int64}) where {T<:Real} = [b*A for b in B]
LinearAlgebra.:+(B::Vector{Matrix{AffExpr}}, A::SparseMatrixCSC{T,Int64}) where {T<:Real} = [b+A for b in B]

LinearAlgebra.:norm(A::AbstractVector{T}, t::Symbol) where T = begin
    if t === :𝓗₂
        return sum([tr(Aₜ'Aₜ) for Aₜ in A]);
    else        
        throw(ArgumentError("The argument '$(t)' is not a valid norm type."));
    end
end