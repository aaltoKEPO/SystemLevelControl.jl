# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines the functions and auxiliary data-types for the 
##  synthesis of SLS controllers 

# FUNCTIONS _____________________________________________________________

function SLS(P::AbstractGeneralizedPlant, S::AbstractVector; 𝓘=nothing, norm=:H2)
    # -- --
    if P isa GeneralizedPlant{T,StateFeedback} where {T}
        # Auxiliary variables
        𝓘 = (𝓘 === nothing) ? [[i] for i in 1:P.Nx] : 𝓘;
        𝓒 = Iterators.partition(𝓘, ceil(Int, length(𝓘)/nworkers()));
        
        # Unpack the internal function arguments
        Sₓ,Sᵤ = S;
        T = length(Sₓ)
        L⁺(Φ,j) = 0;

        let P=P, T=T, Sₓ=Sₓ, Sᵤ=Sᵤ, L⁺=L⁺
            Φ = @distributed (+) for Cⱼ in collect(𝓒)
                _SLS_H₂(Cⱼ, P, T, Sₓ, Sᵤ, L⁺)
            end
            return eachcol(Φ)
        end

    end
# --
end 

function _SLS_H₂(Cⱼ, P::AbstractGeneralizedPlant, T::Integer, 𝓢ₓ::AbstractVector, 𝓢ᵤ::AbstractVector, L⁺::Function)
    # Optimization loop _________________________________________________
    Φ̃ = [[spzeros(P.Nx,P.Nx) for _ in 1:T] [spzeros(P.Nu,P.Nx) for _ in 1:T]];      # SLS Mappings
    for cⱼ in Cⱼ
        (P̃,Ĩ,iiₓ,sₓ,sᵤ) = sparsity_dim_reduction(P, cⱼ, [𝓢ₓ,𝓢ᵤ]);   #  Obtains a reduced-order system based on the sparsity in 𝓢

        if length(cⱼ) == 1
            Φ̃ += _SLS_H₂_ECQP(cⱼ, Ĩ, P̃, T, 𝓢ₓ, 𝓢ᵤ, sₓ, sᵤ)
        else
            # Retrieves the reduced-order state-space matrices
            Ã,B̃₁,B̃₂, C̃₁,D̃₁₁,D̃₁₂, C̃₂,D̃₂₁,D̃₂₂  = P̃;
            B̃₁ = isempty(B̃₁) ? B̃₁ : B̃₁[iiₓ,:];
            D̃₂₁ = isempty(D̃₂₁) ? D̃₂₁ : D̃₂₁[iiₓ,:];

            # Designs and solves the OCP associated with subsystem P̃
            problem = Model(Ipopt.Optimizer); set_silent(problem)
            Φ̃ₓ = [@variable(problem, [1:P̃.Nx,1:P̃.Nw]) for _ in 1:T];
            Φ̃ᵤ = [@variable(problem, [1:P̃.Nu,1:P̃.Nw]) for _ in 1:T];
            
            H_w2z = _create_SLS_ref_operator(problem, [C̃₁ D̃₁₂], Φ̃ₓ, Φ̃ᵤ, [B̃₁; D̃₂₁], D̃₁₁);
    
            @objective(problem,      Min,      norm(H_w2z, :𝓗₂) + L⁺([Φ̃ₓ,Φ̃ᵤ],cⱼ)); # <~ L^+ is not parallelized
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
    end
    # ___________________________________________________________________
    return Φ̃
# --
end 

function _SLS_H₂_ECQP(cⱼ, Ĩ::AbstractMatrix, P::AbstractGeneralizedPlant, T::Integer, 𝓢ₓ::AbstractVector, 𝓢ᵤ::AbstractVector, sₓ::AbstractVector, sᵤ::AbstractVector)
    ## Creates the Hessian matrix 
    H = P.B₁[Bool.(Ĩ)][1]^2 * blockdiag(kron(I(T), P.C₁'P.C₁), kron(I(T), P.D₁₂'P.D₁₂));

    ## Creates the constraint matrix 
    # Dynamical constraints
    G_dyn_A =  I-kron(spdiagm(-1 => ones(T)), sparse(P.A));
    G_dyn_B = 0I-kron(spdiagm(-1 => ones(T)), sparse(P.B₂));
    G_dyn = [G_dyn_A[:,1:(P.Nx*T)]  G_dyn_B[:,1:(P.Nu*T)]];
    
    # Sparsity constraints 
    S_idx = [vcat([         (t-1)*P.Nx .+ findall(iszero, 𝓢ₓ[t][sₓ,cⱼ[1]]) for t in 2:T]...);
             vcat([T*P.Nx + (t-1)*P.Nu .+ findall(iszero, 𝓢ᵤ[t][sᵤ,cⱼ[1]]) for t in 1:T]...)]

    G_sp = sparse(1:length(S_idx), S_idx, ones(length(S_idx)), length(S_idx), (P.Nx+P.Nu)*T); 

    # -
    G = [G_dyn; G_sp];
    g = [Ĩ; zeros(size(G,1)-P.Nx)];

    # Solves system of equations 
    Φ = qr([H G'; G 0I]) \ Array([zeros(size(H,1)); g]);
    
    Φₓ = [sparse(sₓ, repeat(cⱼ,P.Nx), Φ[(1:P.Nx).+(t-1)*P.Nx], size(𝓢ₓ,1), size(𝓢ₓ,2)) for t in 1:T];
    Φᵤ = [sparse(sᵤ, repeat(cⱼ,P.Nu), Φ[(1:P.Nu).+(t-1)*P.Nu.+T*P.Nx], size(𝓢ᵤ,1), size(𝓢ᵤ,2)) for t in 1:T];

    # ___________________________________________________________________
    return [Φₓ Φᵤ]
# --
end 


# OPERATOR OVERLOADS / AUXILIARY FUNCTIONS ______________________________
function _create_SLS_ref_operator(problem::Model, L::AbstractMatrix, Φ̃ₓ::Vector{Matrix{VariableRef}}, Φ̃ᵤ::Vector{Matrix{VariableRef}}, R::AbstractMatrix, D::AbstractMatrix)
    return [@expression(problem, L*[Φ[1];Φ[2]]*R + D) for Φ in zip(Φ̃ₓ,Φ̃ᵤ)]
end

LinearAlgebra.:norm(A::AbstractVector{T}, t::Symbol) where T = begin
    if t === :𝓗₂
        return sum([tr(Aₜ'Aₜ) for Aₜ in A]) / (2π);
    else        
        throw(ArgumentError("The argument '$(t)' is not a valid norm type."));
    end
end