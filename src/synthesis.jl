# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines the functions and auxiliary data-types for the 
##  synthesis of SLS controllers 

# FUNCTIONS _____________________________________________________________

"""
    SLS(P::AbstractGeneralizedPlant, S::AbstractVector)

Synthetise a controller for the plant decribed by `P` according to the spatiotemporal
localization constraints encoded by `S`. Currently (v0), only supports the design of 
Localized \$\\mathcal{H}_2\$ state/output-feedback controllers.

# Examples 
```julia-repl
julia> 
```
"""
function SLS(P::AbstractGeneralizedPlant, S::AbstractVector; J=[nothing], norm=:H2)
    if norm == :H2
        return SLS_H2(P, S, J);
    else        
        throw(ArgumentError("Incorrent norm argument. Currently, the synthesis is implemented for norms: {:H2}"));
    end
end # -- End of SLS

# _______________________________________________________________________
# INTERNAL FUNCTIONS ____________________________________________________

function SLS_H2(P::GeneralizedPlant{<:Number,Ts}, S::AbstractVector, J::AbstractVector) where {Ts<:StateFeedback}
    # Auxiliary variables
    J = (J[1] === nothing) ? [[i] for i in 1:P.Nx] : J;
    𝓒 = Iterators.partition(J, ceil(Int, length(J)/nworkers()));
    
    # Unpack the internal function arguments
    Sₓ,Sᵤ = S;
    T = length(Sₓ)
    L⁺(Φ,j) = 0;

    let P=P, T=T, Sₓ=Sₓ, Sᵤ=Sᵤ, L⁺=L⁺
    return @distributed (+) for Cⱼ in collect(𝓒)
        _SLS_H2(Cⱼ, P, T, Sₓ, Sᵤ, L⁺)
    end
    end
end # -- End of SLS_H2 / StateFeedback

function SLS_H2(P::GeneralizedPlant{<:Number,Ts}, S::AbstractVector, J::AbstractVector) where {Ts<:OutputFeedback}
    # Auxiliary variables
    J = (J === nothing) ? [[i] for i in 1:(P.Nx+P.Ny)] : J;
    𝓒 = Iterators.partition(J, ceil(Int, length(J)/nworkers()));
    
    # Unpack the internal function arguments
    Sₓₓ,Sᵤₓ,Sₓᵧ,Sᵤᵧ = S;
    T = length(Sₓₓ)

    let P=P, T=T, Sₓₓ=Sₓₓ, Sᵤₓ=Sᵤₓ
        # ADMM
    end
end # -- End of SLS_H2 / OutputFeedback


function _SLS_H2(Cⱼ, P::AbstractGeneralizedPlant, T::Integer, 𝓢ₓ::AbstractVector, 𝓢ᵤ::AbstractVector, L⁺::Function)
    # Allocates the SLS mappings
    Φ̃ = [[spzeros(P.Nx,P.Nx) for _ in 1:T], [spzeros(P.Nu,P.Nx) for _ in 1:T]];      
    
    # Optimization loop _________________________________________________
    for cⱼ in Cⱼ
        #  Obtains a reduced-order system based on the sparsity in 𝓢
        (P̃,Ĩ,iiₓ,sₓ,sᵤ) = sparsity_dim_reduction(P, cⱼ, [𝓢ₓ,𝓢ᵤ]);   

        # Solves the reduced-order SLS problem either by solving the KKT system (ECQP)
        #  or by shipping the optimization directly to the general solver (JuMP-based)
        if length(cⱼ) == 1
            _SLS_H2_ECQP!(Φ̃, cⱼ, P̃, Ĩ, iiₓ, T, 𝓢ₓ, 𝓢ᵤ, sₓ, sᵤ)
        else
            _SLS_H2_General!(Φ̃, cⱼ, P̃, Ĩ, iiₓ, T, 𝓢ₓ, 𝓢ᵤ, sₓ, sᵤ)
        end
    end
    # ___________________________________________________________________
    return Φ̃
# --
end 

function _SLS_H2_ECQP!(Φ̃::AbstractVector, cⱼ::AbstractVector, P::AbstractGeneralizedPlant, Ĩ::AbstractMatrix, iiₓ::BitArray,  T::Integer, 𝓢ₓ::AbstractVector, 𝓢ᵤ::AbstractVector, sₓ::AbstractVector, sᵤ::AbstractVector)
    ## Creates the Hessian matrix 
    H = P.B₁[iiₓ][1]^2 * blockdiag(kron(I(T), P.C₁'P.C₁), kron(I(T), P.D₁₂'P.D₁₂));

    ## Creates the constraint matrix 
    # Dynamical constraints
    G_dyn_A =  I - kron(spdiagm(-1 => ones(T)), sparse(P.A));
    G_dyn_B = 0I - kron(spdiagm(-1 => ones(T)), sparse(P.B₂));
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
    
    for t in 1:T 
        Φ̃[1][t][sₓ,cⱼ] .= Φ[(1:P.Nx).+(t-1)*P.Nx];
        Φ̃[2][t][sᵤ,cⱼ] .= Φ[(1:P.Nu).+(t-1)*P.Nu.+T*P.Nx];
    end
# --
end 

function _SLS_H2_General!(Φ̃::AbstractVector, cⱼ::AbstractVector, P::AbstractGeneralizedPlant, Ĩ::AbstractMatrix, iiₓ::BitArray, T::Integer, 𝓢ₓ::AbstractVector, 𝓢ᵤ::AbstractVector, sₓ::AbstractVector, sᵤ::AbstractVector)
    # Retrieves the reduced-order state-space matrices
    A,B₁,B₂, C₁,D₁₁,D₁₂, C₂,D₂₁,D₂₂  = P;
    B₁ = isempty(B₁) ? B₁ : B₁[iiₓ,:];
    D₂₁ = isempty(D₂₁) ? D₂₁ : D₂₁[iiₓ,:];

    # Designs and solves the OCP associated with subsystem P̃
    problem = Model(Ipopt.Optimizer); set_silent(problem)
    Φₓ = [@variable(problem, [1:P.Nx,1:P.Nw]) for _ in 1:T];
    Φᵤ = [@variable(problem, [1:P.Nu,1:P.Nw]) for _ in 1:T];
    
    T_zw = _create_SLS_ref_operator(problem, [C₁ D₁₂], Φₓ, Φᵤ, [B₁; D₂₁], D₁₁);

    @objective(problem,      Min,      norm(T_zw, :𝓗₂));
    @constraint(problem,                Φₓ[1]   .== Ĩ);
    @constraint(problem, [t = 1:(T-1)], Φₓ[t+1] .== A*Φₓ[t] + B₂*Φᵤ[t]);
    @constraint(problem,                   0    .== A*Φₓ[T] + B₂*Φᵤ[T]);

    for t in 1:T
        fix.(Φₓ[t][𝓢ₓ[t][sₓ,cⱼ] .≠ 1], 0.0, force=true);
        fix.(Φᵤ[t][𝓢ᵤ[t][sᵤ,cⱼ] .≠ 1], 0.0, force=true);
    end
    
    optimize!(problem)

    for t in 1:T 
        Φ̃[1][t][sₓ,cⱼ] = value.(Φₓ[t]) .* 𝓢ₓ[t][sₓ,cⱼ];
        Φ̃[2][t][sᵤ,cⱼ] = value.(Φᵤ[t]) .* 𝓢ᵤ[t][sᵤ,cⱼ];
    end
end

# _______________________________________________________________________
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