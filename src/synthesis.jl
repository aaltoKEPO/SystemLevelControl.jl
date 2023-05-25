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
    if norm === :H2
        return SLS_H2(P, S, J);
    else        
        throw(ArgumentError("Incorrent norm argument. Currently, the synthesis is implemented for norms: {:H2}"));
    end
end # -- End of SLS

# _______________________________________________________________________
# INTERNAL FUNCTIONS ____________________________________________________

function SLS_H2(P::GeneralizedPlant{<:Number,Ts}, S::AbstractVector, J::AbstractVector) where {Ts<:StateFeedback}
    # Unpack the internal function arguments
    Sₓ,Sᵤ = S;
    T = length(Sₓ)

    # Auxiliary variables
    J = (J[1] === nothing) ? [[i] for i in 1:P.Nx] : J;
    𝓒 = Iterators.partition(J, ceil(Int, length(J)/nworkers()));

    let P=P, T=T, Sₓ=Sₓ, Sᵤ=Sᵤ
    return @distributed (+) for Cⱼ in collect(𝓒)
        _SLS_H2(Cⱼ, P, T, Sₓ, Sᵤ, [nothing], 0)
    end
    end
end # -- End of SLS_H2 / StateFeedback

function SLS_H2(P::GeneralizedPlant{<:Number,Ts}, S::AbstractVector, J::AbstractVector) where {Ts<:OutputFeedback}
    # Unpack the internal function arguments
    Sₓₓ,Sᵤₓ,Sₓᵧ,Sᵤᵧ = S;
    T = length(Sₓₓ)
    
    # Auxiliary variables
    if J[1] === nothing 
        J = [[[i] for i in 1:(P.Nx+P.Ny)],
             [[i] for i in 1:(P.Nx+P.Nu)]]
    end

    𝓒ₓ = Iterators.partition(J[1], ceil(Int, length(J[1])/nworkers()));
    𝓒ᵧ = Iterators.partition(J[2], ceil(Int, length(J[2])/nworkers()));
    
    # Adjusts the sparsity constraints for the primal and dual problems
    Sₓ,Sᵤ = (hcat.(Sₓₓ,Sₓᵧ), hcat.(Sᵤₓ,Sᵤᵧ));
    Sₓ_a,Sᵤ_a = (vec(hcat.(Sₓₓ',Sᵤₓ')), vec(hcat.(Sₓᵧ',Sᵤᵧ')));
    
    ρ = 50;
    Λ = [[spzeros(size(_Sₓ)) for _Sₓ in Sₓ], [spzeros(size(_Sᵤ)) for _Sᵤ in Sᵤ]];
    Ψ = [[spzeros(size(_Sₓ)) for _Sₓ in Sₓ], [spzeros(size(_Sᵤ)) for _Sᵤ in Sᵤ]];

    # let P=P, T=T, Sₓ=Sₓ, Sᵤ=Sᵤ, Sₓ_a=Sₓ_a, Sᵤ_a=Sᵤ_a, Λ=Λ, Ψ=Ψ
    for _ in 1:100
        ν = vec(vcat.(Ψ...)) - vcat.(Λ...);
        Φ = @distributed (+) for Cⱼ in collect(𝓒ₓ)
            _SLS_H2(Cⱼ, P, T, Sₓ, Sᵤ, ν, ρ)
        end

        ν = vec(vcat.(Φ...)' - vcat.(Λ...)');
        Ψ = @distributed (+) for Cⱼ in collect(𝓒ᵧ) 
            _SLS_H2(Cⱼ, P', T, Sₓ_a, Sᵤ_a, ν, ρ)
        end

        Ψ = [[[Ψ̃[1][:,1:P.Nx]'     Ψ̃[2][:,1:P.Nx]']     for Ψ̃ in zip(Ψ[1],Ψ[2])], 
             [[Ψ̃[1][:,1+P.Nx:end]' Ψ̃[2][:,1+P.Nx:end]'] for Ψ̃ in zip(Ψ[1],Ψ[2])]]

        Λ += (Φ - Ψ)
    end
    # end # ADMM

    return Λ
end # -- End of SLS_H2 / OutputFeedback


function _SLS_H2(Cⱼ::AbstractVector, P::AbstractGeneralizedPlant, T::Integer, Sₓ::AbstractVector, Sᵤ::AbstractVector, ν::AbstractVector, ρ::Real)
    # Allocates the SLS mappings
    Φ̃ = [[spzeros(size(_Sₓ)) for _Sₓ in Sₓ], [spzeros(size(_Sᵤ)) for _Sᵤ in Sᵤ]];
    
    # Optimization loop _________________________________________________
    for cⱼ in Cⱼ
        # Obtains a reduced-order system based on the sparsity in 𝓢 = [Sₓ, Sᵤ]
        (P̃,Ĩ,iᵣ,sₓ,sᵤ) = sparsity_dim_reduction(P, cⱼ, Sₓ, Sᵤ);  

        # Slices the ADMM constant term (if needed)
        if ν[1] === nothing
            ν̃ = [[spzeros(P̃.Nx,P̃.Nw) for _ in Sₓ], 
                 [spzeros(P̃.Nu,P̃.Nw) for _ in Sᵤ]];
        else
            ν̃ = [[_ν[sₓ,cⱼ] for _ν in ν], 
                 [_ν[P.Nx.+sᵤ,cⱼ] for _ν in ν]];
        end

        # Solves the reduced-order SLS problem either by solving the KKT system (ECQP)
        #  or by shipping the optimization directly to the general solver (JuMP-based)
        if length(cⱼ) == 1
            _SLS_H2_ECQP!(Φ̃, cⱼ, P̃, Ĩ, iᵣ, T, Sₓ, Sᵤ, sₓ, sᵤ, ν̃, ρ)
        else
            _SLS_H2_General!(Φ̃, cⱼ, P̃, Ĩ, iᵣ, T, Sₓ, Sᵤ, sₓ, sᵤ, ν̃, ρ)
        end
    end
    # ___________________________________________________________________
    return Φ̃
# --
end 

function _SLS_H2_ECQP!(Φ̃::AbstractVector, cⱼ::AbstractVector, P::AbstractGeneralizedPlant, Ĩ::AbstractMatrix, iᵣ::BitArray, T::Integer, Sₓ::AbstractVector, Sᵤ::AbstractVector, sₓ::AbstractVector, sᵤ::AbstractVector, ν::AbstractVector, ρ::Real)
    ## Creates the Hessian matrix 
    σ = [P.B₁; P.D₂₁][iᵣ][1]
    H = σ^2 * blockdiag(kron(I(T), P.C₁'P.C₁), kron(I(T), P.D₁₂'P.D₁₂));

    # Unpacks the ADMM constant term (or creates a vector of zeros, if state-feedback)
    ν = vec([vcat(ν[1]...); vcat(ν[2]...)]);

    ## Creates the constraint matrix 
    # Dynamical constraints
    G_dyn_A =  I - kron(spdiagm(-1 => ones(T)), sparse(P.A));
    G_dyn_B = 0I - kron(spdiagm(-1 => ones(T)), sparse(P.B₂));
    G_dyn = [G_dyn_A[:,1:(P.Nx*T)]  G_dyn_B[:,1:(P.Nu*T)]];
    
    # Sparsity constraints 
    S_idx = [vcat([         (t-1)*P.Nx .+ findall(iszero, Sₓ[t][sₓ,cⱼ[1]]) for t in 2:T]...);
             vcat([T*P.Nx + (t-1)*P.Nu .+ findall(iszero, Sᵤ[t][sᵤ,cⱼ[1]]) for t in 1:T]...)]

    G_sp = sparse(1:length(S_idx), S_idx, ones(length(S_idx)), length(S_idx), (P.Nx+P.Nu)*T); 

    # -
    G = [G_dyn; G_sp];
    g = [Ĩ; zeros(size(G,1)-P.Nx)];

    # Solves system of equations 
    Φ = qr([H+ρ*I G'; G 0I]) \ Array([ρ*ν; g]);
    
    for t in 1:T 
        Φ̃[1][t][sₓ,cⱼ] .= Sₓ[t][sₓ,cⱼ] .* Φ[(1:P.Nx).+(t-1)*P.Nx] ;
        Φ̃[2][t][sᵤ,cⱼ] .= Sᵤ[t][sᵤ,cⱼ] .* Φ[(1:P.Nu).+(t-1)*P.Nu.+T*P.Nx] ;
    end
# --
end 

function _SLS_H2_General!(Φ̃::AbstractVector, cⱼ::AbstractVector, P::AbstractGeneralizedPlant, Ĩ::AbstractMatrix, iᵣ::BitArray, T::Integer, 𝓢ₓ::AbstractVector, 𝓢ᵤ::AbstractVector, sₓ::AbstractVector, sᵤ::AbstractVector, ν::AbstractVector, ρ::Real)
    # Retrieves the reduced-order state-space matrices
    A,B₁,B₂, C₁,D₁₁,D₁₂, C₂,D₂₁,D₂₂  = P;

    # Designs and solves the OCP associated with subsystem P̃
    problem = Model(Ipopt.Optimizer); set_silent(problem)
    Φₓ = [@variable(problem, [1:P.Nx,1:P.Nw]) for _ in 1:T];
    Φᵤ = [@variable(problem, [1:P.Nu,1:P.Nw]) for _ in 1:T];
    
    T_zw = _create_SLS_ref_operator(problem, [C₁ D₁₂], Φₓ, Φᵤ, [B₁; D₂₁][iᵣ,:], D₁₁);

    @objective(problem,      Min,       norm(T_zw, :𝓗₂) + 0.5ρ*norm([Φₓ,Φᵤ]-ν, :𝓗₂));
    @constraint(problem,                Φₓ[1]   .== Ĩ);
    @constraint(problem, [t = 1:(T-1)], Φₓ[t+1] .== A*Φₓ[t] + B₂*Φᵤ[t]);
    @constraint(problem,                   0    .== A*Φₓ[T] + B₂*Φᵤ[T]);

    for t in 1:T
        fix.(Φₓ[t][𝓢ₓ[t][sₓ,cⱼ] .≠ 1], 0.0, force=true);
        fix.(Φᵤ[t][𝓢ᵤ[t][sᵤ,cⱼ] .≠ 1], 0.0, force=true);
    end
    
    optimize!(problem)

    for t in 1:T 
        Φ̃[1][t][sₓ,cⱼ] = 𝓢ₓ[t][sₓ,cⱼ] .* value.(Φₓ[t]);
        Φ̃[2][t][sᵤ,cⱼ] = 𝓢ᵤ[t][sᵤ,cⱼ] .* value.(Φᵤ[t]);
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