# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
## This file defines some auxiliary functions and data-types mainly
##  to make code less verbose, or slightly more idiomatic

function generateTree(E, S=()->rand(Uniform(0.5,1)))
    N = size(E,1); Vᵥ = spzeros(N); Eₜ = spzeros(N,N); 
    Vᵢ = rand(1:N);
    while sum(Vᵥ) < N
        Vⱼ = rand( findall( E[:,Vᵢ] .== 1 ) );
        if Vᵥ[Vⱼ] == 0 
            Eₜ[Vⱼ,Vᵢ] = S();
            Vᵥ[Vⱼ] = 1; 
        end 
        Vᵢ = Vⱼ;
    end
    return 0.5(Eₜ + Eₜ')
end