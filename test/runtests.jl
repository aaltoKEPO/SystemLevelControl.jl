# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
using SystemLevelControl
using Test, SafeTestsets

@time begin 
    @time @safetestset "(Types) GeneralizedPlant composite type" begin include("types_GeneralizedPlant_test.jl") end
    @time @safetestset "(Types) Operations on AbstractGeneralizedPlants" begin include("types_operations_test.jl") end
    #
    @time @safetestset "Control synthesis methods" begin include("synthesis_test.jl") end
    #
    @time @safetestset "dimensionality reduction methods" begin include("reduction_test.jl") end
end
