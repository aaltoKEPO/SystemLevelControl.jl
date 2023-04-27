# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed by
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------
"""
    AbstractFeedbackStructure

An abstract type for specifying feedback structures. 
Its main usage is to specialize type annotations of `AbstractGeneralizedPlant`'s.
"""
abstract type AbstractFeedbackStructure end

"""
    OutputFeedback

A wrapper type to specify output-feedback structures.
"""
struct OutputFeedback <: AbstractFeedbackStructure end

"""
    StateFeedback

A wrapper type to specify state-feedback structures.
"""
struct StateFeedback <: AbstractFeedbackStructure end
# --