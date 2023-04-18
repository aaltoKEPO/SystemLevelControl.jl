# -----------------------------------------------------------------------
## Copyright (C) 2023- by Otacilio 'Minho' Neto, <otacilio.neto@aalto.fi>
# This code is part of the 'SystemLevelControl.jl' package, licensed
# the MIT License (see <https://spdx.org/licenses/MIT.html> )                
# -----------------------------------------------------------------------

# GENERALIZED PLANT CONVERSIONS _________________________________________ 

# AUXILIARY CONVERSIONS _________________________________________________
# The following code is adapted from ControlSystems.jl (Copyright (c) 2014-2018: Jim Crist, Mattias Fält, Fredrik Bagge Carlson and other contributors)
to_sparse_matrix(T, A::Number) = convert(SparseMatrixCSC{T,Int}, fill(A,1,1))
to_sparse_matrix(T, A::AbstractVector) = convert(SparseMatrixCSC{T,Int}, reshape(A,:,1))
to_sparse_matrix(T, A::AbstractMatrix) = convert(SparseMatrixCSC{T,Int}, A)

fix_feedthrough(D::AbstractMatrix{T}, ny::Integer, nu::Integer) where {T} = (D == 0I) ? spzeros(T, ny, nu) : D;
