module BitPermutations

using Combinatorics

export MBitVector
export PermutationNetwork, BenesNetwork, GRPNetwork
export bitsize, deltaswap, grpswap, invgrpswap
export bitpermute, invbitpermute



"""
    bitsize(T::DataType)
    bitsize(obj::T)

Number of bits in the binary representations of `T`.
"""
bitsize(::Type{T}) where T = 8*sizeof(T)
bitsize(::T) where T = bitsize(T)

abstract type PermutationNetwork{T} end

include("bitarrays.jl")
include("intrinsics.jl")
include("benes.jl")
include("grp.jl")

end # module
