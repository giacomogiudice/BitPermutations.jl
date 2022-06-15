module BitPermutations

using Combinatorics

export MBitVector
export BenesNetwork, GRPNetwork
export bitsize, deltaswap, grpswap, invgrpswap
export bitpermute, invbitpermute


"""
    bitsize(T::DataType)
    bitsize(obj::T)

Number of bits in the binary representations of `T`.
"""
bitsize(::Type{T}) where T = 8*sizeof(T)
bitsize(::T) where T = bitsize(T)

include("bitarrays.jl")
include("benes.jl")
include("intrinsics.jl")
include("grp.jl")

end # module
