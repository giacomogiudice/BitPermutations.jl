module BitPermutations

using Combinatorics

export MBitVector
export BenesNetwork
export bitsize, deltaswap
export bitpermute, invbitpermute

bitsize(::Type{T}) where T = 8*sizeof(T)
bitsize(::T) where T = bitsize(T)

include("bitarrays.jl")
include("benes.jl")

end # module
