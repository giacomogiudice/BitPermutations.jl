module BitPermutations

using Combinatorics

export MBitVector
export BitPermutationAlgorithm, BenesNetwork, GRPNetwork
export AbstractBitPermutation, BitPermutation, AdjointBitPermutation, CompiledBitPermutation

export bitsize, deltaswap, grpswap, invgrpswap
export bitpermute, invbitpermute

export @bitpermutation

include("intrinsics.jl")
include("primitives.jl")
include("bitarrays.jl")
include("networks.jl")
include("permutations.jl")
include("macros.jl")

end # module
