module BitPermutations

using Combinatorics
using Base: Broadcast

export Bits
export PermutationNetwork, BenesNetwork, GRPNetwork
export AbstractBitPermutation, BitPermutation, AdjointBitPermutation

export bitsize
export bitpermute, invbitpermute, bitpermute!, invbitpermute!
export cycles, order

include("intrinsics.jl")
include("primitives.jl")
include("bitarrays.jl")
include("networks.jl")
include("permutations.jl")

end # module
