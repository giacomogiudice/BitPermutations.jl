module BitPermutations

using Combinatorics
using SIMD

export Bits
export PermutationNetwork, BenesNetwork, GRPNetwork, AVXCopyGather
export AbstractBitPermutation, BitPermutation, AdjointBitPermutation

export bitsize
export bitpermute, invbitpermute, bitpermute!, invbitpermute!
export cycles, order

include("primitives.jl")
include("bitarrays.jl")
include("networks.jl")
include("permutations.jl")

end # module
