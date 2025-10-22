module BitPermutations

using Combinatorics
using Preferences
using SIMD

export Bits
export PermutationBackend, BenesNetwork, GRPNetwork, AVXCopyGather
export AbstractBitPermutation, BitPermutation, AdjointBitPermutation

export bitsize
export bitpermute, invbitpermute, bitpermute!, invbitpermute!
export cycles, order

include("primitives.jl")
include("bitarrays.jl")
include("backends.jl")
include("permutations.jl")

end # module
