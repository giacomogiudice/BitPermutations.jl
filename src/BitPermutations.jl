module BitPermutations

using Combinatorics
using Base: Broadcast

export PermutationNetwork, BenesNetwork, GRPNetwork
export AbstractBitPermutation, BitPermutation, AdjointBitPermutation

export bitsize
export bitpermute, invbitpermute
export cycles, order

export CompiledBitPermutation
export @bitpermutation

include("intrinsics.jl")
include("primitives.jl")
include("bitarrays.jl")
include("networks.jl")
include("permutations.jl")
include("compiled.jl")

end # module
