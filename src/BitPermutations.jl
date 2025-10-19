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

function __init__()
    @has_preference("use_bmi2") || @set_preferences!("use_bmi2" => use_bmi2())
    @has_preference("use_avx512") || @set_preferences!("use_avx512" => use_avx512())
    return nothing
end

end # module
