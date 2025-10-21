using Base.BinaryPlatforms: CPUID
using Base: llvmcall

# const ENV_KEY = "BIT_PERMUTATIONS_USE_INTRINSICS"

# """
#     use_bmi2()

# Used to set `USE_BMI2` flag at precompilation.
# """
# function use_bmi2()
#     haskey(ENV, ENV_KEY) && return parse(Bool, ENV[ENV_KEY])
#     return Sys.ARCH == :x86_64 && CPUID.test_cpu_feature(CPUID.JL_X86_bmi2)
# end

# """
#     use_avx512()

# Used to set `USE_AVX512` flag at precompilation.
# """
# function use_avx512()
#     haskey(ENV, ENV_KEY) && return parse(Bool, ENV[ENV_KEY])
#     return Sys.ARCH == :x86_64 && CPUID.test_cpu_feature(CPUID.JL_X86_avx512bitalg)
# end

# const USE_BMI2 = false
# const USE_AVX512 = false

"""
    bitsize(::Type{T})
    bitsize(obj::T)

Number of bits in the binary representations of any primitive type `T`.
"""
function bitsize(x::Type{T}) where {T}
    isprimitivetype(x) || throw(ArgumentError("Argument of `bitsize` must be a primitive type"))
    return 8 * sizeof(x)
end

bitsize(x::T) where {T} = bitsize(T)

# """
#     shift_safe(::Type{T}, s::Integer)

# Changes the shifting amount for bitshifts to guarantee to the compiler that the shift amount will
# not exceed `bitsize(T)`.
# Because bit shifting behaves slighly differently in Julia vs. LLVM, this help the compiler emit
# less code.

# See also: https://github.com/JuliaLang/julia/issues/30674.
# """
# @inline shift_safe(::Type{T}, s::Integer) where {T} = s & (bitsize(T) - 1)

# """
#     deltaswap(x::T, m::T, shift::Integer)
#     deltaswap(x::AbstractVector, m::AbstractVector{Bool}, shift::Integer)

# Swaps bits in `x` selected by mask `m` with ones to the left by an amount specifield
# by `shift`. The `AbstractArray` version is not optimized to be fast.
# """
# function deltaswap(x::T, m::T, shift::Integer) where {T<:Unsigned}
#     shift = shift_safe(T, shift)
#     t = ((x >> shift) ⊻ x) & m
#     return x ⊻ t ⊻ (t << shift)
# end

# function deltaswap(x::AbstractVector, m::AbstractVector{Bool}, shift::Integer)
#     @assert length(x) === length(m)
#     y = copy(x)
#     for (i, mᵢ) in enumerate(m)
#         if !iszero(mᵢ)
#             y[i + shift], y[i] = x[i], x[i + shift]
#         end
#     end
#     return y
# end

# """
#     pdep(x::T, m::T)

# Place first bits of `x` at locations specified by `m` in the return value.
# Opposite operation of `pdep`.

# See also [`pext`](@ref).
# """
# pdep(x::T, m::T) where {T<:Union{UInt32,UInt64}} = USE_BMI2 ? _pdep(x, m) : _pdep_fallback(x, m)
# pdep(x::T, m::T) where {T} = _pdep_fallback(x, m)

# """
#     pext(x::T, m::T)

# Select bits of `x` with `m` and place them at the beginning of the return value.
# Opposite operation of `pdep`.

# See also [`pdep`](@ref).
# """
# pext(x::T, m::T) where {T<:Union{UInt32,UInt64}} = USE_BMI2 ? _pext(x, m) : _pext_fallback(x, m)
# pext(x::T, m::T) where {T} = _pext_fallback(x, m)

# @inline _pdep(x::UInt32, m::UInt32) = ccall("llvm.x86.bmi.pdep.32", llvmcall, UInt32, (UInt32, UInt32), x, m)
# @inline _pdep(x::UInt64, m::UInt64) = ccall("llvm.x86.bmi.pdep.64", llvmcall, UInt64, (UInt64, UInt64), x, m)

# @inline _pext(x::UInt32, m::UInt32) = ccall("llvm.x86.bmi.pext.32", llvmcall, UInt32, (UInt32, UInt32), x, m)
# @inline _pext(x::UInt64, m::UInt64) = ccall("llvm.x86.bmi.pext.64", llvmcall, UInt64, (UInt64, UInt64), x, m)

# # Fallbacks adapted from https://www.chessprogramming.org/BMI2
# function _pdep_fallback(x::T, mask::T) where {T}
#     r = zero(T)
#     bb = one(T)
#     while !iszero(mask)
#         !iszero(x & bb) && (r |= mask & -mask)
#         mask &= mask - one(T)
#         bb += bb
#     end
#     return r
# end

# function _pext_fallback(x::T, mask::T) where {T}
#     r = zero(T)
#     bb = one(T)
#     while !iszero(mask)
#         !iszero(x & mask & -mask) && (r |= bb)
#         mask &= mask - one(T)
#         bb += bb
#     end
#     return r
# end

# """
#     grpswap(x::T, m::T, [shift::Integer], [m̄::T])
#     grpswap(x::AbstractVector, m::AbstractVector{Bool})

# Moves the bits in `x` selected by the mask `m` to the left, the rest get moved to the right.
# The `shift` (number of 0s in m) and inverse mask `m̄` can be provided so they don't have to be
# recomputed.
# The `AbstractArray` version is not optimized to be fast.

# See also [`invgrpswap`](@ref).
# """
# function grpswap(x::T, m::T, shift::Integer=count_zeros(m), m̄::T=(~m)) where {T<:Unsigned}
#     shift = shift_safe(T, shift)
#     return pext(x, m) << shift | pext(x, m̄)
# end

# function grpswap(x::AbstractVector, m::AbstractVector{Bool})
#     @assert length(x) === length(m)
#     return vcat(x[.~m], x[m])
# end

# """
#     invgrpswap(x::T, m::T, [shift::Integer], [m̄::T])
#     invgrpswap(x::AbstractVector, m::AbstractVector{Bool})

# Performs the inverse operation of `grpswap`.

# See also [`grpswap`](@ref).
# """
# function invgrpswap(x::T, m::T, shift::Integer=count_zeros(m), m̄::T=(~m)) where {T<:Unsigned}
#     shift = shift_safe(T, shift)
#     return pdep(x >> shift, m) | pdep(x, m̄)
# end

# function invgrpswap(x::AbstractVector, m::AbstractVector{Bool})
#     @assert length(x) === length(m)
#     s = sum(!, m)
#     y = similar(x)
#     y[.~m] = x[begin:(begin + s - 1)]
#     y[m] = x[(begin + s):end]
#     return y
# end

# """
#     avx_bit_shuffle(x::T, mask::Vec{W,UInt8})

# Perform a bit-shuffle of the bits in `x` using AVX instructions.
# The mask `mask` stores the zero-based indices of the locations of each bit in the permuted vector.
# """
# function avx_bit_shuffle(x::T, m::Vec{W,UInt8}) where {T<:Union{UInt16,UInt32,UInt64},W}
#     return USE_AVX512 ? _avx_bit_shuffle(x, m) : _avx_bit_shuffle_fallback(x, m)
# end

# avx_bit_shuffle(x::T, m::Vec{W,UInt8}) where {T<:Unsigned,W} = _avx_bit_shuffle_fallback(x, m)

# const LLVM_TYPES = Dict(UInt8 => "i8", UInt16 => "i16", UInt32 => "i32", UInt64 => "i64")

# @generated function _avx_bit_shuffle(x::T, mask::Vec{W,UInt8}) where {T<:Union{UInt16,UInt32,UInt64},W}
#     @assert USE_AVX512 || error("AVX512 instructions are not supported on this processor.")
#     # `W` is the size of the each block, `N` is the total number of bits in the SIMD register
#     @assert W == bitsize(T)
#     N = 8 * W
#     R = LLVM_TYPES[T]
#     # Use `vpshufbitqmb` intrisic to perform the permutation, adapted from 
#     # https://lemire.me/blog/2023/06/29/dynamic-bit-shuffle-using-avx-512/
#     decl = """
#         define $R @bit_shuffle(<$W x i8> %a, <$W x i8> %b) {
#             %mask = call <$W x i1> @llvm.x86.avx512.vpshufbitqmb.$N(<$W x i8> %a, <$W x i8> %b)
#             %res = bitcast <$W x i1> %mask to $R
#             ret $R %res
#        }

#        declare <$W x i1> @llvm.x86.avx512.vpshufbitqmb.$N(<$W x i8>, <$W x i8>)
#     """
#     instr = "bit_shuffle"

#     return quote
#         # Copy value of `x` 8 times to fill SIMD register, reinterpret it to match LLVM definition
#         rep = reinterpret(Vec{$W,UInt8}, Vec{8,$T}(x))
#         # Call intrinsic via LLVM API
#         llvmcall(($decl, $instr), $T, Tuple{Vec{$W,UInt8},Vec{$W,UInt8}}, rep, mask)
#     end
# end

# # Fallback: perform slow permutation by checking the bit at each mask location and setting it on the output
# function _avx_bit_shuffle_fallback(x::T, mask::Vec{W,UInt8}) where {T<:Unsigned,W}
#     @assert W == bitsize(T)
#     u = one(T)
#     out = zero(T)
#     @inbounds for i in 1:W
#         bit = !iszero(x & u << mask[i])
#         out |= T(bit) << shift_safe(T, i - 1)
#     end
#     return out
# end
