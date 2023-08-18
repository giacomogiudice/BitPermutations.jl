using Base.BinaryPlatforms: CPUID
using Base: llvmcall

const ENV_KEY = "BIT_PERMUTATIONS_USE_INTRINSICS"

"""
    use_bmi2()

Used to set `USE_BMI2` flag at precompilation.
"""
function use_bmi2()
    haskey(ENV, ENV_KEY) && return parse(Bool, ENV[ENV_KEY])
    return Sys.ARCH == :x86_64 && CPUID.test_cpu_feature(CPUID.JL_X86_bmi2)
end

"""
    use_avx512()

Used to set `USE_AVX512` flag at precompilation.
"""
function use_avx512()
    haskey(ENV, ENV_KEY) && return parse(Bool, ENV[ENV_KEY])
    return Sys.ARCH == :x86_64 && CPUID.test_cpu_feature(CPUID.JL_X86_avx512bitalg)
end

const USE_BMI2 = use_bmi2()
const USE_AVX512 = use_avx512()

"""
    pdep(x::T, m::T)

Place first bits of `x` at locations specified by `m` in the return value.
Opposite operation of `pdep`.

See also [`pext`](@ref).
"""
pdep(x::T, m::T) where {T<:Union{UInt32,UInt64}} = USE_BMI2 ? _pdep(x, m) : _pdep_fallback(x, m)
pdep(x::T, m::T) where {T} = _pdep_fallback(x, m)

"""
    pext(x::T, m::T)

Select bits of `x` with `m` and place them at the beginning of the return value.
Opposite operation of `pdep`.

See also [`pdep`](@ref).
"""
pext(x::T, m::T) where {T<:Union{UInt32,UInt64}} = USE_BMI2 ? _pext(x, m) : _pext_fallback(x, m)
pext(x::T, m::T) where {T} = _pext_fallback(x, m)

@inline _pdep(x::UInt32, m::UInt32) = ccall("llvm.x86.bmi.pdep.32", llvmcall, UInt32, (UInt32, UInt32), x, m)
@inline _pdep(x::UInt64, m::UInt64) = ccall("llvm.x86.bmi.pdep.64", llvmcall, UInt64, (UInt64, UInt64), x, m)

@inline _pext(x::UInt32, m::UInt32) = ccall("llvm.x86.bmi.pext.32", llvmcall, UInt32, (UInt32, UInt32), x, m)
@inline _pext(x::UInt64, m::UInt64) = ccall("llvm.x86.bmi.pext.64", llvmcall, UInt64, (UInt64, UInt64), x, m)

# Fallbacks adapted from https://www.chessprogramming.org/BMI2
function _pdep_fallback(x::T, mask::T) where {T}
    r = zero(T)
    bb = one(T)
    while !iszero(mask)
        !iszero(x & bb) && (r |= mask & -mask)
        mask &= mask - one(T)
        bb += bb
    end
    return r
end

function _pext_fallback(x::T, mask::T) where {T}
    r = zero(T)
    bb = one(T)
    while !iszero(mask)
        !iszero(x & mask & -mask) && (r |= bb)
        mask &= mask - one(T)
        bb += bb
    end
    return r
end

"""
    avx_bit_shuffle(x::T, mask::Vec{W,UInt8})

Perform a bit-shuffle of the bits in `x` using AVX instructions.
The mask `mask` stores the zero-based indices of the locations of each bit in the permuted vector.
"""
function avx_bit_shuffle(x::T, m::Vec{W,UInt8}) where {T<:Union{UInt16,UInt32,UInt64},W}
    return USE_AVX512 ? _avx_bit_shuffle(x, m) : _avx_bit_shuffle_fallback(x, m)
end

avx_bit_shuffle(x::T, m::Vec{W,UInt8}) where {T<:Unsigned,W} = _avx_bit_shuffle_fallback(x, m)

const LLVM_TYPES = Dict(UInt8 => "i8", UInt16 => "i16", UInt32 => "i32", UInt64 => "i64")

@generated function _avx_bit_shuffle(x::T, mask::Vec{W,UInt8}) where {T<:Union{UInt16,UInt32,UInt64},W}
    @assert USE_AVX512 || error("AVX512 instructions are not supported on this processor.")
    # `W` is the size of the each block, `N` is the total number of bits
    @assert W == 8 * sizeof(T)
    N = 8 * W
    R = LLVM_TYPES[T]
    # Use `vpshufbitqmb` intrisic to perform the permutation, adapted from 
    # https://lemire.me/blog/2023/06/29/dynamic-bit-shuffle-using-avx-512/
    decl = """
        define $R @bit_shuffle(<$W x i8> %a, <$W x i8> %b) {
            %mask = call <$W x i1> @llvm.x86.avx512.vpshufbitqmb.$N(<$W x i8> %a, <$W x i8> %b)
            %res = bitcast <$W x i1> %mask to $R
            ret $R %res
       }

       declare <$W x i1> @llvm.x86.avx512.vpshufbitqmb.$N(<$W x i8>, <$W x i8>)
    """
    instr = "bit_shuffle"

    return quote
        # Copy value of `x` 8 times to fill SIMD register, reinterpret it to match LLVM definition
        rep = reinterpret(Vec{$W,UInt8}, Vec{8,$T}(x))
        # Call intrinsic via LLVM API
        llvmcall(($decl, $instr), $T, Tuple{Vec{$W,UInt8},Vec{$W,UInt8}}, rep, mask)
    end
end

# Fallback: perform slow permutation by checking the bit at each mask location and setting it on the output
function _avx_bit_shuffle_fallback(x::T, mask::Vec{W,UInt8}) where {T<:Unsigned,W}
    # `W` is the number of bits in `T`
    @assert W == 8 * sizeof(T)
    u = one(T)
    out = zero(T)
    @inbounds for i in 1:W
        bit = !iszero(x & u << mask[i])
        out |= T(bit) << (i - 1)
    end
    return out
end
