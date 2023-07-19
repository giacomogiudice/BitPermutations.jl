using Base: llvmcall

"""
    has_bmi2()

Tests if the processor supports hardware BMI2 instructions.
"""
function has_bmi2()
    Sys.ARCH ≡ :x86_64 || return false

    # Explanation:
    # https://stackoverflow.com/questions/32214843/compiler-macro-to-detect-bmi2-instruction-set
    CPUInfo = zeros(Int32, 4)
    ccall(:jl_cpuidex, Cvoid, (Ptr{Cint}, Cint, Cint), CPUInfo, 7, 0)
    return !iszero(CPUInfo[2] & 0x100)
end

"""
    use_bmi2()

Used to set `USE_BMI2` flag at precompilation.
"""
function use_bmi2()
    flag = get(ENV, "BP_USE_BMI2", true)
    return (flag isa Bool ? flag : parse(Bool, flag)) && has_bmi2()
end

const USE_BMI2 = use_bmi2()

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
Dictionary for the lookup of the LLVM type for a Julia type.
"""
const LLVM_TYPES = Dict(UInt8 => "i8", UInt16 => "i16", UInt32 => "i32", UInt64 => "i64")

"""
    _avx_bit_shuffle(x::T, m::Vec{W,UInt8})

Perform a bit-shuffle using AVX instructions.
"""
@generated function _avx_bit_shuffle(x::T, m::Vec{W,UInt8}) where {T<:Unsigned,W}
    @assert W == 8 * sizeof(T)
    N = 8 * W
    R = LLVM_TYPES[T]
    # Use `vpshufbitqmb` intrisic to perform the permutation
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
        llvmcall(($decl, $instr), $T, Tuple{Vec{$W,UInt8},Vec{$W,UInt8}}, rep, m)
    end
end
