"""
    has_bmi2()

Tests if the processor supports hardware BMI2 instructions.
"""
function has_bmi2()
    CPUInfo = zeros(Int32, 4)
    # Explanation:
    # https://stackoverflow.com/questions/32214843/compiler-macro-to-detect-bmi2-instruction-set
    if Sys.ARCH === :x86_64
        try
            ccall(:jl_cpuidex, Cvoid, (Ptr{Cint}, Cint, Cint), CPUInfo, 7, 0)
            return !iszero(CPUInfo[2] & 0x100)
        catch err
            @warn "Failed to query CPU information. $(err.msg)"
        end
    end
    return false
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
pdep(x::T, m::T) where T<:Union{UInt32,UInt64} = USE_BMI2 ? _pdep(x, m) : _pdep_fallback(x, m)
pdep(x::T, m::T) where T = _pdep_fallback(x, m)

"""
    pext(x::T, m::T)

Select bits of `x` with `m` and place them at the beginning of the return value.
Opposite operation of `pdep`.

See also [`pdep`](@ref).
"""
pext(x::T, m::T) where T<:Union{UInt32,UInt64} = USE_BMI2 ? _pext(x, m) : _pext_fallback(x, m)
pext(x::T, m::T) where T = _pext_fallback(x, m)

@inline _pdep(x::UInt32, m::UInt32) = ccall("llvm.x86.bmi.pdep.32", llvmcall, UInt32, (UInt32, UInt32), x, m)
@inline _pdep(x::UInt64, m::UInt64) = ccall("llvm.x86.bmi.pdep.64", llvmcall, UInt64, (UInt64, UInt64), x, m)

@inline _pext(x::UInt32, m::UInt32) = ccall("llvm.x86.bmi.pext.32", llvmcall, UInt32, (UInt32, UInt32), x, m)
@inline _pext(x::UInt64, m::UInt64) = ccall("llvm.x86.bmi.pext.64", llvmcall, UInt64, (UInt64, UInt64), x, m)

# Fallbacks adapted from https://www.chessprogramming.org/BMI2
function _pdep_fallback(x::T, mask::T) where T
    r = zero(T)
    bb = one(T)
    while !iszero(mask)
        !iszero(x & bb) && (r |= mask & -mask)
        mask &= mask - one(T)
        bb += bb
    end
    return r
end

function _pext_fallback(x::T, mask::T) where T
    r = zero(T)
    bb = one(T)
    while !iszero(mask)
        !iszero(x & mask & -mask) && (r |= bb)
        mask &= mask - one(T)
        bb += bb
    end
    return r
end
