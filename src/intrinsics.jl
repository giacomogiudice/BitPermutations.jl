function has_bmi2()
    CPUInfo = zeros(Int32, 4)
    ccall(:jl_cpuidex, Cvoid, (Ptr{Cint}, Cint, Cint), CPUInfo, 7, 0)
    return !iszero(CPUInfo[2] & 0x100)
end

const USE_BMI2 = parse(Bool, get(ENV, "BP_USE_BMI2", "true")) && has_bmi2()

pdep(x::T, m::T) where T<:Union{UInt32,UInt64} = USE_BMI2 ? _pdep(x, m) : _pdep_fallback(x, m)
pdep(x::T, m::T) where T = _pdep_fallback(x, m)

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
