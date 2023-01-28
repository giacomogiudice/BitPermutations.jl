"""
    bitsize(T::Unsigned)
    bitsize(obj::T)

Number of bits in the binary representations of `T`.
"""
bitsize(::Type{T}) where {T<:Unsigned} = 8 * sizeof(T)
bitsize(::T) where {T<:Unsigned} = bitsize(T)

"""
    mask_shift(::Type{T}, s::Integer)

Changes the mask to guarantee it doesn't shift beyond `bitsize(T)`.
Because bit shifting behaves slighly differently in Julia vs. LLVM, this help the compiler emit
less code.

See also: https://github.com/JuliaLang/julia/issues/30674.
"""
@inline mask_shift(::Type{T}, s::Integer) where {T} = s & (bitsize(T) - 1)

"""
    deltaswap(x::T, m::T, shift::Integer)
    deltaswap(x::AbstractVector, m::AbstractVector{Bool}, shift::Integer)

Swaps bits in `x` selected by mask `m` with ones to the left by an amount specifield
by `shift`. The `AbstractArray` version is not optimized to be fast.
"""
function deltaswap(x::T, m::T, shift::Integer) where {T<:Unsigned}
    shift_masked = mask_shift(T, shift)
    t = ((x >> shift_masked) ⊻ x) & m
    return x ⊻ t ⊻ (t << shift_masked)
end

function deltaswap(x::AbstractVector, m::AbstractVector{Bool}, shift::Integer)
    @assert length(x) === length(m)
    y = copy(x)
    for (i, mᵢ) in enumerate(m)
        if !iszero(mᵢ)
            y[i + shift], y[i] = x[i], x[i + shift]
        end
    end
    return y
end

"""
    grpswap(x::T, m::T, [shift::Integer], [m̄::T])
    grpswap(x::AbstractVector, m::AbstractVector{Bool})

Moves the bits in `x` selected by the mask `m` to the left, the rest get moved to the right.
The `shift` (number of 0s in m) and inverse mask `m̄` can be provided so they don't have to be
recomputed.
The `AbstractArray` version is not optimized to be fast.

See also [`invgrpswap`](@ref).
"""
function grpswap(x::T, m::T, shift::Integer=count_zeros(m), m̄::T=~m) where {T<:Unsigned}
    masked_shift = mask_shift(T, shift)
    return pext(x, m) << masked_shift | pext(x, m̄)
end

function grpswap(x::AbstractVector, m::AbstractVector{Bool})
    @assert length(x) === length(m)
    return vcat(x[.~m], x[m])
end

"""
    invgrpswap(x::T, m::T, [shift::Integer], [m̄::T])
    invgrpswap(x::AbstractVector, m::AbstractVector{Bool})

Performs the inverse operation of `grpswap`.

See also [`grpswap`](@ref).
"""
function invgrpswap(x::T, m::T, shift::Integer=count_zeros(m), m̄::T=~m) where {T<:Unsigned}
    masked_shift = mask_shift(T, shift)
    return pdep(x >> masked_shift, m) | pdep(x, m̄)
end

function invgrpswap(x::AbstractVector, m::AbstractVector{Bool})
    @assert length(x) === length(m)
    s = sum(!, m)
    y = similar(x)
    y[.~m] = x[begin:(begin + s - 1)]
    y[m] = x[(begin + s):end]
    return y
end
