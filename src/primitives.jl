"""
    bitsize(T::DataType)
    bitsize(obj::T)

Number of bits in the binary representations of `T`.
"""
bitsize(::Type{T}) where T<:Base.BitInteger = 8*sizeof(T)
bitsize(::T) where T = bitsize(T)

"""
    deltaswap(x::T, m::T, shift::Int)
    deltaswap(x::AbstractVector, m::AbstractVector{Bool}, shift::Int)

Swaps bits in `x` selected by mask `m` with ones to the left by an amount specifield 
by `shift`. The `AbstractArray` version is not optimized to be fast.
"""
@inline function deltaswap(x::T, m::T, shift::Int) where T<:Base.BitInteger
    t = ((x >> shift) ⊻ x) & m
    return x ⊻ t ⊻ (t << shift)
end

function deltaswap(x::AbstractVector, m::AbstractVector{Bool}, shift::Int)
    @assert length(x) === length(m)
    y = copy(x)
    for (i, mᵢ) in enumerate(m)
        if !iszero(mᵢ)
            y[i+shift], y[i] = x[i], x[i+shift]
        end
    end
    return y
end

"""
    grpswap(x::T, m::T, [shift::Int], [m̄::T])
    grpswap(x::AbstractVector, m::AbstractVector{Bool})

Moves the bits in `x` selected by the mask `m` to the left, the rest get moved to the right.
The `shift` (number of 0s in m) and inverse mask `m̄` can be provided so they don't have to be recomputed.
The `AbstractArray` version is not optimized to be fast.

See also [`invgrpswap`](@ref).
"""
@inline function grpswap(x::T, m::T, shift::Int=count_zeros(m), m̄::T=~m) where T<:Base.BitInteger
    return pext(x, m) << shift | pext(x, m̄)
end

function grpswap(x::AbstractVector, m::AbstractVector{Bool})
    return vcat(x[.~m], x[m])
end

"""
    invgrpswap(x::T, m::T, [shift::Int], [m̄::T])
    invgrpswap(x::AbstractVector, m::AbstractVector{Bool})

Performs the inverse operation of `grpswap`.

See also [`grpswap`](@ref).
"""
@inline function invgrpswap(x::T, m::T, shift::Int=count_zeros(m), m̄::T=~m) where T<:Base.BitInteger
    return pdep(x >> shift, m) | pdep(x, m̄)
end

function invgrpswap(x::AbstractVector, m::AbstractVector{Bool})
    s = sum(.~m)
    y = similar(x)
    y[.~m] = x[begin:(begin+s-1)]
    y[m] = x[(begin+s):end]
    return y
end
