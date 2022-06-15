"""
    MBitVector{T} <: AbstractArray{Bool,1}

A mutable, statically-sized bit vector encoded in the bits in `T`.

An `MBitVector` should behave similarly to a `BitVector`, but it is slightly faster as 
the size is known at compile-time to be `bitsize(T)`.
"""
mutable struct MBitVector{T} <: AbstractArray{Bool,1}
    chunk::T

    MBitVector(input::T) where T<:Integer = new{T}(input)
    
    function MBitVector{T}(input::A) where {T<:Integer,A}
        # Input is the data chunk
        isprimitivetype(A) && isbits(input) && return new{T}(convert(T, input))
        # Iterate over input
        a = zero(T)
        n = 0
        # Use bitrotate to cycle through input elements
        for x in input
            convert(Bool, x) && (a |= one(T))
            a = bitrotate(a, -1)
            n += 1
        end
        n ≤ bitsize(T) || throw(OverflowError("Input type $A does not fit in MBitVector{$T}"))
        a = bitrotate(a, n)
        return new{T}(a)
    end
end

# Getter
@inline chunk(m::MBitVector) = m.chunk

# AbstractArray interface
Base.size(::MBitVector{T}) where T = (bitsize(T),)

function Base.getindex(m::MBitVector{T}, i::Int) where T
    @boundscheck checkbounds(m, i)
    a = chunk(m)
    u = one(T) << (i - 1)
    return !iszero(a & u)
end

function Base.setindex!(m::MBitVector{T}, x, i::Int) where T
    @boundscheck checkbounds(m, i)
    a = chunk(m)
    u = one(T) << (i - 1)
    m.chunk = ifelse(convert(Bool, x), a | u, a & ~u)
    return m
end

Base.IndexStyle(::Type{<:MBitVector}) = IndexLinear()

# Comparison
Base.isequal(m₁::MBitVector{T}, m₂::MBitVector{T}) where T = chunk(m₁) == chunk(m₂)

# Initializers
Base.zero(m::MBitVector{T}) where T = MBitVector{T}(zero(T))
Base.similar(m::MBitVector{T}) where T = zero(m)

# Fast Base array operations
Base.any(m::MBitVector) = !iszero(chunk(m))
Base.all(m::MBitVector) = chunk(m) === ~zero(T)
Base.sum(m::MBitVector) = count_ones(chunk(m))

# Conversions
Base.convert(::Type{T}, m::MBitVector) where T<:Integer = convert(T, chunk(m))

# Fast iteration
function Base.iterate(m::MBitVector)
    val, rest = _peel(chunk(m))
    return val, (rest, 1)
end

function Base.iterate(m::MBitVector{T}, state::Tuple{T,Int}) where T
    rest, n = state
    n === length(m) && return nothing
    val, rest = _peel(rest)
    return val, (rest, n+1)
end

_peel(a::Integer) = (Bool(a & one(a)), a >> 1)

# Logical operations
Base.:~(m::MBitVector{T}) where T = MBitVector{T}(~chunk(m))

Base.:&(m₁::MBitVector{T}, m₂::MBitVector{T}) where T = MBitVector{T}(chunk(m₁) & chunk(m₂))

Base.:|(m₁::MBitVector{T}, m₂::MBitVector{T}) where T = MBitVector{T}(chunk(m₁) | chunk(m₂))

Base.xor(m₁::MBitVector{T}, m₂::MBitVector{T}) where T = MBitVector{T}(xor(chunk(m₁), chunk(m₂)))

Base.circshift(m::MBitVector{T}, n::Int) where T = MBitVector{T}(bitrotate(chunk(m), n))

