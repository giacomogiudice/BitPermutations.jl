"""
    Bits{T} <: AbstractArray{Bool,1}

A mutable, statically-sized bit vector encoded in the bits in an element of type
`T <: Unsigned`.

An `Bits` should behave similarly to a `BitVector`, but it is slightly faster as 
the size is known at compile-time to be `bitsize(T)`.
"""
mutable struct Bits{T} <: AbstractArray{Bool,1}
    chunk::T

    # Constructor for Integers
    Bits{T}(input::T) where T<:Unsigned = new{T}(input)
    
    # Constructor for iterable inputs
    function Bits{T}(input) where T<:Unsigned
        # Iterate over input
        val = zero(T)
        mask = one(T)
        n = 0
        # Set nth bit using mask 
        for x in input
            convert(Bool, x) && (val |= mask)
            mask <<= 1
            n += 1
        end
        n ≤ bitsize(T) || throw(OverflowError("Input $input does not fit in Bits{$T}"))
        return new{T}(val)
    end
end

"""
    Bits(input::T)
    Bits{T}(input::T)

Construct an `Bits` from the bits in `input`.
"""
Bits(input::T) where T = Bits{T}(input)

"""
    Bits{T}(input)

Construct an `Bits` encoding the bits of some array, generator, or iterable in the bits of type `T`.
""" 
Bits(input::Bits{T}) where T = Bits{T}(chunk(input))

# Getter
@inline chunk(m::Bits) = m.chunk

# AbstractArray interface
Base.size(::Bits{T}) where T = (bitsize(T),)

function Base.getindex(m::Bits{T}, i::Int) where T
    @boundscheck checkbounds(m, i)
    a = chunk(m)
    u = one(T) << (i - 1)
    return !iszero(a & u)
end

function Base.setindex!(m::Bits{T}, x, i::Int) where T
    @boundscheck checkbounds(m, i)
    a = chunk(m)
    u = one(T) << (i - 1)
    m.chunk = ifelse(convert(Bool, x), a | u, a & ~u)
    return m
end

Base.IndexStyle(::Type{<:Bits}) = IndexLinear()

# Comparison
Base.isequal(m₁::Bits{T}, m₂::Bits{T}) where T = m₁ == m₂
Base.:(==)(m₁::Bits{T}, m₂::Bits{T}) where T = chunk(m₁) == chunk(m₂)
Base.isless(m₁::Bits{T}, m₂::Bits{T}) where T = chunk(m₁) < chunk(m₂)

# Initializers
Base.zero(m::Bits{T}) where T = Bits{T}(zero(T))
Base.similar(m::Bits{T}) where T = zero(m)

# Fast Base array operations
Base.any(m::Bits) = !iszero(chunk(m))
Base.all(m::Bits) = chunk(m) === ~zero(chunk(m))
Base.sum(m::Bits) = count_ones(chunk(m))

# Conversions
Base.convert(::Type{T}, m::Bits) where T<:Integer = convert(T, chunk(m))

# Fast iteration
function Base.iterate(m::Bits)
    val, rest = _peel(chunk(m))
    return val, (rest, 1)
end

function Base.iterate(m::Bits{T}, state::Tuple{T,Int}) where T
    rest, n = state
    n === length(m) && return nothing
    val, rest = _peel(rest)
    return val, (rest, n+1)
end

_peel(a::Integer) = (Bool(a & one(a)), a >> 1)

# Specialize broadcasting of logical operations,
# arbitrary ones are harder because of variable length issue
Base.broadcasted(::Broadcast.DefaultArrayStyle{1}, ::typeof(~), m::Bits) = Bits(~chunk(m))
Base.broadcasted(bs::Broadcast.DefaultArrayStyle{1}, ::typeof(!), m::Bits) = Base.broadcasted(bs, ~, m)

for op in (:&, :|, :xor)
    @eval function Base.broadcasted(::Broadcast.DefaultArrayStyle{1}, ::typeof($op), m₁::Bits{T}, m₂::Bits{T}) where T 
        return Bits($op(chunk(m₁), chunk(m₂)))
    end
end

# Fast circshift 
Base.circshift(m::Bits, n::Int) = Bits(bitrotate(chunk(m), n))
