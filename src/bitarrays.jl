"""
    Bits{T} <: AbstractArray{Bool,1}

A mutable, statically-sized bit vector encoded in the bits in an element of type `T <: Integer`.

A `Bits` should behave similarly to a `BitVector`, but it is slightly faster as the size is known
at compile-time to be `bitsize(T)`.
"""
mutable struct Bits{T} <: AbstractArray{Bool,1}
    chunk::T

    # Constructor for Integers
    Bits{T}(input::Integer) where {T<:Integer} = new{T}(convert(T, input))

    # Constructor for iterable inputs
    function Bits{T}(input) where {T<:Integer}
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
    Bits(input::Integer)
    Bits{T}(input::Integer)
    Bits{T}(input)

Construct a `Bits` from the bits in `input`, passed either as a `Integer` or from an array, or any
other generic iterable.
It is basically a view into the bits of `input`, which can then be manipulated using vector-like
syntax.
"""
Bits(input::T) where {T} = Bits{T}(input)
Bits(input::Bits{T}) where {T} = Bits{T}(chunk(input))

# Getter
@inline chunk(m::Bits) = m.chunk

# AbstractArray interface
Base.size(::Bits{T}) where {T} = (bitsize(T),)

function Base.getindex(m::Bits{T}, i::Int) where {T}
    @boundscheck checkbounds(m, i)
    a = chunk(m)
    u = one(T) << shift_safe(T, i - 1)
    return !iszero(a & u)
end

function Base.setindex!(m::Bits{T}, x, i::Int) where {T}
    @boundscheck checkbounds(m, i)
    a = chunk(m)
    v = convert(T, convert(Bool, x))
    s = shift_safe(T, i - 1)
    m.chunk = (a & ~(one(T) << s)) | (v << s)
    return m
end

Base.IndexStyle(::Type{<:Bits}) = IndexLinear()

# Comparison
Base.isequal(m₁::Bits{T}, m₂::Bits{T}) where {T} = m₁ == m₂
Base.:(==)(m₁::Bits{T}, m₂::Bits{T}) where {T} = chunk(m₁) == chunk(m₂)
Base.isless(m₁::Bits{T}, m₂::Bits{T}) where {T} = chunk(m₁) < chunk(m₂)

# Initializers
Base.zero(m::Bits{T}) where {T} = Bits{T}(zero(T))
Base.similar(m::Bits{T}) where {T} = zero(m)

# Fast Base array operations
Base.any(m::Bits) = !iszero(chunk(m))
Base.all(m::Bits) = chunk(m) === ~zero(chunk(m))
Base.sum(m::Bits) = count_ones(chunk(m))

# Conversions
Base.convert(::Type{T}, m::Bits) where {T<:Integer} = convert(T, chunk(m))

# Fast iteration
function Base.iterate(m::Bits)
    val, rest = _peel(chunk(m))
    return val, (rest, 1)
end

function Base.iterate(m::Bits{T}, state::Tuple{T,Int}) where {T}
    rest, n = state
    n === length(m) && return nothing
    val, rest = _peel(rest)
    return val, (rest, n + 1)
end

_peel(a::Integer) = (Bool(a & one(a)), a >> 1)

# Specialize broadcasting of logical operations,
# arbitrary ones are harder because of variable length issue
Base.broadcasted(::Broadcast.DefaultArrayStyle{1}, ::typeof(~), m::Bits) = Bits(~chunk(m))
Base.broadcasted(bs::Broadcast.DefaultArrayStyle{1}, ::typeof(!), m::Bits) = Base.broadcasted(bs, ~, m)

for op in (:&, :|, :xor)
    @eval function Base.broadcasted(::Broadcast.DefaultArrayStyle{1}, ::typeof($op), m₁::Bits{T}, m₂::Bits{T}) where {T}
        return Bits($op(chunk(m₁), chunk(m₂)))
    end
end

# Fast circshift 
Base.circshift(m::Bits, n::Int) = Bits(bitrotate(chunk(m), n))
