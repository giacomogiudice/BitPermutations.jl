"""
    AbstractBitPermutation{T,A}

Abstract bit permutation supertype.
"""
abstract type AbstractBitPermutation{T} end

(perm::AbstractBitPermutation)(x) = bitpermute(x, perm)

"""
    bitpermute(x::Number, p::AbstractBitPermutation)
    bitpermute(x::Number, net::PermutationNetwork)
    bitpermute(x::AbstractArray{T}, p::AbstractBitPermutation{T})

Perform a permutation the bits in `x` using an `AbstractBitPermutation`, using the permutation
specified in `p`.
The syntax `p(x)` can be used as well.
Internally, the permutation operations are stored as a `PermutationNetwork`, so the permutation can
be called using the network directly.
If the input is an `AbstractArray` subtype, the permutation is performed element-wise.
This can be called as well with `p.(x)` or `bitpermute.(x, p)`.
For cases in which the input array is not needed afterwards, the in-place version `bitpermute!` is
preferred.

See also [`invbitpermute`](@ref), [`bitpermute!`](@ref).
"""
bitpermute(x::AbstractArray{T}, P::AbstractBitPermutation{T}) where {T} = bitpermute_elementwise(x, P)

"""
    bitpermute!(x::AbstractArray{T}, P::AbstractBitPermutation{T})

Perform the  permutation of the bits in each element of `x` using an `AbstractBitPermutation` in
place, i.e. by mutating the array `x`.

See also [`bitpermute`](@ref).
"""
bitpermute!(x::AbstractArray{T}, P::AbstractBitPermutation{T}) where {T} = bitpermute_elementwise!(x, P)

"""
    invbitpermute(x::Number, p::AbstractBitPermutation)
    invbitpermute(x::Number, net::PermutationNetwork)
    invbitpermute(x::AbstractArray{T}, P::AbstractBitPermutation{T})

Perform the inverse permutation of the bits in `x` using an `AbstractBitPermutation`, using the
inverse of the permutation specified in `p`.
The syntax `p'(x)` can be used as well.
Internally, the permutation operations are stored as a `PermutationNetwork`, so the permutation can
be called using the network directly.
If the input is an `AbstractArray` subtype, the inverse permutation is performed element-wise.
This can be called as well with `p'.(x)` or `invbitpermute.(x, p)`.
For cases in which the input array is not needed afterwards, the in-place version `invbitpermute!`
is preferred.

See also [`bitpermute`](@ref), [`invbitpermute!`](@ref).
"""
invbitpermute(x::AbstractArray{T}, P::AbstractBitPermutation{T}) where {T} = invbitpermute_elementwise(x, P)

"""
    invbitpermute!(x::AbstractArray{T}, P::AbstractBitPermutation{T})

Perform the inverse permutation of the bits in each element of `x` using an
`AbstractBitPermutation` in place, i.e. by mutating the array `x`.

See also [`invbitpermute`](@ref).
"""
invbitpermute!(x::AbstractArray{T}, P::AbstractBitPermutation{T}) where {T} = invbitpermute_elementwise!(x, P)

Base.broadcasted(P::AbstractBitPermutation{T}, x::AbstractArray{T}) where {T} = bitpermute_elementwise(x, P)

function Base.broadcasted(::typeof(bitpermute), x::AbstractArray{T}, P::AbstractBitPermutation{T}) where {T}
    return bitpermute_elementwise(x, P)
end

function Base.broadcasted(::typeof(invbitpermute), x::AbstractArray{T}, P::AbstractBitPermutation{T}) where {T}
    return invbitpermute_elementwise(x, P)
end

"""
    BitPermutation{T,A}

Represents a bit permutation of type `T`, performed with a `PermutationNetwork` of type `A`.
"""
struct BitPermutation{T,P<:PermutationNetwork{T}} <: AbstractBitPermutation{T}
    vector::Vector{Int}
    network::P

    function BitPermutation{T}(perm::Vector{Int}, net::P) where {T,P<:PermutationNetwork{T}}
        return new{T,P}(perm, net)
    end
end

"""
    BitPermutation{T}(p::AbstractVector{Int}; type=[DEFAULT_TYPE], [...])

Construct a `BitPermutation{T}` from permutation vector `p`.
The type of algorithm can be specified by the `type` keyword argument.
If none is provided, `BenesNetwork` is chosen, unless `T<:Union{UInt32,UInt64}` and
BMI2 instructions are supported by the processor, in which case a `GRPNetwork` is chosen.

See also: [`BenesNetwork`](@ref), [`GRPNetwork`](@ref).
"""
function BitPermutation{T}(p::AbstractVector{Int}; type::Type=_default_type(T), kwargs...) where {T}
    perm = collect(p)
    net = (type){T}(perm; kwargs...)
    return BitPermutation{T}(perm, net)
end

_default_type(::Type{T}) where {T<:UInt16} = ifelse(USE_AVX512, AVXCopyGather, BenesNetwork)

function _default_type(::Type{T}) where {T<:Union{UInt32,UInt64}}
    return ifelse(USE_AVX512, AVXCopyGather, ifelse(USE_BMI2, GRPNetwork, BenesNetwork))
end

_default_type(_::Type{T}) where {T} = BenesNetwork

Base.show(io::IO, ::BitPermutation{T}) where {T} = print(io, "BitPermutation{$T}")

function Base.show(io::IO, ::MIME"text/plain", P::BitPermutation{T}) where {T}
    println(io, "BitPermutation{$T} with network $(P.network):")
    foreach(cycles(P)) do cycle
        length(cycle) === 1 && return nothing
        print(io, "(")
        join(io, cycle, " ")
        print(io, ")")
    end
    return nothing
end

bitpermute(x::Number, P::BitPermutation{T}) where {T} = bitpermute(convert(T, x), P.network)
invbitpermute(x::Number, P::BitPermutation{T}) where {T} = invbitpermute(convert(T, x), P.network)

bitpermute_elementwise(x::AbstractArray, P::BitPermutation) = bitpermute_elementwise(x, P.network)
invbitpermute_elementwise(x::AbstractArray, P::BitPermutation) = invbitpermute_elementwise(x, P.network)
bitpermute_elementwise!(x::AbstractArray, P::BitPermutation) = bitpermute_elementwise!(x, P.network)
invbitpermute_elementwise!(x::AbstractArray, P::BitPermutation) = invbitpermute_elementwise!(x, P.network)

Base.Vector(P::BitPermutation) = P.vector

"""
    AdjointBitPermutation{T,P}

Represents the adjoint permutation of type `P`, where the regular and inverse permutation are
swapped.
"""
struct AdjointBitPermutation{T,P<:AbstractBitPermutation{T}} <: AbstractBitPermutation{T}
    parent::P
end

Base.adjoint(perm::AbstractBitPermutation) = AdjointBitPermutation(perm)

bitpermute(x::Number, P::AdjointBitPermutation) = invbitpermute(x, P.parent)
invbitpermute(x::Number, P::AdjointBitPermutation) = bitpermute(x, P.parent)

bitpermute_elementwise(x::AbstractArray, P::AdjointBitPermutation) = invbitpermute_elementwise(x, P.parent)
invbitpermute_elementwise(x::AbstractArray, P::AdjointBitPermutation) = bitpermute_elementwise(x, P.parent)
bitpermute_elementwise!(x::AbstractArray, P::AdjointBitPermutation) = invbitpermute_elementwise!(x, P.parent)
invbitpermute_elementwise!(x::AbstractArray, P::AdjointBitPermutation) = bitpermute_elementwise!(x, P.parent)

Base.Vector(P::AdjointBitPermutation) = invperm(Vector(P.parent))

"""
    Cycles

Iterator over the cycles of a permutation vector.
"""
struct Cycles
    vector::Vector{Int}
end

Base.IteratorEltype(::Type{<:Cycles}) = Base.HasEltype()
Base.IteratorSize(::Type{<:Cycles}) = Base.SizeUnknown()

Base.eltype(::Cycles) = Vector{Int}

Base.iterate(iter::Cycles) = _nextcycle!(trues(length(iter.vector)), iter.vector)
Base.iterate(iter::Cycles, state::BitVector) = _nextcycle!(state, iter.vector)

function _nextcycle!(state::BitVector, p::AbstractVector{Int})
    # Find first not visited index
    i = findfirst(state)
    isnothing(i) && return nothing

    state[i] = false
    cycle = [i]
    # Iterate until you get back to i
    j = p[i]
    while j ≠ i
        push!(cycle, j)
        state[j] = false
        j = p[j]
    end
    return cycle, state
end

"""
    cycles(P::AbstractBitPermutation)

Returns an iterator over the cycles of the permutation `P`.
"""
cycles(P::AbstractBitPermutation) = Cycles(Vector(P))

"""
    order(P::AbstractBitPermutation)

Return the order of the permutation `P`, i.e. the amount of times you have to apply the permutation
to obtain the trivial permutation.

See also [`isodd`](@ref), [`iseven`](@ref) to compute the parity.
"""
order(P::AbstractBitPermutation) = mapreduce(length, lcm, cycles(P); init=1)

"""
    sign(P::AbstractBitPermutation)

Returns the sign of the permutation `P`, +1 (-1) if it is of even (odd) parity.

See also [`isodd`](@ref), [`iseven`](@ref).
"""
Base.sign(P::AbstractBitPermutation) = ifelse(isodd(P), -1, +1)

"""
    isodd(P::AbstractBitPermutation)

Returns a `Bool` corresponding to whether or not the parity of the permutation is odd.

See also [`iseven`](@ref).
"""
function Base.isodd(P::BitPermutation{T,BenesNetwork{T}}) where {T}
    return isodd(count_ones(mapreduce(first, ⊻, P.network.params; init=zero(T))))
end
function Base.isodd(P::BitPermutation{T,GRPNetwork{T}}) where {T}
    return isodd(count_ones(mapreduce(first, ⊻, P.network.params; init=zero(T))))
end
Base.isodd(P::AdjointBitPermutation) = isodd(P.parent)
Base.isodd(P::AbstractBitPermutation) = mapreduce(isodd ∘ length, ⊻, cycles(P); init=false)

"""
    iseven(P::AbstractBitPermutation)

Returns a `Bool` corresponding to whether or not the parity of the permutation is even.

See also [`isodd`](@ref).
"""
Base.iseven(P::AbstractBitPermutation) = !isodd(P)
