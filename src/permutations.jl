"""
    AbstractBitPermutation{T,A}

Abstract bit permutation supertype.
"""
abstract type AbstractBitPermutation{T} end

(perm::AbstractBitPermutation)(x) = bitpermute(x, perm)

Base.broadcasted(P::AbstractBitPermutation, x) = Base.broadcasted(bitpermute, x, P)

"""
    BitPermutation{T,A}

Represents a bit permutation of type `T`, performed with a `PermutationNetwork` of type `A`.
"""
struct BitPermutation{T,P<:PermutationNetwork{T}} <: AbstractBitPermutation{T}
    vector::Vector{Int}
    network::P

    function BitPermutation{T}(perm::Vector{Int}, net::P) where {T,P<:PermutationNetwork{T}}
        new{T,P}(perm, net)
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
function BitPermutation{T}(p::AbstractVector{Int}; type::Type=_default_type(T), kwargs...) where T
    perm = collect(p)
    net = (type){T}(perm; kwargs...)
    return BitPermutation{T}(perm, net) 
end

_default_type(::Type{T}) where T<:Union{UInt32,UInt64} = ifelse(USE_BMI2, GRPNetwork, BenesNetwork)
_default_type(_::Type{T}) where T = BenesNetwork

Base.show(io::IO, ::BitPermutation{T}) where T = print(io, "BitPermutation{$T}")

function Base.show(io::IO, ::MIME"text/plain", P::BitPermutation{T}) where T
    println(io, "BitPermutation{$T} with network $(P.network):")
    foreach(cycles(P)) do cycle
        length(cycle) === 1 && return
        print(io, "(")
        join(io, cycle, " ")
        print(io, ")")
    end
    return nothing
end

"""
    bitpermute(x, p::AbstractBitPermutation)
    bitpermute(x, n::PermutationNetwork)

Permute bits in `x` using a an `AbstractBitPermutation` or a `PermutationNetwork`.
Callable as well as `p(x)`.

See also [`invbitpermute`](@ref).
"""
@inline bitpermute(x, P::BitPermutation) = bitpermute(x, P.network)

"""
    invbitpermute(x, p::AbstractBitPermutation)
    invbitpermute(x, n::PermutationNetwork)

Permute bits in `x` using a an `AbstractBitPermutation` or a `PermutationNetwork`.
Callable as well as `p'(x)`.

See also [`bitpermute`](@ref).
"""
@inline invbitpermute(x, P::BitPermutation) = invbitpermute(x, P.network)

@inline Base.broadcasted(f::Union{typeof(bitpermute),typeof(invbitpermute)}, x::AbstractArray, P::BitPermutation) = Base.broadcasted(f, x, P.network)

Base.Vector(P::BitPermutation) = P.vector

"""
    AdjointBitPermutation{T,P}

Represents the adjoint permutation of type `P`, where the regular and inverse permutation are swapped.
"""
struct AdjointBitPermutation{T,P<:AbstractBitPermutation{T}} <: AbstractBitPermutation{T}
    parent::P
end

Base.adjoint(perm::AbstractBitPermutation) = AdjointBitPermutation(perm)

@inline bitpermute(x, P::AdjointBitPermutation) = invbitpermute(x, P.parent)
@inline invbitpermute(x, P::AdjointBitPermutation) = bitpermute(x, P.parent)

Base.broadcasted(::typeof(bitpermute), x::AbstractArray, P::AdjointBitPermutation) = Base.broadcasted(invbitpermute, x, P.parent)
Base.broadcasted(::typeof(invbitpermute), x::AbstractArray, P::AdjointBitPermutation) = broadcasted(bitpermute, x, P.parent)

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

Return the order of the permutation `P`, i.e. the amount of times you have to apply the permutation to obtain the trivial permutation.

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
Base.isodd(P::BitPermutation{T,BenesNetwork{T}}) where T = isodd(count_ones(mapreduce(first, ⊻, P.network.params; init=zero(T))))
Base.isodd(P::BitPermutation{T,GRPNetwork{T}}) where T = isodd(count_ones(mapreduce(first, ⊻, P.network.params; init=zero(T))))
Base.isodd(P::AdjointBitPermutation) = isodd(P.parent)
Base.isodd(P::AbstractBitPermutation) = mapreduce(isodd ∘ length, ⊻, cycles(P); init=false)

"""
    iseven(P::AbstractBitPermutation)

Returns a `Bool` corresponding to whether or not the parity of the permutation is even.

See also [`isodd`](@ref).
"""
Base.iseven(P::AbstractBitPermutation) = !isodd(P)
