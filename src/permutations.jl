"""
    AbstractBitPermutation{T,A}

Abstract bit permutation supertype.
"""
abstract type AbstractBitPermutation{T<:Base.BitInteger} end

(perm::AbstractBitPermutation{T})(x::T) where T = bitpermute(x, perm)

"""
    BitPermutation{T,A}

Represents a bit permutation of type `T`, performed with a `BitPermutationAlgorithm` of type `A`.
"""
struct BitPermutation{T,A<:BitPermutationAlgorithm{T}} <: AbstractBitPermutation{T}
    alg::A

    function BitPermutation{T}(alg::A) where {T,A<:BitPermutationAlgorithm{T}}
        new{T,A}(alg)
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
    alg = (type){T}(p; kwargs...)
    return BitPermutation{T}(alg) 
end

_default_type(::Type{T}) where T<:Union{UInt32,UInt64} = ifelse(USE_BMI2, GRPNetwork, BenesNetwork)
_default_type(_::Type{T}) where T = BenesNetwork

Base.show(io::IO, ::BitPermutation{T}) where T = print(io, "BitPermutation{$T}")

function Base.show(io::IO, ::MIME"text/plain", perm::BitPermutation{T}) where T
    println(io, "BitPermutation{$T} with algorithm:")
    print(io, " ", perm.alg)
    return nothing
end

"""
    bitpermute(x::T, p::AbstractBitPermutation{T})
    bitpermute(x::T, n::BitPermutationAlgorithm{T})

Permute bits in `x` using a an `AbstractBitPermutation` or a `BitPermutationAlgorithm`.
Callable as well as `p(x)`.

See also [`invbitpermute`](@ref).
"""
@inline bitpermute(x::T, P::BitPermutation{T}) where T = bitpermute(x, P.alg)

"""
    invbitpermute(x::T, p::AbstractBitPermutation{T})
    invbitpermute(x::T, n::BitPermutationAlgorithm{T})

Permute bits in `x` using a an `AbstractBitPermutation` or a `BitPermutationAlgorithm`.
Callable as well as `p'(x)`.

See also [`bitpermute`](@ref).
"""
@inline invbitpermute(x::T, P::BitPermutation{T}) where T = invbitpermute(x, P.alg)

"""
    AdjointBitPermutation{T,P}

Represents the adjoint permutation of type `P`, where the regular and inverse permutation are swapped.
"""
struct AdjointBitPermutation{T,P<:AbstractBitPermutation{T}} <: AbstractBitPermutation{T}
    parent::P
end

Base.adjoint(perm::AbstractBitPermutation) = AdjointBitPermutation(perm)

@inline bitpermute(x::T, P::AdjointBitPermutation{T}) where T = invbitpermute(x, P.parent)
@inline invbitpermute(x::T, P::AdjointBitPermutation{T}) where T = bitpermute(x, P.parent)

"""
    CompiledBitPermutation{T,F,F̄}

Represents a bit permutation of type `T`, where the regular and inverse permutations are functions of type `F` and `F̄` respectively.
"""
struct CompiledBitPermutation{T,F<:Function,F̄<:Function} <: AbstractBitPermutation{T}
    regular::F
    inverse::F̄

    function CompiledBitPermutation{T}(regular::Function, inverse::Function) where T
        return new{T,typeof(regular),typeof(inverse)}(regular, inverse)
    end
end

@inline bitpermute(x::T, P::CompiledBitPermutation{T}) where T = P.regular(x)
@inline invbitpermute(x::T, P::CompiledBitPermutation{T}) where T = P.inverse(x)

