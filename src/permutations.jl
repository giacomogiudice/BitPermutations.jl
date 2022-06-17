
abstract type AbstractPermutation{T<:Base.BitInteger} end

(perm::AbstractPermutation{T})(x::T) where T = bitpermute(x, perm)

struct BitPermutation{T,A<:PermutationAlgorithm{T}} <: AbstractPermutation{T}
    alg::A

    function BitPermutation{T}(alg::A) where {T,A<:PermutationAlgorithm{T}}
        new{T,A}(alg)
    end
end

_default_type(::Type{T}) where T<:Union{UInt32,UInt64} = ifelse(USE_BMI2, GRPNetwork, BenesNetwork)
_default_type(_::Type{T}) where T = BenesNetwork

function BitPermutation{T}(p::AbstractVector{Int}; type::Type=_default_type(T)) where T
    alg = (type){T}(p)
    return BitPermutation{T}(alg) 
end

function Base.show(io::IO, perm::BitPermutation)
    print(typeof(perm))
    return nothing
end

@inline bitpermute(x::T, P::BitPermutation{T}) where T = bitpermute(x, P.alg)
@inline invbitpermute(x::T, P::BitPermutation{T}) where T = invbitpermute(x, P.alg)

struct AdjointBitPermutation{T,P<:AbstractPermutation{T}} <: AbstractPermutation{T}
    parent::P
end

Base.adjoint(perm::AbstractPermutation) = AdjointBitPermutation(perm)

@inline bitpermute(x::T, P::AdjointBitPermutation{T}) where T = invbitpermute(x, P.parent)
@inline invbitpermute(x::T, P::AdjointBitPermutation{T}) where T = invbitpermute(x, P.parent)

struct CompiledBitPermutation{T,F<:Function,F̄<:Function} <: AbstractPermutation{T}
    regular::F
    inverse::F̄

    function CompiledBitPermutation{T}(regular::Function, inverse::Function) where T
        return new{T,typeof(regular),typeof(inverse)}(regular, inverse)
    end
end

@inline bitpermute(x::T, P::CompiledBitPermutation{T}) where T = P.regular(x)
@inline invbitpermute(x::T, P::CompiledBitPermutation{T}) where T = P.inverse(x)

