"""
    CompiledBitPermutation{T,F,F̄}

Represents a bit permutation of type `T`, where the regular and inverse permutations are functions of type `F` and `F̄` respectively.
"""
struct CompiledBitPermutation{T,F<:Function,F̄<:Function} <: AbstractBitPermutation{T}
    vector::Vector{Int}
    regular::F
    inverse::F̄

    function CompiledBitPermutation{T}(vector, regular::Function, inverse::Function) where T
        return new{T,typeof(regular),typeof(inverse)}(vector, regular, inverse)
    end
end

@inline bitpermute(x::T, P::CompiledBitPermutation{T}) where T = P.regular(x)
@inline invbitpermute(x::T, P::CompiledBitPermutation{T}) where T = P.inverse(x)

@inline Vector(P::CompiledBitPermutation) = P.vector

"""
    @bitpermutation T list
    @bitpermutation(T, list)

Generates as `CompiledBitPermutation{T}` from the permutation defined in `list`, which can be a `Vector` or a `Tuple`.
The permutation is first computed as a `BitPermutation`, and the corresponding instructions for `bitpermute` and `invbitpermute` are converted to code, which the compiler will hopefully optimize.

See also [`CompiledBitPermutation`](@ref), `BitPermutation`](@ref), [`bitpermute`](@ref), [`invbitpermute`](@ref).
"""
macro bitpermutation(typesymb, permexpr)
    # Parse input arguments
    typeof(typesymb) ∈ (Symbol, DataType) || 
        throw(ArgumentError("Expected Symbol in first argument, got $(typeof(typesymb))"))
    T = eval(typesymb)
    typeof(permexpr) === Expr && permexpr.head ∈ (:vect, :tuple) ||
        throw(ArgumentError("Expected Vector or Tuple as second argument"))
    p = convert(Vector{Int}, permexpr.args)

    # Construct network
    perm = BitPermutation{T}(p)

    # Concatenate layer operations with `|>`
    regular, inverse = _ops(perm.network)

    # Return something like:
    # "x -> x |> (x′ -> deltaswap(x′,...)) |> (x′ -> deltaswap(x′,...)) |> ..."
    return esc(:(CompiledBitPermutation{$T}($p, $regular, $inverse)))
end

function _ops(net::BenesNetwork)
    freg = mapfoldl((f₁, f₂) -> :($f₁ |> $f₂), net.params; init=:(x)) do args
        return :(x′ -> deltaswap(x′, $(args...)))
    end
    finv = mapfoldl((f₁, f₂) -> :($f₁ |> $f₂), Iterators.reverse(net.params); init=:(x)) do args
        return :(x′ -> deltaswap(x′, $(args...)))
    end

    return :(x -> $freg), :(x -> $finv)
end

function _ops(net::GRPNetwork)
    freg = mapfoldl((f₁, f₂) -> :($f₁ |> $f₂), net.params; init=:(x)) do args
        return :(x′ -> grpswap(x′, $(args...)))
    end

    finv = mapfoldl((f₁, f₂) -> :($f₁ |> $f₂), Iterators.reverse(net.params); init=:(x)) do args
        return :(x′ -> invgrpswap(x′, $(args...)))
    end

    return :(x -> $freg), :(x -> $finv)
end
