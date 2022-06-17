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
    regular, inverse = _ops(perm.alg)

    # Return something like:
    # "x -> x |> (x′ -> deltaswap(x′,...)) |> (x′ -> deltaswap(x′,...)) |> ..."
    return esc(:(CompiledBitPermutation{$T}($regular, $inverse)))
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
