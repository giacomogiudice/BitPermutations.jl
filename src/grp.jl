struct GRPNetwork{T}
    ops::Vector{Tuple{T,T}}
end

function GRPNetwork{T}(perm::AbstractVector{Int}) where T
    n = length(perm)
    isbitstype(T) && n ≤ bitsize(T) || throw(OverflowError("Permutation is too long for Type $T"))
    isperm(perm) || throw(ArgumentError("Input vector is not a permutation"))

    USE_BMI2 || @warn "Not using BMI2 instructions, performance may be limited" maxlog=1

    # Pad inverted permutation vector
    p = vcat(invperm(perm), n+1:bitsize(T))

    masks = zeros(T, trailing_zeros(bitsize(T)))

    for ind in 1:length(masks)
        k = 2^(ind-1)
        mask = MBitVector(zero(T))
        for (i, pᵢ) in enumerate(p)
            iseven(cld(pᵢ, k)) && (mask[i] = 1)
        end
        p = arraygrpswap(p, mask)
        masks[ind] = chunk(mask)
    end
    @assert issorted(p)
    return GRPNetwork{T}([(m, ~m) for m in masks if !iszero(m)])
end

function Base.show(io::IO, ::MIME"text/plain", net::GRPNetwork)
    println("$(typeof(net)) with $(length(net.ops)) GRP-swaps")
    return nothing
end

function bitpermute(net::GRPNetwork{T}, x::T) where T
    # Shift should always be bitsize(T) ÷ 2
    s = bitsize(T) >> 1
    return foldl(net.ops; init=x) do x′, (m, m̄)
        return grpswap(x′, m, s, m̄)
    end
end

function invbitpermute(net::GRPNetwork{T}, x::T) where T
    # Shift should always be bitsize(T) ÷ 2
    s = bitsize(T) >> 1
    return foldl(Iterators.reverse(net.ops); init=x) do x′, (m, m̄)
        return invgrpswap(x′, m, s, m̄)
    end
end

function grpswap(x::T, m::T, shift::Int=count_zeros(m), m̄::T=~m) where T
    return pext(x, m) << shift | pext(x, m̄)
end

function invgrpswap(x::T, m::T, shift::Int=count_zeros(m), m̄::T=~m) where T
    return pdep(x >> shift, m) | pdep(x, m̄)
end

function arraygrpswap(x::AbstractVector, m::AbstractVector{Bool})
    return vcat(x[.~m], x[m])
end

function arrayinvgrpswap(x::AbstractVector, m::AbstractVector{Bool})
    s = length(m) - sum(m)
    y = similar(x)
    y[.~m] = x[begin:(begin+s-1)]
    y[m] = x[(begin+s):end]
    return y
end
