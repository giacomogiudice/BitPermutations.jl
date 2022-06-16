struct GRPNetwork{T} <: PermutationNetwork{T}
    params::Vector{Tuple{T,Int,T}}
end

# Algorithm explained in https://programming.sirrida.de/bit_perm.html#lee_sag
function GRPNetwork{T}(perm::AbstractVector{Int}) where T
    n = length(perm)
    isbitstype(T) && n ≤ bitsize(T) || throw(OverflowError("Permutation is too long for Type $T"))
    isperm(perm) || throw(ArgumentError("Input vector is not a permutation"))
    USE_BMI2 || @warn "Not using BMI2 instructions, performance may be limited" maxlog=1

    # Pad permutation vector
    p = [collect(perm); n+1:bitsize(T)]

    # Build partial lists
    lists = [Int[]]
    k = 1
    for (i, pᵢ) in enumerate(p)
        # Go to next list if previous item is larger
        if !isempty(lists[k]) && lists[k][end] > pᵢ
            push!(lists, Int[])
            k += 1
        end
        push!(lists[k], pᵢ)
    end

    # Merge lists pairwise until there is only one left
    masks = MBitVector[]
    while length(lists) > 1
        # Make sure `lists` has an even size
        isodd(length(lists)) && push!(lists, Int[])
        nhalf = length(lists) ÷ 2
        mergedlists = Vector{Vector{Int}}(undef, nhalf)
        partialmasks = Vector{BitVector}(undef, nhalf)
        
        # Sequentially merge the lists
        for k in 1:nhalf
            l, r = lists[k], lists[k+nhalf]
            mergedlists[k] = [l; r]
            partialmasks[k] = [falses(length(l)); trues(length(r))]
            inds = sortperm(mergedlists[k])
            mergedlists[k] = mergedlists[k][inds]
            partialmasks[k] = partialmasks[k][inds]
        end

        # Push mask to front since order of masks should be reversed
        lists = mergedlists
        pushfirst!(masks, MBitVector{T}(vcat(partialmasks...)))
    end
    @assert issorted(first(lists))

    return GRPNetwork{T}([(m, count_zeros(m), ~m) for m in Iterators.map(chunk, masks)])
end

function Base.show(io::IO, ::MIME"text/plain", net::GRPNetwork)
    println("$(typeof(net)) with $(length(net.params)) operations")
    return nothing
end

function bitpermute(x::T, net::GRPNetwork{T}) where T
    return foldl(net.params; init=x) do x′, args
        return grpswap(x′, args...)
    end
end

function invbitpermute(x::T, net::GRPNetwork{T}) where T
    return foldl(Iterators.reverse(net.params); init=x) do x′, args
        return invgrpswap(x′, args...)
    end
end

@inline function grpswap(x::T, m::T, shift::Int=count_zeros(m), m̄::T=~m) where T<:Integer
    return pext(x, m) << shift | pext(x, m̄)
end

@inline function invgrpswap(x::T, m::T, shift::Int=count_zeros(m), m̄::T=~m) where T<:Integer
    return pdep(x >> shift, m) | pdep(x, m̄)
end

function grpswap(x::AbstractVector, m::AbstractVector{Bool})
    return vcat(x[.~m], x[m])
end

function invgrpswap(x::AbstractVector, m::AbstractVector{Bool})
    s = sum(.~m)
    y = similar(x)
    y[.~m] = x[begin:(begin+s-1)]
    y[m] = x[(begin+s):end]
    return y
end
