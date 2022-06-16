struct BenesNetwork{T} <: PermutationNetwork{T}
    params::Vector{Tuple{T,Int}}
end

function BenesNetwork{T}(perm::AbstractVector{Int}; verbose::Bool=false, rearrange::Bool=true) where T
    n = length(perm)
    isbitstype(T) && n ≤ bitsize(T) || throw(OverflowError("Permutation is too long for Type $T"))
    isperm(perm) || throw(ArgumentError("Input vector is not a permutation"))

    # Pad permutation vector
    p = [perm; n+1:bitsize(T)]

    shiftset = reverse([2^i for i in 0:trailing_zeros(bitsize(T)÷2)])

    !rearrange && return BenesNetwork(_ops(T, p, shiftset))
    params = nothing
    for shifts in permutations(shiftset)
        newparams = _params(T, p, shifts)
        verbose && @info "Permutation $(shifts) requires $(length(newparams)) operations"
        if isnothing(params) || length(newparams) < length(params)
            params = newparams
        end
    end
    return BenesNetwork(params)
end

function Base.show(io::IO, ::MIME"text/plain", net::BenesNetwork)
    println("$(typeof(net)) with $(length(net.params)) operations")
    return nothing
end

function _params(::Type{T}, p::AbstractVector{Int}, shifts) where T
    # Compute masks
    masks = _masks(T, p, shifts)
    # Merge innermost masks
    m = length(shifts)
    masks[m] = masks[m] ⊻ masks[m+1]
    popat!(masks, m+1)
    deltas = [shifts; reverse(shifts[begin:end-1])]
    return [(chunk(m), δ) for (m, δ) in zip(masks, deltas) if any(m)]
end

function _masks(::Type{T}, p::AbstractVector{Int}, shifts) where T
    masks = Vector{MBitVector{T}}(undef, 2*length(shifts))
    source = collect(1:length(p))
    target = collect(p)
    for (ind, shift) in enumerate(shifts)
        mask_forward, mask_backward = _stagemasks(T, target, shift)
        masks[ind], masks[end-(ind-1)] = mask_forward, mask_backward

        p_forward = deltaswap(source, mask_forward, shift)
        p_backward = deltaswap(target, mask_backward, shift)

        # New permutation for next stage
        target = p_forward[p_backward]
    end
    return masks
end

function _stagemasks(::Type{T}, p::AbstractVector, shift::Int) where T
    # All nodes must be repartitioned in two sets with the same number of elements
    # encoded with {true|false}, corresponding to which partition they are destined to
    n = length(p)
    ax = first(axes(p))
    p̄ = invperm(p)

    # Keeps track of visited nodes
    visited = MBitVector(zero(T))
    # Partition of next stage
    partition = MBitVector{T}(c for _ in 1:(n÷2shift) for c in 0:1 for _ in 1:shift)
    # Instantiate masks
    mask_forward = MBitVector(zero(T))
    mask_backward = MBitVector(zero(T))

    # Heuristic: start from the smallest index which is not moved
    i = 1
    for (j, pⱼ) in enumerate(p)
        if pⱼ === j
            i = j
            break 
        end
    end

    # Encodes the beginning of each cycle
    i₀ = i
    
    while !isnothing(i)
        # Heuristic to have fewer masks: try to not reshuffle indices
        color = partition[i₀]

        # Compute mask of i on forward side:
        # move i if color is different than partition
        mask_forward[i] = partition[i] ≠ color
        visited[i] = true
        # Compute mask of i on backward side:
        # i comes from partition `color`, so move i if on wrong partition
        ī = p̄[i]
        mask_backward[ī] = partition[ī] ≠ color

        # Conjugate on forward side has opposite color
        j = _conjugate(i, partition, shift)
        visited[j] = true

        # Compute mask of j on backward side:
        # j comes from partition `!color`, so move j if on wrong partition
        j̄ = p̄[j]
        mask_backward[j̄] = partition[j̄] ≠ !color

        # Move to conjugate on backward side
        j′ = p[_conjugate(p̄[j], partition, shift)]
        if j′ == i₀
            # Cycle is over
            i = findfirst(iszero, visited)
            i₀ = i
        else
            # Keep following cycle
            i = j′
        end
    end

    # Normalize mask
    mask_forward = (circshift(mask_forward, -shift) | mask_forward) & ~partition
    mask_backward = (circshift(mask_backward, -shift) | mask_backward) & ~partition
    return mask_forward, mask_backward
end

function _conjugate(ind::Int, partition::AbstractVector, shift::Int)
    ax = first(axes(partition))
    @boundscheck checkbounds(ax, ind)
    ret = partition[ind] ? ind - shift : ind + shift
    @boundscheck checkbounds(ax, ret)
    return ret
end


function bitpermute(x::T, net::BenesNetwork{T}) where T
    return foldl(net.params; init=x) do x′, (mask, shift)
        return deltaswap(x′, mask, shift)
    end
end

function invbitpermute(x::T, net::BenesNetwork{T}) where T
    return foldl(Iterators.reverse(net.params); init=x) do x′, (mask, shift)
        return deltaswap(x′, mask, shift)
    end
end

# Swaps bits in x selected by mask m with ones to the left by an amount `shift`
@inline function deltaswap(x::T, m::T, shift::Int) where T<:Integer
    t = ((x >> shift) ⊻ x) & m
    return x ⊻ t ⊻ (t << shift)
end

# Slow for arrays
function deltaswap(x::AbstractVector, m::AbstractVector{Bool}, shift::Int)
    @assert length(x) === length(m)
    y = copy(x)
    for (i, mᵢ) in enumerate(m)
        if !iszero(mᵢ)
            y[i+shift], y[i] = x[i], x[i+shift]
        end
    end
    return y
end
