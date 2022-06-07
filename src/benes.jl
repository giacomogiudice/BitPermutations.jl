struct BenesNetwork{T}
    ops::Vector{Tuple{MBitVector{T},Int}}
end

function BenesNetwork{T}(perm::AbstractVector{Int}; verbose::Bool=false, rearrange::Bool=true) where T
    n = length(perm)
    isbitstype(T) && n ≤ bitsize(T) || throw(OverflowError("Permutation is too long for Type $T"))
    # ispow2(n) || throw(ArgumentError("Permutation must be exactly a power of 2"))

    # Pad permutation vector
    p = vcat(perm, n+1:bitsize(T))

    shiftset = [2^i for i in 0:Int(log2(bitsize(T)÷2))]
    !rearrange && return BenesNetwork(_ops(T, p, shifts))

    ops = nothing
    for shifts in permutations(shiftset)
        newops = _ops(T, p, shifts)
        verbose && @info "Permutation $(shifts) requires $(length(newops)) operations"
        if isnothing(ops) || length(newops) < length(ops)
            ops = newops
        end
    end
    return BenesNetwork(ops)
end

function Base.show(io::IO, ::MIME"text/plain", net::BenesNetwork)
    println("$(typeof(net)) with $(length(net.ops)) deltaswaps")
    return nothing
end

function _ops(::Type{T}, p::AbstractVector{Int}, shifts) where T
    # Compute masks
    masks = _masks(T, p, shifts)
    # Merge innermost masks
    m = length(shifts)
    masks[m] = masks[m] ⊻ masks[m+1]
    popat!(masks, m+1)
    deltas = [shifts; reverse(shifts[begin:end-1])]
    return [(m, δ) for (m, δ) in zip(masks, deltas) if any(m)]
end

function _masks(::Type{T}, p::AbstractVector{Int}, shifts) where T
    masks = Vector{MBitVector{T}}(undef, 2*length(shifts))
    source = collect(1:length(p))
    target = collect(p)
    for (ind, shift) in enumerate(shifts)
        mask_forward, mask_backward = _stagemasks(T, target, shift)
        masks[ind], masks[end-(ind-1)] = mask_forward, mask_backward

        p_forward = arraydeltaswap(source, mask_forward, shift)
        p_backward = arraydeltaswap(target, mask_backward, shift)

        # New permutation for next stage
        target = p_forward[p_backward]
    end
    return masks
end

function _conjugate(ind::Int, partition::AbstractVector, shift::Int)
    ax = first(axes(partition))
    @boundscheck checkbounds(ax, ind)
    ret = partition[ind] ? ind - shift : ind + shift
    @boundscheck checkbounds(ax, ret)
    return ret
end


function _stagemasks(::Type{T}, p::AbstractVector, shift::Int) where T
    n = length(p)
    ax = first(axes(p))
    p̄ = invperm(p)

    # Keeps track of visited nodes
    visited = MBitVector{T}(zero(T))
    # Partition of next stage
    partition = MBitVector{T}(c for _ in 1:(n÷2shift) for c in 0:1 for _ in 1:shift)
    # Instantiate masks
    mask_forward = MBitVector{T}(zero(T))
    mask_backward = MBitVector{T}(zero(T))


    # All nodes must be repartitioned in two sets with the same number of elements
    # encoded with {true|false}, corresponding to which partition they are destined to
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

function bitpermute(net::BenesNetwork{T}, x::T) where T
    return foldl(net.ops; init=x) do x′, (mask, shift)
        return deltaswap(x′, mask, shift)
    end
end

function invbitpermute(net::BenesNetwork{T}, x::T) where T
    return foldl(Iterators.reverse(net.ops); init=x) do x′, (mask, shift)
        return deltaswap(x′, mask, shift)
    end
end

# Swaps bits in x selected by mask m with ones to the left by an amount `shift`
function deltaswap(x::Unsigned, m::Unsigned, shift::Int)
    t = ((x >> shift) ⊻ x) & m
    return x ⊻ t ⊻ (t << shift)
end

deltaswap(x::T, m::MBitVector{T}, shift::Int) where T = deltaswap(x, chunk(m), shift)

function arraydeltaswap(x::AbstractVector, m::MBitVector, shift::Int)
    @assert length(x) === length(m)
    y = copy(x)
    for (i, mᵢ) in enumerate(m)
        if !iszero(mᵢ)
            y[i+shift], y[i] = x[i], x[i+shift]
        end
    end
    return y
end

# Internal debugging
bitvec(x) = reverse(bitstring(x))
