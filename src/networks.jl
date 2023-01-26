"""
    PermutationNetwork{T}

Abstract type for all bit permutation algorithms.
"""
abstract type PermutationNetwork{T} end

"""
    BenesNetwork{T}

Represents a Beneš network, which is a series of `deltaswap` operations with specified shifts and
masks.
"""
struct BenesNetwork{T} <: PermutationNetwork{T}
    params::Vector{Tuple{T,Int}}
end

function BenesNetwork{T}(perm::AbstractVector{Int}; verbose::Bool=false, rearrange::Bool=true) where {T}
    n = length(perm)
    n ≤ bitsize(T) || throw(OverflowError("Permutation is too long for Type $T"))
    isperm(perm) || throw(ArgumentError("Input vector is not a permutation"))

    # Pad permutation vector
    p = [collect(perm); (n + 1):bitsize(T)]

    # Use `trailing_zeros` to compute the log2 quickly
    shifts = reverse!([2^i for i in 0:trailing_zeros(bitsize(T) ÷ 2)])

    # Initialize masks
    masks = Vector{Bits{T}}(undef, 2 * length(shifts) - 1)

    # Fill with garbage so they are all computed the first time
    prevshifts = zeros(Int, length(shifts))
    bestshifts = zeros(Int, length(shifts))
    bestmasks = similar(masks)
    bestscore = -1

    # Loop over all possible permutations
    iterator = rearrange ? permutations(shifts) : (shifts,)
    source = collect(1:length(p))
    for shifts in iterator
        cached = true
        target = p
        # Populate masks with first arrangement of shifts
        for (index, shift) in enumerate(shifts)
            # Avoid recomputing masks if same shift as before
            if !(cached && shift == prevshifts[index])
                _set_stage_masks!(masks, shifts, target, index)
                cached = false
            end
            # Compute permutation for next stage
            if index ≠ lastindex(shifts)
                mask_forward, mask_backward = masks[index], masks[end - index + 1]
                p_forward = deltaswap(source, mask_forward, shift)
                p_backward = deltaswap(target, mask_backward, shift)
                target = p_forward[p_backward]
            end
        end
        # Score based on number of zero masks
        score = sum(!any, masks)
        verbose && @info "Ordering $shifts has score: $score ($(length(masks) - score) ops)"
        if score > bestscore
            bestscore = score
            bestshifts = copy(shifts)
            bestmasks = copy(masks)
        end

        prevshifts = copy(shifts)
    end

    # Extract parameters
    shifts = bestshifts
    masks = bestmasks
    deltas = [shifts; reverse!(shifts[begin:(end - 1)])]

    return BenesNetwork([(chunk(m), δ) for (m, δ) in zip(masks, deltas) if any(m)])
end

function _set_stage_masks!(
    masks::AbstractVector{Bits{T}}, shifts::AbstractVector{Int}, p::AbstractVector{Int}, index::Int
) where {T}
    # All nodes must be repartitioned in two sets with the same number of elements
    # encoded with {true|false}, corresponding to which partition they are destined to
    n = length(p)
    ax = axes(p, 1)
    p̄ = invperm(p)
    shift = shifts[index]

    # Keeps track of visited nodes
    visited = Bits(zero(T))
    # Partition of next stage
    partition = Bits{T}(c for _ in 1:(n ÷ 2shift) for c in 0:1 for _ in 1:shift)
    # Instantiate masks
    mask_forward = Bits(zero(T))
    mask_backward = Bits(zero(T))

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

    @inbounds while !isnothing(i)
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

    # Normalize mask (move masking bits to MSB position)
    mask_forward = _normalize_mask(mask_forward, partition, shift)
    mask_backward = _normalize_mask(mask_backward, partition, shift)

    # Set masks
    if index ≠ lastindex(shifts)
        masks[index] = mask_forward
        masks[end - index + 1] = mask_backward
    else
        # Last index: merge the innermost masks
        masks[index] = mask_forward .⊻ mask_backward
    end

    return nothing
end

# Moves bits in mask to MSB position 
function _normalize_mask(mask::Bits{T}, partition::Bits{T}, shift::Int) where {T}
    m = chunk(mask)
    p = chunk(partition)
    return Bits{T}(((m >> shift) | m) & ~p)
end

# Finds index shifted up or down by `shift`
function _conjugate(ind::Int, partition::AbstractVector, shift::Int)
    ax = first(axes(partition))
    @boundscheck checkbounds(ax, ind)
    @inbounds ret = partition[ind] ? ind - shift : ind + shift
    @boundscheck checkbounds(ax, ret)
    return ret
end

Base.show(io::IO, net::BenesNetwork) = print(io, "$(typeof(net)) with $(length(net.params)) operations")

function bitpermute(x::Number, net::BenesNetwork{T}) where {T}
    return foldl(net.params; init=convert(T, x)) do x′, (mask, shift)
        return deltaswap(x′, mask, shift)
    end
end

function Base.broadcasted(::typeof(bitpermute), x::AbstractArray, net::BenesNetwork)
    return foldl(net.params; init=x) do x′, (mask, shift)
        return deltaswap.(x′, mask, shift)
    end
end

function invbitpermute(x::Number, net::BenesNetwork{T}) where {T}
    return foldl(Iterators.reverse(net.params); init=convert(T, x)) do x′, (mask, shift)
        return deltaswap(x′, mask, shift)
    end
end

function Base.broadcasted(::typeof(invbitpermute), x::AbstractArray, net::BenesNetwork)
    return foldl(Iterators.reverse(net.params); init=x) do x′, (mask, shift)
        return deltaswap.(x′, mask, shift)
    end
end

"""
    GRPNetwork{T}

Represents a GRP network, which is a series of `grpswap` operations with specified shifts and masks.
"""
struct GRPNetwork{T} <: PermutationNetwork{T}
    params::Vector{Tuple{T,Int,T}}
end

# Algorithm explained in https://programming.sirrida.de/bit_perm.html#lee_sag
function GRPNetwork{T}(perm::AbstractVector{Int}) where {T}
    n = length(perm)
    n ≤ bitsize(T) || throw(OverflowError("Permutation is too long for Type $T"))
    isperm(perm) || throw(ArgumentError("Input vector is not a permutation"))
    USE_BMI2 || @warn "Not using BMI2 instructions, performance may be limited" maxlog = 1

    # Pad permutation vector
    p = [collect(perm); (n + 1):bitsize(T)]

    # Build partial lists of indices
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
    masks = Bits[]
    while length(lists) > 1
        # Make sure `lists` has an even size
        isodd(length(lists)) && push!(lists, Int[])
        nhalf = length(lists) ÷ 2
        mergedlists = Vector{Vector{Int}}(undef, nhalf)
        partialmasks = Vector{BitVector}(undef, nhalf)

        # Sequentially merge the lists
        for k in 1:nhalf
            l, r = lists[k], lists[k + nhalf]
            mergedlists[k] = [l; r]
            partialmasks[k] = [falses(length(l)); trues(length(r))]
            inds = sortperm(mergedlists[k])
            mergedlists[k] = mergedlists[k][inds]
            partialmasks[k] = partialmasks[k][inds]
        end

        # Push mask to front since order of masks should be reversed
        lists = mergedlists
        pushfirst!(masks, Bits{T}(vcat(partialmasks...)))
    end
    @assert issorted(first(lists))

    return GRPNetwork{T}([(m, count_zeros(m), ~m) for m in Iterators.map(chunk, masks)])
end

Base.show(io::IO, net::GRPNetwork) = print(io, "$(typeof(net)) with $(length(net.params)) operations")

function bitpermute(x::Number, net::GRPNetwork{T}) where {T}
    return foldl(net.params; init=convert(T, x)) do x′, args
        return grpswap(x′, args...)
    end
end

function Base.broadcasted(::typeof(bitpermute), x::AbstractArray, net::GRPNetwork)
    return foldl(net.params; init=x) do x′, args
        return grpswap.(x′, args...)
    end
end

function invbitpermute(x::Number, net::GRPNetwork{T}) where {T}
    return foldl(Iterators.reverse(net.params); init=convert(T, x)) do x′, args
        return invgrpswap(x′, args...)
    end
end

function Base.broadcasted(::typeof(invbitpermute), x::AbstractArray, net::GRPNetwork)
    return foldl(Iterators.reverse(net.params); init=x) do x′, args
        return invgrpswap.(x′, args...)
    end
end
