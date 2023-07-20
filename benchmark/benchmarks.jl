using BitPermutations
using BenchmarkTools
using Random

# Set seed for reproducibility
Random.seed!(42)

# Check if using instrinsics
@assert BitPermutations.USE_BMI2 && BitPermutations.USE_AVX512

# Struct in which to save benchmark results
suite = BenchmarkGroup()
test_types = (UInt8, UInt16, UInt32, UInt64, UInt128)
avx_types = (UInt16, UInt32, UInt64)
N = 10_000

# Generate results for benchmarking bit permutations on a single bitstring of type `T`
suite["Single"] = BenchmarkGroup()

rand_perm(::Type{T}) where {T} = randperm(bitsize(T))
benes_perm(::Type{T}) where {T} = BitPermutation{T}(randperm(bitsize(T)); type=BenesNetwork)
grp_perm(::Type{T}) where {T} = BitPermutation{T}(randperm(bitsize(T)); type=GRPNetwork)
avx_perm(::Type{T}) where {T} = BitPermutation{T}(randperm(bitsize(T)); type=AVXCopyGather)

for T in test_types
    group = suite["Single"][T] = BenchmarkGroup()

    group["BitVector"] = @benchmarkable begin
        @inbounds for _ in 1:($N)
            permute!(v, p)
        end
        return v
    end setup = begin
        v = BitVector(Bits(rand($T)))
        p = rand_perm($T)
    end
    group["Bits"] = @benchmarkable begin
        @inbounds for _ in 1:($N)
            permute!(b, p)
        end
        return b
    end setup = begin
        b = Bits(rand($T))
        p = rand_perm($T)
    end
    group["BenesNetwork"] = @benchmarkable begin
        @inbounds for _ in 1:($N)
            x = bitpermute(x, p)
        end
        return x
    end setup = begin
        x = rand($T)
        p = benes_perm($T)
    end
    group["GRPNetwork"] = @benchmarkable begin
        @inbounds for _ in 1:($N)
            x = bitpermute(x, p)
        end
        return x
    end setup = begin
        x = rand($T)
        p = grp_perm($T)
    end
    T in avx_types || continue
    group["AVXCopyGather"] = @benchmarkable begin
        @inbounds for _ in 1:($N)
            x = bitpermute(x, p)
        end
        return x
    end setup = begin
        x = rand($T)
        p = avx_perm($T)
    end
end

suite["Broadcasted"] = BenchmarkGroup()

for T in test_types
    group = suite["Broadcasted"][T] = BenchmarkGroup()

    group["BitVector"] = @benchmarkable begin
        foreach(v -> permute!(v, p), vs)
        return vs
    end setup = begin
        vs = [BitVector(Bits(rand($T))) for _ in 1:($N)]
        p = rand_perm($T)
    end
    group["Bits"] = @benchmarkable begin
        foreach(b -> permute!(b, p), bs)
        return bs
    end setup = begin
        bs = [Bits(rand($T)) for _ in 1:($N)]
        p = rand_perm($T)
    end
    group["BenesNetwork"] = @benchmarkable bitpermute!(xs, p) setup = begin
        xs = rand($T, $N)
        p = benes_perm($T)
    end
    group["GRPNetwork"] = @benchmarkable bitpermute!(xs, p) setup = begin
        xs = rand($T, $N)
        p = grp_perm($T)
    end
    T in avx_types || continue
    group["AVXCopyGather"] = @benchmarkable bitpermute!(xs, p) setup = begin
        xs = rand($T, $N)
        p = avx_perm($T)
    end
end

results = run(suite; verbose=true)

if !isempty(ARGS)
    fname = first(ARGS)
    @info "Saving median of results as $fname"
    BenchmarkTools.save(fname, median(results))
end
