using BitPermutations
using BenchmarkTools
using Random

# Set seed for reproducibility
Random.seed!(42)

# Check if using BMI2
@info "USE_BMI2 = $(BitPermutations.USE_BMI2)"

# Struct in which to save benchmark results
suite = BenchmarkGroup()
test_types = (UInt8, UInt16, UInt32, UInt64, UInt128)
N = 10_000

# Generate results for benchmarking bit permutations on a single bitstring of type `T`
suite["Single"] = BenchmarkGroup()

rand_perm(::Type{T}) where {T} = shuffle!(collect(1:bitsize(T)))
benes_perm(::Type{T}) where {T} = BitPermutation{T}(rand_perm(T); type=BenesNetwork)
grp_perm(::Type{T}) where {T} = BitPermutation{T}(rand_perm(T); type=GRPNetwork)

for T in test_types
    group = suite["Single"][T] = BenchmarkGroup()

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
    group["Bits"] = @benchmarkable begin
        @inbounds for _ in 1:($N)
            permute!(b, p)
        end
        return b
    end setup = begin
        b = Bits(rand($T))
        p = rand_perm($T)
    end
    group["BitVector"] = @benchmarkable begin
        @inbounds for _ in 1:($N)
            permute!(v, p)
        end
        return v
    end setup = begin
        v = BitVector(Bits(rand($T)))
        p = rand_perm($T)
    end
end

suite["Broadcasted"] = BenchmarkGroup()

for T in test_types
    group = suite["Broadcasted"][T] = BenchmarkGroup()

    group["BenesNetwork"] = @benchmarkable map!(p, xs) setup = begin
        xs = rand($T, $N)
        p = benes_perm($T)
    end
    group["GRPNetwork"] = @benchmarkable map!(p, xs) setup = begin
        xs = rand($T, $N)
        p = grp_perm($T)
    end
    group["Bits"] = @benchmarkable begin
        foreach(b -> permute!(b, p), bs)
        return bs
    end setup = begin
        bs = [Bits(rand($T)) for _ in 1:($N)]
        p = rand_perm($T)
    end
    group["BitVector"] = @benchmarkable begin
        foreach(v -> permute!(v, p), vs)
        return vs
    end setup = begin
        vs = [BitVector(Bits(rand($T))) for _ in 1:($N)]
        p = rand_perm($T)
    end
end

results = run(suite; verbose=true)

if !isempty(ARGS)
    fname = first(ARGS)
    @info "Saving median of results as $fname"
    BenchmarkTools.save(fname, median(results))
end
