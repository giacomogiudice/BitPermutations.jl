using BitPermutations
using Test
using Random

# Set seed for reproducibility
Random.seed!(42)

# For permutations it is more natural to think of bits from left to right
bitstr(x) = reverse(bitstring(x)) 

# Check if using BMI2
@info "USE_BMI2 = $(BitPermutations.USE_BMI2), ENV[\"BP_USE_BMI2\"] = $(get(ENV, "BP_USE_BMI2", nothing))"

# Test types
testtypes = (UInt8, UInt16, UInt32, UInt64, UInt128)

@testset "MBitVector{$T}" for T in testtypes
    x = rand(T)
    @test (@inferred bitsize(x)) === 8*sizeof(x)

    # Construction
    v = @inferred MBitVector(x)
    @test @inferred MBitVector{T}(x for x in v) == v
    @test @inferred all(((i, x),) -> v[i] === x, enumerate(v))
    bitstr(x) == string(v)

    # Comparison
    w = MBitVector(v)
    @test (@inferred w == v) && (@inferred isequal(w, v))
    v₀ = MBitVector(zero(T))
    @test (v ≥ v₀) && !any(v₀) && !all(v₀)
    @test any(v) ? (v₀ < v) : (v === v₀)
    @test (@inferred sum(v)) === count_ones(x) && (@inferred sum(.~v)) === count_zeros(x)

    # Logic
    y = rand(T)
    w = MBitVector(y)
    @test convert(T, ~v) === ~x
    @test convert(T, v & w) === x & y
    @test convert(T, v | w) === x | y
    @test convert(T, v ⊻ w) === x ⊻ y
    i = rand(axes(x, 1))
    @test convert(T, circshift(v, i)) === bitrotate(x, i)

    # Assignment
    for _ in 1:20
        i = rand(axes(x, 1))
        b = rand(Bool)
        v[i] = b
        @test v[i] === b
    end
end

@testset "Primitives for type $T" for T in testtypes
    # Simple deltaswaps
    a, b, c, d, e = T(1), T(2), T(3), T(5), T(12)
    @test [bitstr(x)[1:4] for x ∈ (a, b, c, d, e)] == ["1000","0100","1100","1010","0011"]
    @test (@inferred deltaswap(a, a, 1)) === b     # 1̄000... -> 0100...
    @test (@inferred deltaswap(d, b, 1)) === c     # 10̄10... -> 1100...
    @test (@inferred deltaswap(c, c, 2)) === e     # 1̄1̄00... -> 0011...
    
    # Compare with array version
    for _ in 1:20
        x = rand(T)
        s = rand(1:4)
        m = rand(T); m = (m & ~(m << s)) >> s
        v = @inferred deltaswap(MBitVector(x), MBitVector(m), s)
        v′ = @inferred MBitVector(deltaswap(x, m, s))
        @test v == v′
        @test (@inferred deltaswap(deltaswap(x, m, s), m, s)) === x
    end

    # GRP swaps
    a, b, c, d, e = T(1), T(2), T(3), T(6), T(14)
    @test [bitstr(x)[1:4] for x ∈ (a, b, c, d, e)] == ["1000","0100","1100","0110","0111"]
    @test (@inferred grpswap(a, ~c)) === a     # 1̄000... -> 1000...
    @test (@inferred grpswap(b, ~b)) === a     # 1̄000... -> 1000...
    @test (@inferred grpswap(d, ~e)) === c     # 1̄11̄1... -> 1100...

    # Compare with array version
    for _ in 1:20
        x = rand(T)
        m = rand(T)
        v = @inferred grpswap(MBitVector(x), MBitVector(m))
        v′ = @inferred MBitVector(grpswap(x, m))
        @test v == v′
        @test (@inferred invgrpswap(grpswap(x, m), m)) === x
        v′ = @inferred invgrpswap(v, MBitVector(m))
        @test v′ == MBitVector(x)
    end
end

@testset "BenesNetwork{$T}" for T in testtypes
    for _ in 1:10
        n = rand(1:bitsize(T))
        p = randperm(n)
        p′ = collect(1:bitsize(T))
        p′[1:n] = p′[p]
        net = BenesNetwork{T}(p)
        for _ in 1:100
            x = rand(T)
            @test bitstr(x)[p′] === bitstr(@inferred bitpermute(x, net))
            @test bitstr(x)[invperm(p′)] === bitstr(@inferred invbitpermute(x, net))
        end
    end
end

@testset "GRPNetwork{$T}" for T in testtypes
    for _ in 1:10
        n = rand(1:bitsize(T))
        p = randperm(n)
        p′ = collect(1:bitsize(T))
        p′[1:n] = p′[p]
        net = GRPNetwork{T}(p)
        for _ in 1:100
            x = rand(T)
            @test bitstr(x)[p′] === bitstr(@inferred bitpermute(x, net))
            @test bitstr(x)[invperm(p′)] === bitstr(@inferred invbitpermute(x, net))
        end
    end
end

@testset "BitPermutation{$T}" for T in testtypes
    for _ in 1:10
        n = rand(1:bitsize(T))
        p = randperm(n)
        P₁ = BitPermutation{T}(p; type=BenesNetwork)
        P₂ = BitPermutation{T}(p; type=GRPNetwork)

        # Make sure printing returns something
        buf = IOBuffer()
        show(buf, P₁)
        @test !isempty(String(take!(buf)))
        show(buf, P₂)
        @test !isempty(String(take!(buf)))
        show(buf, "text/plain", P₁)
        @test !isempty(String(take!(buf)))
        show(buf, "text/plain", P₂)
        @test !isempty(String(take!(buf)))

        for _ in 1:20
            x = rand(T)
            @test P₁(x) === P₂(x) === bitpermute(x, P₁) === bitpermute(x, P₂)
            @test P₁'(x) === P₂'(x) === invbitpermute(x, P₁) === invbitpermute(x, P₂)
            @test P₁(x) === invbitpermute(x, P₁') === P₂(x) === invbitpermute(x, P₂')
        end
    end
end

@testset "CompiledBitPermutation{$T}" for T in testtypes
    P₁ = @eval @bitpermutation $T (2, 6, 5, 8, 4, 7, 1, 3)
    P₂ = BitPermutation{T}([2, 6, 5, 8, 4, 7, 1, 3])
    for _ in 1:100
        x = rand(T)
        @test P₁(x) === P₂(x) === bitpermute(x, P₁) === bitpermute(x, P₂)
        @test P₁'(x) === P₂'(x) === invbitpermute(x, P₁) === invbitpermute(x, P₂)
    end
end
