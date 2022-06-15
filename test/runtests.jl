using BitPermutations
using Test
using Random

# Set seed for reproducibility
Random.seed!(42)

# For permutations it is more natural to think of bits from left to right
bitvec(x) = reverse(bitstring(x)) 

# Test types
testtypes = (UInt8, UInt16, UInt32, UInt64, UInt128)

@testset "MBitVector{$T}" for T in testtypes
    x = rand(T)
    @test (@inferred bitsize(x)) === 8*sizeof(x)

    v = @inferred MBitVector(x)
    @test @inferred MBitVector{T}(x for x in v) == v
    @test @inferred all(((i, x),) -> v[i] === x, enumerate(v))
    bitvec(x) == string(v)

    y = rand(T)
    w = MBitVector(y)

    @test convert(T, ~v) === ~x
    @test convert(T, v & w) === x & y
    @test convert(T, v | w) === x | y
    @test convert(T, v ⊻ w) === x ⊻ y
    i = rand(axes(x, 1))
    @test convert(T, circshift(v, i)) === bitrotate(x, i)

    # Assignment
    for _ in 1:10
        i = rand(axes(x, 1))
        b = rand(Bool)
        v[i] = b
        @test v[i] === b
    end
end

@testset "Primitives for type $T" for T in testtypes
    # Simple deltaswaps
    a, b, c, d, e = T(1), T(2), T(3), T(5), T(12)
    @test [bitvec(x)[1:4] for x ∈ (a, b, c, d, e)] == ["1000","0100","1100","1010","0011"]
    @test (@inferred deltaswap(a, a, 1)) === b     # 1̄000... -> 0100...
    @test (@inferred deltaswap(d, b, 1)) === c     # 10̄10... -> 1100...
    @test (@inferred deltaswap(c, c, 2)) === e     # 1̄1̄00... -> 0011...
    
    # Compare with array version
    for _ in 1:20
        x = rand(T)
        s = rand(1:4)
        m = rand(T); m = (m & ~(m << s)) >> s
        v = @inferred BitPermutations.arraydeltaswap(MBitVector(x), MBitVector(m), s)
        v′ = @inferred MBitVector(deltaswap(x, m, s))
        @test v == v′
        @test (@inferred deltaswap(deltaswap(x, m, s), m, s)) === x
    end

    # GRP swaps
    a, b, c, d, e = T(1), T(2), T(3), T(6), T(14)
    @test [bitvec(x)[1:4] for x ∈ (a, b, c, d, e)] == ["1000","0100","1100","0110","0111"]
    @test (@inferred grpswap(a, ~c)) === a     # 1̄000... -> 1000...
    @test (@inferred grpswap(b, ~b)) === a     # 1̄000... -> 1000...
    @test (@inferred grpswap(d, ~e)) === c     # 1̄11̄1... -> 1100...

    # Compare with array version
    for _ in 1:20
        x = rand(T)
        m = rand(T)
        v = @inferred BitPermutations.arraygrpswap(MBitVector(x), MBitVector(m))
        v′ = @inferred MBitVector(grpswap(x, m))
        @test v == v′
        @test (@inferred invgrpswap(grpswap(x, m), m)) === x
        v′ = @inferred BitPermutations.arrayinvgrpswap(v, MBitVector(m))
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
            @test bitvec(x)[p′] === bitvec(@inferred bitpermute(net, x))
            @test bitvec(x)[invperm(p′)] === bitvec(@inferred invbitpermute(net, x))
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
            @test bitvec(x)[p′] === bitvec(@inferred bitpermute(net, x))
            @test bitvec(x)[invperm(p′)] === bitvec(@inferred invbitpermute(net, x))
        end
    end
end
