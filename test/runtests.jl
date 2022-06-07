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

    v = MBitVector(x)
    @test all(((i, x),) -> v[i] === x, enumerate(v))
    bitvec(x) == string(v)

    @test convert(T, ~v) === ~x

    y = rand(T)
    w = MBitVector(y)

    @test convert(T, x ⊻ y) === x ⊻ y
    
    i = rand(axes(x, 1))
    v[i] = 1
    @test v[i] === true
end
# @testset "Simple deltaswaps" begin
#     # Simple deltaswaps
#     a, b, c, d, e = T(1), T(2), T(3), T(5), T(12)
#     @test [bitvec(x)[1:4] for x ∈ (a, b, c, d, e)] == ["1000","0100","1100","1010","0011"]
#     @test (@inferred deltaswap(a, a, 1)) === b     # 1̄000... -> 0100...
#     @test (@inferred deltaswap(d, b, 1)) === c     # 10̄10... -> 1100...
#     @test (@inferred deltaswap(c, c, 2)) === e     # 1̄1̄00... -> 0011...
# end


@testset "BenesNetwork{$T}" for T in testtypes
    for _ in 1:100
        n = rand(1:bitsize(T))
        p = randperm(n)
        p′ = collect(1:bitsize(T))
        p′[1:n] = p′[p]
        net = BenesNetwork{T}(p)
        for _ in 1:100
            x = rand(T)
            @test bitvec(x)[p′] === bitvec(@inferred bitpermute(net, x))
        end
    end
end
