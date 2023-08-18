using BitPermutations
using BitPermutations: USE_BMI2, USE_AVX512
using BitPermutations: deltaswap, grpswap, invgrpswap
using BitIntegers
using Random
using Test

# Set seed for reproducibility
Random.seed!(42)

# For permutations it is more natural to think of bits from left to right
bitstr(x) = reverse(bitstring(x))

# Check if using instrinsics
@info """
ENV[\"BIT_PERMUTATIONS_USE_INTRINSICS\"] = $(get(ENV, "BIT_PERMUTATIONS_USE_INTRINSICS", nothing))
USE_BMI2 = $(USE_BMI2)
USE_AVX512 = $(USE_AVX512)
"""

# Test types from Base
base_types = (UInt8, UInt16, UInt32, UInt64, UInt128)

# Types defined in BitIntegers
custom_types = (UInt256, UInt512, UInt1024)

@testset "Bits{$T}" for T in (base_types..., custom_types...)
    x = rand(T)
    @test (@inferred bitsize(x)) === 8 * sizeof(x)

    # Construction
    v = @inferred Bits(x)
    @test @inferred Bits{T}(x for x in v) == v
    @test @inferred all(((i, x),) -> v[i] === x, enumerate(v))
    bitstr(x) == string(v)

    # Comparison
    w = Bits(v)
    @test (@inferred w == v) && (@inferred isequal(w, v))
    v₀ = Bits(zero(T))
    @test (v ≥ v₀) && !any(v₀) && !all(v₀)
    @test any(v) ? (v₀ < v) : (v === v₀)
    @test (@inferred sum(v)) === count_ones(x) && (@inferred sum(!, v)) === count_zeros(x)

    # Broadcasting logic
    y = rand(T)
    w = Bits(y)
    @test convert(T, .~v) === ~x
    @test convert(T, v .& w) === x & y
    @test convert(T, v .| w) === x | y
    @test convert(T, v .⊻ w) === x ⊻ y
    @test convert(T, .!v) === ~x

    # Bitrotate not defined for BitIntegers
    if T in base_types
        i = rand(axes(x, 1))
        @test convert(T, circshift(v, i)) === bitrotate(x, i)
    end

    # Assignment
    for _ in 1:20
        i = rand(axes(x, 1))
        b = rand(Bool)
        v[i] = b
        @test v[i] === b
    end
end

@testset "Primitives for type $T" for T in (base_types..., custom_types...)
    # Simple deltaswaps
    a, b, c, d, e = T(1), T(2), T(3), T(5), T(12)
    @test [bitstr(x)[1:4] for x in (a, b, c, d, e)] == ["1000", "0100", "1100", "1010", "0011"]
    @test (@inferred deltaswap(a, a, 1)) === b     # 1̄000... -> 0100...
    @test (@inferred deltaswap(d, b, 1)) === c     # 10̄10... -> 1100...
    @test (@inferred deltaswap(c, c, 2)) === e     # 1̄1̄00... -> 0011...

    # Compare with array version
    for _ in 1:20
        x = rand(T)
        s = rand(1:4)
        m = rand(T)
        m = (m & ~(m << s)) >> s
        v = @inferred deltaswap(Bits(x), Bits(m), s)
        v′ = @inferred Bits(deltaswap(x, m, s))
        @test v == v′
        @test (@inferred deltaswap(deltaswap(x, m, s), m, s)) === x
    end

    # GRP swaps
    a, b, c, d, e = T(1), T(2), T(3), T(6), T(14)
    @test [bitstr(x)[1:4] for x in (a, b, c, d, e)] == ["1000", "0100", "1100", "0110", "0111"]
    @test (@inferred grpswap(a, ~c)) === a     # 1̄000... -> 1000...
    @test (@inferred grpswap(b, ~b)) === a     # 1̄000... -> 1000...
    @test (@inferred grpswap(d, ~e)) === c     # 1̄11̄1... -> 1100...

    # Compare with array version
    for _ in 1:20
        x = rand(T)
        m = rand(T)
        v = @inferred grpswap(Bits(x), Bits(m))
        v′ = @inferred Bits(grpswap(x, m))
        @test v == v′
        @test (@inferred invgrpswap(grpswap(x, m), m)) === x
        v′ = @inferred invgrpswap(v, Bits(m))
        @test v′ == Bits(x)
    end
end

@testset "BenesNetwork{$T}" for T in base_types
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

@testset "GRPNetwork{$T}" for T in base_types
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

@testset "AVXCopyGather{$T}" for T in base_types
    for _ in 1:10
        n = rand(1:bitsize(T))
        p = randperm(n)
        p′ = collect(1:bitsize(T))
        p′[1:n] = p′[p]
        net = AVXCopyGather{T}(p)
        for _ in 1:100
            x = rand(T)
            @test bitstr(x)[p′] === bitstr(@inferred bitpermute(x, net))
            @test bitstr(x)[invperm(p′)] === bitstr(@inferred invbitpermute(x, net))
        end
    end
end

@testset "BitPermutation{$T}" for T in base_types
    p₀ = [2, 6, 5, 8, 4, 7, 1, 3]

    backends = [BenesNetwork, GRPNetwork, AVXCopyGather]

    @testset "$backend_type" for backend_type in backends
        P = BitPermutation{T}(p₀; type=backend_type)

        # Make sure printing returns something
        buf = IOBuffer()
        show(buf, P)
        @test !isempty(String(take!(buf)))
        show(buf, "text/plain", P)
        @test !isempty(String(take!(buf)))

        # Test permutation properties
        @test (@inferred order(P)) === 4
        @test collect(cycles(P)) == [[1, 2, 6, 7], [3, 5, 4, 8]]
        @test isodd(P) === isodd(P') === false
        @test iseven(P) === iseven(P') === true
        @test sign(P) === sign(P') === 1

        # Test conversion to `Vector`
        @test Vector(P) == p₀

        # Test random permutations
        for _ in 1:10
            p = randperm(bitsize(T))
            P = BitPermutation{T}(p; type=backend_type)

            for _ in 1:20
                x = rand(T)
                v = Bits(x)
                y = convert(T, permute!(copy(v), p))
                ȳ = convert(T, invpermute!(copy(v), p))
                @test P(x) === bitpermute(x, P) === invbitpermute(x, P') === y
                @test P'(x) === invbitpermute(x, P) === bitpermute(x, P') === ȳ
            end
        end

        # Test permutations of arrays
        for _ in 1:10
            p = randperm(bitsize(T))
            P = BitPermutation{T}(p; type=backend_type)

            arr = rand(T, 1000)
            arr_copy = copy(arr)
            @test P.(arr) == bitpermute.(arr, P) == invbitpermute.(arr, P') == [bitpermute(x, P) for x in arr]
            @test P'.(arr) == invbitpermute.(arr, P) == bitpermute.(arr, P') == [invbitpermute(x, P) for x in arr]

            arr = copy(arr_copy)
            arr_permuted = P.(arr)
            @test arr_permuted == (@inferred bitpermute(arr, P))
            @test arr_permuted == (@inferred bitpermute!(arr, P))
            arr = copy(arr_copy)
            @test arr == (@inferred invbitpermute(arr_permuted, P))
            @test arr == (@inferred invbitpermute!(arr_permuted, P))

            arr = copy(arr_copy)
            arr_permuted = P'.(arr)
            @test arr_permuted == (@inferred bitpermute(arr, P'))
            @test arr_permuted == (@inferred bitpermute!(arr, P'))
            arr = copy(arr_copy)
            @test arr == (@inferred invbitpermute(arr_permuted, P'))
            @test arr == (@inferred invbitpermute!(arr_permuted, P'))
        end
    end
end

@testset "BitPermutation{$T}" for T in custom_types
    p = randperm(bitsize(T))
    # Rearranging takes too long 
    P = BitPermutation{T}(p; type=BenesNetwork, rearrange=false)

    for _ in 1:100
        x = rand(T)
        v = bitstr(x)
        @test bitstr(@inferred P(x)) == v[p]
        @test bitstr(@inferred P'(x)) == v[invperm(p)]
    end

    # Test permutations of arrays
    for _ in 1:10
        p = randperm(bitsize(T))
        P = BitPermutation{T}(p; type=BenesNetwork, rearrange=false)
        arr = rand(T, 1000)
        arr_copy = copy(arr)
        @test P.(arr) == bitpermute.(arr, P) == invbitpermute.(arr, P') == [bitpermute(x, P) for x in arr]
        @test P'.(arr) == invbitpermute.(arr, P) == bitpermute.(arr, P') == [invbitpermute(x, P) for x in arr]

        arr = copy(arr_copy)
        arr_permuted = P.(arr)
        @test arr_permuted == (@inferred bitpermute(arr, P))
        @test arr_permuted == (@inferred bitpermute!(arr, P))
        arr = copy(arr_copy)
        @test arr == (@inferred invbitpermute(arr_permuted, P))
        @test arr == (@inferred invbitpermute!(arr_permuted, P))

        arr = copy(arr_copy)
        arr_permuted = P'.(arr)
        @test arr_permuted == (@inferred bitpermute(arr, P'))
        @test arr_permuted == (@inferred bitpermute!(arr, P'))
        arr = copy(arr_copy)
        @test arr == (@inferred invbitpermute(arr_permuted, P'))
        @test arr == (@inferred invbitpermute!(arr_permuted, P'))
    end
end
