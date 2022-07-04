# BitPermutations.jl
[![CI][ci-img]][ci-url]
[![codecov][codecov-img]][codecov-url]

[ci-img]: https://github.com/giacomogiudice/BitPermutations.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/giacomogiudice/BitPermutations.jl/actions?query=workflow%3ACI

[codecov-img]: https://codecov.io/gh/giacomogiudice/BitPermutations.jl/branch/master/graph/badge.svg?token=G71Y9FQH6K
[codecov-url]: https://codecov.io/gh/giacomogiudice/BitPermutations.jl?token=G71Y9FQH6K

Efficient routines for repeated bit permutations.

## Introduction
Permutations of *n* bits can be performed in *O(log(n))* operations which reshuffle the individual bits in parallel.
To have each reshuffling layer efficient, they have to be chosen from sets of operations which are efficient on the CPU.
Precomputing these operations is slighly non-trivial, so this package may be useful only if you need to compute the application of a given permutation to a large number of words.

## Usage
To define a permutation of the bits in type `T`, construct a `BitPermutation{T}`.

Here is an example with `T = UInt8`.
Bits are ordered from LSB to MSB.
```julia
using BitPermutations

v = [2,6,5,8,4,7,1,3]

p = BitPermutation{UInt8}(v)
```
we can then apply the permutation on any bitstring of type `UInt8` by using `bitpermute`.
If a functional approach is your jam, you can equivalently call the `BitPermutation` directly
```julia
x = 0x08          # 0 0 0 1 0 0 0 0 (LSB -> MSB)

bitpermute(x, p)  # 0 0 0 0 1 0 0 0, or 0x10

p(x)              # idem
```

To inspect the result, we can use `bitstring`, or we can use a `MBitVector` (defined by the package).
It is basically a faster `BitVector`, since its size if fixed (but is mutable).
```julia
using BitPermutations: MBitVector

mb = MBitVector(x)

[mb[v], MBitVector(bitpermute(x, p))]
```

The neat thing of the underlying network is that the inverse permutation can be computed at the same cost.
This can be performed using `invbitpermute` or calling the adjoint permutation
```julia
invbitpermute(x, p) === p'(x)

[mb[invperm(v)], MBitVector(invbitpermute(x, p))]
```

Internally, a `Vector` of bit masks and shifts are stored and then applied sequentially at each call to `bitpermute`.
If you know the permutation beforehand, you can use the `@bitpermutation` macro.
It will minimize the overhead of doing the reduction and compile down to a set of bit moving instructions.
```julia
p = @bitpermutation UInt32 (2, 6, 5, 8, 4, 7, 1, 3)

x = rand(UInt32)

p(x), p'(x)
```

## Benchmarks
As discussed later on, choosing types `UInt32` or `UInt64` can lead to significant speedups on processors which support BMI2 instructions.
Here are some benchmarks on an Intel Haswell processor.
```julia
using BitPermutations, Random, BenchmarkTools
using BitPermutations: MBitVector

T = UInt64
p = shuffle!(collect(1:bitsize(T)))
x = rand(T)

benes = BitPermutation{T}(p; type=BenesNetwork)

grp = BitPermutation{T}(p; type=GRPNetwork)
```
Notice that random permutations are typically the worst-case scenario, as the networks are usually the deepest.
Nevertheless
```
julia> @btime bitpermute($x, $benes);
  30.981 ns (0 allocations: 0 bytes)

julia> @btime bitpermute($x, $grp);
  10.872 ns (0 allocations: 0 bytes)
```
If you don't have hardware acceleration, the `GRPNetwork` permutation will still work but will be very slow.
Let's now compare to naively permuting `MBitVector`s and `BitVector`s from `Base`.
```julia 
mv = MBitVector(x)
bv = BitVector(mv)
```
They are around 5x (15x) and 10x (30x) slower compared to the `BenesNetwork` (`GRPNetwork`). 
```
julia> @btime permute!($mv, $p);
  166.468 ns (1 allocation: 576 bytes)

julia> @btime permute!($bv, $p);
  303.579 ns (1 allocation: 576 bytes) 
```
If your permutation is not random, it is likely to have even larger speedups as the networks will have fewer layers.

## Details
For a more in-depth explanation, the wonderful [https://programming.sirrida.de/bit_perm.html](https://programming.sirrida.de/bit_perm.html) is well worth reading.

Two different ways are performing the permutation are implemented: rearranged **Beneš networks** and **GRP networks**.
The latter is only faster on CPUs which support the [BMI2](https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set) instruction set.
Hence, the permutation is constructed using a `BenesNetwork{T}`, unless `T<:Union{UInt32,UInt64}` and BMI2 instructions are supported, in which case it uses a `GRPNetwork{T}`.
BMI2 intrinsics can be disabled by setting `ENV["PB_USE_BMI2"] = false` before loading the package or setting
```bash
export PB_USE_BMI2=false
```
before launching Julia.

### Beneš networks
A Beneš network is a series of **butterfly** or **delta-swap** operations, in which each node can be swapped or not with the corresponding node shifted by a fixed amount *δ*.
These operations are arranged in pairs, in a nested fashion, with the shifts chosen to be powers of two.

<img src="./network.png" alt="Beneš network" width=320/>

Above is an example of a network with 8 nodes, with 3 different stages (pairs of reshuffling) which have as *δ* respectively 4, 2, 1.
The two innermost operations can always be fused into a single one.

Somewhat remarkably, in this way one can perform any arbitrary permutation.
It is also relatively efficient, as each delta-swap should take around 6 cycles on modern processors.
The construction of the masks for the swaps is explained in: Donald E. Knuth, *The Art of Computer Programming*, Vol. 4, Fascicle 1A, ([Addison-Wesley, Stanford, 2008](https://www-cs-faculty.stanford.edu/~knuth/taocp.html)), available as pre-fascicle [here](http://www-cs-faculty.stanford.edu/%7Eknuth/fasc1a.ps.gz).
The general idea is to view each stage (pair of masks) as a repartition of the nodes into two sets.
We can then construct a graph in which each edge corresponds to the constraint of having two nodes in different partitions, both on the input and output side.
If we can 2-color the graph, we can then route each node to their corresponding partition, or color. 
Fortunately, the graph is bipartite and it is very easy to do so: we iterate through the cycles of the graph and assign alternating colors to the nodes we visit.

One would wish to optimize the network in such a way that most masks are trivial (i. e. no swaps).
Unfortunately I do not know of any other way that exhaustively checking all possible *log(n)!* arrangements.
This can be disabled by setting the keyword argument `rearrage=false` to the constructor.

### GRP networks
GRP networks work in a similar way to Beneš network, except that each layer is a different reshuffling, known as *GRP* or [sheeps-and-goats](https://programming.sirrida.de/bit_perm.html#sag) (see also TAOCP).
GRP networks are typically shallower, but the reshuffling operation is only efficient if specific instructions are available in hardware, as they can be performed in 8 cycles.

The [PEXT/PDEP](https://www.chessprogramming.org/BMI2) instructions used for the GRP reshuffling is supported by Intel starting from the Haswell architecture (released in 2013) and by AMD from the Zen 3 architecture (released in 2020).
On older AMD architectures, PEXT/PDEP is implemented in microcode and are reportedly slower.
On such machines you may want to experiment which method is faster and possibly disable calls to BMI2 with `ENV["PB_USE_BMI2"] = false`.

Fallback operations are implemented but are typically much slower then butterfly operations.
The construction of the masks follows the algorithm in: R. Lee, Z. Shi, X. Yang, *Efficient permutation instructions for fast software cryptography*, [IEEE Micro](https://doi.org/10.1109/40.977759) (2001); which is well explained [here](https://programming.sirrida.de/bit_perm.html#lee_sag).

# Enhancements
Several improvements could be made.
Here I just list the first ones off the top of my head:

- **Preconditioning** A simple way of reducing the depth of the network is to try cheap operations like `bitreverse` and `bitshift` before or after the network to try to save masks. This is what is done [here](https://programming.sirrida.de/calcperm.php).
- **Lookup tables** Small permutations could benefit from doing sub-permutations with a precomputed table. One could use `pext` to extract say 8 bits at a time, look up their permutation in a table of 256 elements, and join the results with `|`. This approach is fast but scales linearly with size, both in time and space, so it is interesting for permutations on `UInt8`s, `UInt16`s and possibly `UInt32`s.
- **PSHUFB** The [PSHUFB](https://www.chessprogramming.org/SSSE3#PSHUFB) instruction is part of SSE3 and can perform arbitrary byte reshufflings. It could be used to compress a bunch of layers or combined with lookup tables for some very promising speedups.
- **Rearrangement** Finding the best possible arrangement of the stages of a Beneš network (such that the most masks are trivial) is currently done by exhaustive search. It is likely there is a better way of constructing the masks, but I am not aware of any.
- **Better fallbacks** Faster software fallbacks for `pext/pdep` exist, like [zp7](https://github.com/zwegner/zp7).
- **Code refactoring** It should be possible to take a more functional approach and define a permutation as a series of transformations `T` ↦ `T`, but I'm not sure how to do that while preserving type inference and performance. This would allow for more generic algorithms and extensions.

# Compatibility
This package is compatible with Julia 1.5 and above.
