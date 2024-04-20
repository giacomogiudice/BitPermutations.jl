# Benchmarking

To gather the benchmarking data, run the following from the project root:

```
julia --project --threads=1 --optimize=3 --check-bounds=no benchmark/benchmarks.jl benchmark/data/results.json
```

To plot the results, run
```
julia benchmark/plots.jl benchmark/data/results.json
```
