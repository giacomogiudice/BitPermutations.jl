using BenchmarkTools
using Plots
using Plots.Measures

bench_group = ["BenesNetwork", "GRPNetwork", "AVXCopyGather", "Bits", "BitVector"]
base_types = (UInt8, UInt16, UInt32, UInt64, UInt128)
avx_types = (UInt16, UInt32, UInt64)
N = 10_000

fname = nothing
if !isempty(ARGS)
    fname = first(ARGS)
    results, = BenchmarkTools.load(fname)
end

gap = 2
tick_locations = eachindex(base_types) * (length(bench_group) + gap) .+ length(bench_group) / 2
marker_size = 4

get_time(group, key) = haskey(group, key) ? time(group[key]) : NaN

p_single = plot(; title="Single repeated permutation", legend=:bottomright)
spans = (eachindex(base_types) .+ 1) * (length(bench_group) + gap) .- gap / 2
vspan!(spans; linecolor=:grey, fillcolor=:grey, alpha=0.05, label=nothing)
xs = [eachindex(base_types) * (length(bench_group) + gap) .+ k for k in eachindex(bench_group)]
ys = [[get_time(results["Single"][T], g) / N for T in base_types] for g in bench_group]
scatter!(
    xs, ys; seriescolor=permutedims(1:length(bench_group)), markersize=marker_size, labels=permutedims(bench_group)
)

p_broadcasted = plot(; labels=permutedims(bench_group), title="In-place broadcasted permutation", legend=:none)
spans = (eachindex(base_types) .+ 1) * (length(bench_group) + gap) .- gap / 2
vspan!(spans; linecolor=:grey, fillcolor=:grey, alpha=0.05, label=nothing)
xs = [eachindex(base_types) * (length(bench_group) + gap) .+ k for k in eachindex(bench_group)]
ys = [[get_time(results["Broadcasted"][T], g) / N for T in base_types] for g in bench_group]
scatter!(
    xs, ys; seriescolor=permutedims(1:length(bench_group)), markersize=marker_size, labels=permutedims(bench_group)
)

p = plot(
    p_single,
    p_broadcasted;
    layout=(1, 2),
    size=(920, 360),
    framestyle=:box,
    markerstrokewidth=0.0,
    grid=:y,
    minorgrid=:y,
    link=:both,
    yscale=:log10,
    yticks=(10) .^ (0:3),
    ylabel="Time per permutation [ns]",
    xticks=(tick_locations, string.(base_types)),
    xtick_direction=:none,
    foreground_color_legend=nothing,
    margin=10pt,
    guidefontsize=10,
    titlefontsize=10,
    tickfontsize=9,
    format=:svg,
)

display(p)
