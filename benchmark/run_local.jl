using Pkg: Pkg

Pkg.activate(".")
Pkg.develop(; path="./..")

include("benchmark.jl")

tune!(SUITE)
result = run(SUITE; verbose=true)

BenchmarkTools.save("bench.json", result)
