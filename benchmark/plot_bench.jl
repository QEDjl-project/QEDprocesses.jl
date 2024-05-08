using Pkg: Pkg
Pkg.activate(".")

using BenchmarkTools
using BenchmarkPlots, StatsPlots

plotpath = "plots"
if !isdir(plotpath)
    mkdir(plotpath)
end

result = BenchmarkTools.load("bench.json")[1]

P = plot(result["stateful"]["accessor"]; yscale=:log10, ylim=(1e-1, 100))
plot!(P; title="benchmark: particle stateful accessor")
savefig(joinpath(plotpath, "stateful_accessor.pdf"))

P = plot(result["stateful"]["construction"]; yscale=:log10, ylim=(1, 1e10))
plot!(P; title="benchmark: construction particle stateful")
savefig(joinpath(plotpath, "stateful_construction.pdf"))

P = plot(result["phase space point"]["generate"]; yscale=:log10, ylim=(1e2, 1e5))
plot!(P; title="benchmark: generate phase space points")
savefig(joinpath(plotpath, "phase_space_generation.pdf"))

P = plot(result["phase space point"]["construction"]; yscale=:log10, ylim=(1e1, 1e3))
plot!(P; title="benchmark: construction phase space points")
savefig(joinpath(plotpath, "phase_space_construction.pdf"))

P = plot(result["phase space point"]["momentum"]; yscale=:log10, ylim=(1e0, 1e3))
plot!(P; title="benchmark: phase space points momentum access")
savefig(joinpath(plotpath, "phase_space_momentum_access.pdf"))

P = plot(result["phase space point"]["indexing"]; yscale=:log10, ylim=(1e0, 1e3))
plot!(P; title="benchmark: phase space points indexingaccess")
savefig(joinpath(plotpath, "phase_space_momentum_indexing.pdf"))
