using Pkg: Pkg
Pkg.activate(".")

using BenchmarkTools
using BenchmarkPlots, StatsPlots

include("utils.jl")

plotpath = "plots"
if !isdir(plotpath)
    mkdir(plotpath)
end

result = BenchmarkTools.load("bench.json")[1]

data = result["stateful"]["accessor"]
P = plot(data; yscale=:log10,ylim=_find_y_lims(data))
plot!(P; title="benchmark: particle stateful accessor")
savefig(joinpath(plotpath, "stateful_accessor.pdf"))

data = result["stateful"]["construction"]
P = plot(data; yscale=:log10,ylim=_find_y_lims(data))
plot!(P; title="benchmark: construction particle stateful")
savefig(joinpath(plotpath, "stateful_construction.pdf"))

data = result["phase space point"]["generate"]
P = plot(data; yscale=:log10,ylim=_find_y_lims(data))
plot!(P; title="benchmark: generate phase space points")
savefig(joinpath(plotpath, "phase_space_generation.pdf"))

data = result["phase space point"]["construction"]
P = plot(data; yscale=:log10,ylim=_find_y_lims(data))
plot!(P; title="benchmark: construction phase space points")
savefig(joinpath(plotpath, "phase_space_construction.pdf"))

data = result["phase space point"]["momentum"]
P = plot(data; yscale=:log10,ylim=_find_y_lims(data))
plot!(P; title="benchmark: phase space points momentum access")
savefig(joinpath(plotpath, "phase_space_momentum_access.pdf"))

data = result["phase space point"]["indexing"]
P = plot(data; yscale=:log10,ylim=_find_y_lims(data))
plot!(P; title="benchmark: phase space points indexingaccess")
savefig(joinpath(plotpath, "phase_space_momentum_indexing.pdf"))
