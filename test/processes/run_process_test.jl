@time @safetestset "generic process" begin
    include("generic_process/process.jl")
end

@time @safetestset "general one photon compton" begin
    include("one_photon_compton/process.jl")
end

@time @safetestset "perturbative one photon compton" begin
    include("one_photon_compton/perturbative.jl")
end
