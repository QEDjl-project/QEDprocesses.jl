@time @safetestset "generic process" begin
    include("generic_process/process.jl")
end

@time @safetestset "general one photon compton" begin
    include("one_photon_compton/process.jl")
end

@testset "perturbative one photon compton" begin
    @safetestset "kinematics" begin
        include("one_photon_compton/perturbative/kinematics.jl")
    end

    @safetestset "cross section" begin
        include("one_photon_compton/perturbative/cross_section.jl")
    end
end
