using QEDprocesses
using Test
using SafeTestsets

begin
    # modules
    @time @safetestset "propagators" begin
        include("propagators.jl")
    end
    @time @safetestset "cross section & probability" begin
        include("cross_sections.jl")
    end

    @time @safetestset "phase spaces" begin
        include("phase_spaces.jl")
    end

    # scattering processes
    include("processes/run_process_test.jl")
end
