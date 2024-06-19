using QEDprocesses
using Test
using SafeTestsets

begin
    # # Interfaces
    @time @safetestset "model interface" begin
        include("interfaces/model_interface.jl")
    end
    @time @safetestset "process interface" begin
        include("interfaces/process_interface.jl")
    end

    # TODO: remove after restructuring
    #@time @safetestset "computation setup interface" begin
    #    include("interfaces/setup_interface.jl")
    #end

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
