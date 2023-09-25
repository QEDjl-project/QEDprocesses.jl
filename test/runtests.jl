using QEDprocesses
using Test
using SafeTestsets

begin
    # Interfaces
    @time @safetestset "particle interface" begin
        include("interfaces/particle_interface.jl")
    end
    @time @safetestset "model interface" begin
        include("interfaces/model_interface.jl")
    end
    @time @safetestset "process interface" begin
        include("interfaces/process_interface.jl")
    end

    # modules
    @time @safetestset "particles types" begin
        include("particle_types.jl")
    end
end
