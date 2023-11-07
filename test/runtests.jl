using QEDprocesses
using Test
using SafeTestsets

begin
    # Interfaces
    @time @safetestset "model interface" begin
        include("interfaces/model_interface.jl")
    end
    @time @safetestset "process interface" begin
        include("interfaces/process_interface.jl")
    end
    @time @safetestset "computation setup interface" begin
        include("interfaces/setup_interface.jl")
    end

    # modules
    @time @safetestset "propagators" begin
        include("propagators.jl")
    end
end
