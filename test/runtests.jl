using Test
using SafeTestsets

begin
    # scattering processes
    include("processes/run_process_test.jl")

    @time @safetestset "GPU tests" begin
        include("gpu/process_interface.jl")
    end
end
