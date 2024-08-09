using Test
using SafeTestsets

begin
    # scattering processes
    include("processes/run_process_test.jl")

    include("gpu/process_interface.jl")
end
