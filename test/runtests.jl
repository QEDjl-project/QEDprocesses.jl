using QEDprocesses
using Test
using SafeTestsets

begin
    # modules
    @time @safetestset "cross section & probability" begin
        include("cross_sections.jl")
    end

    # scattering processes
    include("processes/run_process_test.jl")
end
