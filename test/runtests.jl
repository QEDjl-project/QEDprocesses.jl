using Test
using SafeTestsets

# check if we run CPU tests (yes by default)
cpu_tests = tryparse(Bool, get(ENV, "TEST_CPU", "1"))

if cpu_tests
    # scattering processes
    include("processes/run_process_test.jl")
else
    @info "Skipping CPU tests"
end

begin
    @time @safetestset "GPU testing" begin
        include("gpu/runtests.jl")
    end
end
