using Test
using SafeTestsets

include("utils.jl")

# check if we run CPU tests (yes by default)
cpu_tests = _is_test_platform_active(["CI_QED_TEST_CPU", "TEST_CPU"], true)

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
