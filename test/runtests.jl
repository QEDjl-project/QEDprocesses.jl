using QEDprocesses
using Test 
using SafeTestsets

begin
    # Interfaces
    @time @safetestset "particle interface" begin include("interfaces/particle_interface.jl") end

    # modules
    @time @safetestset "particles types" begin include("particle_types.jl") end
end
