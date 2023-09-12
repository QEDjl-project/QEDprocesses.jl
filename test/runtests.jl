using QEDprocesses
using Test 
using SafeTestsets

begin
    @safetestset "particle interface" begin include("interfaces/particle_interface.jl") end
end
