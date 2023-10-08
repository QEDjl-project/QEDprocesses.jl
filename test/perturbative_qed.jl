using Random
using QEDbase
using QEDprocesses

@testset "interaction type" begin
    @test fundamental_interaction_type(PerturbativeQED()) == :electromagnetic
end
