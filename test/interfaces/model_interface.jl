using QEDprocesses

include("../test_implementation/TestImplementation.jl")

@testset "hard interface" begin
    TESTMODEL = TestImplementation.TestModel()
    @test fundamental_interaction_type(TESTMODEL) == :test_interaction
end

@testset "interface fail" begin
    TESTMODEL_FAIL = TestImplementation.TestModel_FAIL()
    @test_throws MethodError fundamental_interaction_type(TESTMODEL_FAIL)
end
