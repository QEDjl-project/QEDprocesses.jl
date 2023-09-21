using QEDprocesses

struct TestModel <: AbstractModelDefinition end
fundamental_interaction_type(::TestModel) = :test_interaction

struct TestModel_FAIL <: AbstractModelDefinition end

@testset "hard interface" begin
    @test fundamental_interaction_type(TestModel()) == :test_interaction
end

@testset "interface fail" begin
    @test_throws MethodError fundamental_interaction_type(TestModel_FAIL())
end
