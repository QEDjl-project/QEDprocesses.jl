

using QEDprocesses

@testset "default interface" begin
    struct TestParticle <: AbstractParticle end
    @test !is_fermion(TestParticle())
    @test !is_boson(TestParticle())
    @test is_particle(TestParticle())
    @test !is_anti_particle(TestParticle())
end

@testset "hard interface" begin
    struct TestParticle <: AbstractParticle end
    @test_throws Exception charge(TestParticle())
    @test_throws Exception mass(TestParticle())
end
