
using Random
using QEDprocesses

RNG = MersenneTwister(137137137)

RND_MASS = rand(RNG)
RND_CHARGE = rand(RNG)

struct TestParticle <: AbstractParticle end
QEDprocesses.mass(::TestParticle) = RND_MASS 
QEDprocesses.charge(::TestParticle) = RND_CHARGE

@testset "default interface" begin
    @test !is_fermion(TestParticle())
    @test !is_boson(TestParticle())
    @test is_particle(TestParticle())
    @test !is_anti_particle(TestParticle())
end

@testset "hard interface" begin
    @test mass(TestParticle()) == RND_MASS
    @test charge(TestParticle()) == RND_CHARGE
end

struct TestParticle_FAIL <: AbstractParticle end
@testset "interface fail" begin
    @test_throws MethodError charge(TestParticle_FAIL())
    @test_throws MethodError mass(TestParticle_FAIL())
end
