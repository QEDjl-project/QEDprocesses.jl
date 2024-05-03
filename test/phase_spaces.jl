using Random
using QEDbase
using QEDprocesses

# can be removed when QEDbase exports them
import QEDbase.is_incoming, QEDbase.is_outgoing

@testset "Phasespace Point" begin
    DIRECTIONS = [Incoming(), Outgoing()]
    PARTICLES = [Electron(), Positron()] #=, Muon(), AntiMuon(), Tauon(), AntiTauon()=#
    SPINANDPOLS = [AllSpin(), SpinUp(), SpinDown(), AllPol(), PolX(), PolY()]

    for (particle, dir, spinorpol) in Iterators.product(PARTICLES, DIRECTIONS, SPINANDPOLS)
        momentum = rand(SFourMomentum)

        if (is_fermion(particle) && (spinorpol isa AbstractSpin)) ||
            (is_boson(particle) && (spinorpol isa AbstractPolarization))
            particle_stateful = ParticleStateful(momentum, particle, spinorpol, dir)

            @test particle_stateful.mom == momentum
            @test is_fermion(particle_stateful) == is_fermion(particle)
            @test is_boson(particle_stateful) == is_boson(particle)
            @test is_particle(particle_stateful) == is_particle(particle)
            @test is_anti_particle(particle_stateful) == is_anti_particle(particle)
            @test is_incoming(particle_stateful) == is_incoming(dir)
            @test is_outgoing(particle_stateful) == is_outgoing(dir)
        else
            @test_throws "Tried to construct a stateful $(is_fermion(particle) ? "fermion" : "boson") with a $(spinorpol) instead of a $(is_fermion(particle) ? "spin" : "polarization")" ParticleStateful(
                momentum, particle, spinorpol, dir
            )
        end
    end
end
