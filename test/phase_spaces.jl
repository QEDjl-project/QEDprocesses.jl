using Random
using StaticArrays
using QEDbase
using QEDprocesses

# can be removed when QEDbase exports them
import QEDbase.is_incoming, QEDbase.is_outgoing

include("test_implementation/TestImplementation.jl")
TESTMODEL = TestImplementation.TestModel()
TESTPSDEF = TestImplementation.TestPhasespaceDef()

RNG = Random.MersenneTwister(727)

@testset "Stateful Particle" begin
    DIRECTIONS = [Incoming(), Outgoing()]
    PARTICLES = [Electron(), Positron()] #=, Muon(), AntiMuon(), Tauon(), AntiTauon()=#
    SPINANDPOLS = [AllSpin(), SpinUp(), SpinDown(), AllPol(), PolX(), PolY()]

    for (particle, dir, spinorpol) in Iterators.product(PARTICLES, DIRECTIONS, SPINANDPOLS)
        momentum = rand(RNG, SFourMomentum)

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
            if (VERSION >= v"1.8")
                # julia versions before 1.8 did not have support for regex matching in @test_throws
                @test_throws "Tried to construct a stateful $(is_fermion(particle) ? "fermion" : "boson") with a $(spinorpol) instead of a $(is_fermion(particle) ? "spin" : "polarization")" ParticleStateful(
                    momentum, particle, spinorpol, dir
                )
            end
            @test_throws InvalidInputError ParticleStateful(
                momentum, particle, spinorpol, dir
            )
        end
    end
end

@testset "Phasespace Point" begin
    in_el = ParticleStateful(rand(RNG, SFourMomentum), Electron(), AllSpin(), Incoming())
    in_ph = ParticleStateful(rand(RNG, SFourMomentum), Photon(), AllPol(), Incoming())
    out_el = ParticleStateful(rand(RNG, SFourMomentum), Electron(), AllSpin(), Outgoing())
    out_ph = ParticleStateful(rand(RNG, SFourMomentum), Photon(), AllPol(), Outgoing())

    in_particles_valid = SVector(in_el, in_ph)
    in_particles_invalid = SVector(in_el, out_ph)

    out_particles_valid = SVector(out_el, out_ph)
    out_particles_invalid = SVector(out_el, in_ph)

    model = TESTMODEL
    process = TestImplementation.TestProcess(
        SVector{2,AbstractParticle}(Electron(), Photon()),
        SVector{2,AbstractParticle}(Electron(), Photon()),
    )
    phasespace_def = TESTPSDEF

    PhaseSpacePoint(process, model, phasespace_def, in_particles_valid, out_particles_valid)

    if (VERSION >= v"1.8")
        # julia versions before 1.8 did not have support for regex matching in @test_throws
        @test_throws r"Stateful particle (.*) is given as an incoming particle but is outgoing" PhaseSpacePoint(
            process, model, phasespace_def, in_particles_invalid, out_particles_valid
        )

        @test_throws r"Stateful particle (.*) is given as an outgoing particle but is incoming" PhaseSpacePoint(
            process, model, phasespace_def, in_particles_valid, out_particles_invalid
        )

        @test_throws r"Process given particle type \((.*)Electron\(\)\) does not match stateful particle type \((.*)Photon\(\)\)" PhaseSpacePoint(
            process, model, phasespace_def, SVector(in_ph, in_el), out_particles_valid
        )

        @test_throws r"Process given particle type \((.*)Electron\(\)\) does not match stateful particle type \((.*)Photon\(\)\)" PhaseSpacePoint(
            process, model, phasespace_def, in_particles_valid, SVector(out_ph, out_el)
        )
    end

    @test_throws InvalidInputError PhaseSpacePoint(
        process, model, phasespace_def, in_particles_invalid, out_particles_valid
    )

    @test_throws InvalidInputError PhaseSpacePoint(
        process, model, phasespace_def, in_particles_valid, out_particles_invalid
    )

    @test_throws InvalidInputError PhaseSpacePoint(
        process, model, phasespace_def, SVector(in_ph, in_el), out_particles_valid
    )

    @test_throws InvalidInputError PhaseSpacePoint(
        process, model, phasespace_def, in_particles_valid, SVector(out_ph, out_el)
    )
end
