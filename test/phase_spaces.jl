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
    SPECIES = [Electron(), Positron()] #=, Muon(), AntiMuon(), Tauon(), AntiTauon()=#
    SPINANDPOLS = [AllSpin(), SpinUp(), SpinDown(), AllPol(), PolX(), PolY()]

    for (species, dir, spin_or_pol) in Iterators.product(SPECIES, DIRECTIONS, SPINANDPOLS)
        momentum = rand(RNG, SFourMomentum)

        if (is_fermion(species) && (spin_or_pol isa AbstractSpin)) ||
            (is_boson(species) && (spin_or_pol isa AbstractPolarization))
            particle_stateful = ParticleStateful(dir, species, momentum, spin_or_pol)

            @test particle_stateful.mom == momentum
            @test is_fermion(particle_stateful) == is_fermion(species)
            @test is_boson(particle_stateful) == is_boson(species)
            @test is_particle(particle_stateful) == is_particle(species)
            @test is_anti_particle(particle_stateful) == is_anti_particle(species)
            @test is_incoming(particle_stateful) == is_incoming(dir)
            @test is_outgoing(particle_stateful) == is_outgoing(dir)
            @test mass(particle_stateful) == mass(species)
            @test charge(particle_stateful) == charge(species)

            if (is_fermion(species))
                @test spin(particle_stateful) == spin_or_pol
                @test_throws MethodError polarization(particle_stateful)
            else
                @test polarization(particle_stateful) == spin_or_pol
                @test_throws MethodError spin(particle_stateful)
            end
        else
            if (VERSION >= v"1.8")
                # julia versions before 1.8 did not have support for regex matching in @test_throws
                @test_throws "MethodError: no method matching ParticleStateful" ParticleStateful(
                    dir, species, momentum, spin_or_pol
                )
            end
            @test_throws MethodError ParticleStateful(dir, species, momentum, spin_or_pol)
        end
    end
end

@testset "Phasespace Point" begin
    in_el_mom = rand(RNG, SFourMomentum)
    in_ph_mom = rand(RNG, SFourMomentum)
    out_el_mom = rand(RNG, SFourMomentum)
    out_ph_mom = rand(RNG, SFourMomentum)

    in_el = ParticleStateful(Incoming(), Electron(), in_el_mom)
    in_ph = ParticleStateful(Incoming(), Photon(), in_ph_mom)
    out_el = ParticleStateful(Outgoing(), Electron(), out_el_mom)
    out_ph = ParticleStateful(Outgoing(), Photon(), out_ph_mom)

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

    psp = PhaseSpacePoint(
        process, model, phasespace_def, in_particles_valid, out_particles_valid
    )

    @testset "Accessor" begin
        @test momentum(psp, Incoming(), 1) == in_el.mom
        @test momentum(psp, Incoming(), 2) == in_ph.mom
        @test momentum(psp, Outgoing(), 1) == out_el.mom
        @test momentum(psp, Outgoing(), 2) == out_ph.mom

        @test psp[Incoming(), 1] == in_el
        @test psp[Incoming(), 2] == in_ph
        @test psp[Outgoing(), 1] == out_el
        @test psp[Outgoing(), 2] == out_ph
    end

    @testset "Error handling" begin
        if (VERSION >= v"1.8")
            # julia versions before 1.8 did not have support for regex matching in @test_throws
            @test_throws r"stateful particle (.*) is given as an incoming particle but is outgoing" PhaseSpacePoint(
                process, model, phasespace_def, in_particles_invalid, out_particles_valid
            )

            @test_throws r"stateful particle (.*) is given as an outgoing particle but is incoming" PhaseSpacePoint(
                process, model, phasespace_def, in_particles_valid, out_particles_invalid
            )

            @test_throws r"process given particle species \((.*)Electron\(\)\) does not match stateful particle species \((.*)Photon\(\)\)" PhaseSpacePoint(
                process, model, phasespace_def, SVector(in_ph, in_el), out_particles_valid
            )

            @test_throws r"process given particle species \((.*)Electron\(\)\) does not match stateful particle species \((.*)Photon\(\)\)" PhaseSpacePoint(
                process, model, phasespace_def, in_particles_valid, SVector(out_ph, out_el)
            )
        end

        @test_throws BoundsError momentum(psp, Incoming(), -1)
        @test_throws BoundsError momentum(psp, Outgoing(), -1)
        @test_throws BoundsError momentum(psp, Incoming(), 4)
        @test_throws BoundsError momentum(psp, Outgoing(), 4)

        @test_throws BoundsError psp[Incoming(), -1]
        @test_throws BoundsError psp[Outgoing(), -1]
        @test_throws BoundsError psp[Incoming(), 4]
        @test_throws BoundsError psp[Outgoing(), 4]

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

    @testset "Generation" begin
        test_psp = generate_phase_space(
            process, model, phasespace_def, [in_el_mom, in_ph_mom], [out_el_mom, out_ph_mom]
        )

        @test test_psp.proc == process
        @test test_psp.model == model
        @test test_psp.ps_def == phasespace_def

        @test test_psp[Incoming(), 1] == in_el
        @test test_psp[Incoming(), 2] == in_ph
        @test test_psp[Outgoing(), 1] == out_el
        @test test_psp[Outgoing(), 2] == out_ph
    end
end
