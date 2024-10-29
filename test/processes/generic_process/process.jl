using QEDprocesses
using Random
using QEDcore

include("groundtruths.jl")

const RNG = MersenneTwister(77697185)

PARTICLES = [Electron(), Positron(), Photon()]
POLS = [PolX(), PolY(), AllPol(), SyncedPolarization(1)]
SPINS = [SpinUp(), SpinDown(), AllSpin(), SyncedSpin(1)]
POL_AND_SPIN_COMBINATIONS = Iterators.product(SPINS, POLS, SPINS, POLS)
BUF = IOBuffer()

@testset "constructor" begin
    @testset "default" begin
        proc = ScatteringProcess((Photon(), Electron()), (Photon(), Electron()))

        @test QEDprocesses.spin_pols(proc, Incoming())[1] == AllPol()
        @test QEDprocesses.spin_pols(proc, Incoming())[2] == AllSpin()
        @test QEDprocesses.spin_pols(proc, Outgoing())[1] == AllPol()
        @test QEDprocesses.spin_pols(proc, Outgoing())[2] == AllSpin()

        print(BUF, proc)
        @test String(take!(BUF)) == "generic QED process \"ke -> ke\""

        show(BUF, MIME"text/plain"(), proc)
        @test String(take!(BUF)) ==
            "generic QED process\n    incoming: photon ($(AllPol())), electron ($(AllSpin()))\n    outgoing: photon ($(AllPol())), electron ($(AllSpin()))\n"

        @test isphysical(proc, PerturbativeQED())
    end

    @testset "all spins+pols" begin
        @testset "$in_spin, $in_pol, $out_spin, $out_pol" for (
            in_spin, in_pol, out_spin, out_pol
        ) in POL_AND_SPIN_COMBINATIONS
            proc = ScatteringProcess(
                (Photon(), Electron()),
                (Photon(), Electron()),
                (in_pol, in_spin),
                (out_pol, out_spin),
            )

            @test QEDprocesses.spin_pols(proc, Incoming())[1] == in_pol
            @test QEDprocesses.spin_pols(proc, Incoming())[2] == in_spin
            @test QEDprocesses.spin_pols(proc, Outgoing())[1] == out_pol
            @test QEDprocesses.spin_pols(proc, Outgoing())[2] == out_spin

            print(BUF, proc)
            @test String(take!(BUF)) == "generic QED process \"ke -> ke\""

            show(BUF, MIME"text/plain"(), proc)
            @test String(take!(BUF)) ==
                "generic QED process\n    incoming: photon ($(in_pol)), electron ($(in_spin))\n    outgoing: photon ($(out_pol)), electron ($(out_spin))\n"

            @test isphysical(proc, PerturbativeQED())
        end
    end

    @testset "invalid types passed" begin
        struct INVALID_PARTICLE end

        @test_throws MethodError ScatteringProcess(
            (INVALID_PARTICLE(), Electron()), (Photon(), Electron())
        )

        @test_throws MethodError ScatteringProcess(
            (Photon(), Electron()), (Photon(), INVALID_PARTICLE())
        )
    end

    @testset "incompatible spin/pols" begin
        @test_throws InvalidInputError(
            "particle \"electron\" is a fermion and should have a spin, but has \"all polarizations\"",
        ) ScatteringProcess(
            (Electron(), Photon()),
            (Electron(), Photon()),
            (AllPol(), AllPol()),
            (AllSpin(), AllPol()),
        )
        @test_throws InvalidInputError(
            "particle \"photon\" is a boson and should have a polarization, but has \"all spins\"",
        ) ScatteringProcess(
            (Electron(), Photon()),
            (Electron(), Photon()),
            (AllSpin(), AllSpin()),
            (AllSpin(), AllPol()),
        )
    end

    @testset "incompatible number of spins/pols" begin
        IN_PARTICLES = Tuple(rand(RNG, PARTICLES, 2))
        OUT_PARTICLES = Tuple(rand(RNG, PARTICLES, 2))
        @testset "2 particles, 1 spin/pol" begin
            @test_throws MethodError ScatteringProcess(
                IN_PARTICLES,
                OUT_PARTICLES,
                _random_spin_pols(RNG, IN_PARTICLES)[1:1],
                _random_spin_pols(RNG, OUT_PARTICLES),
            )
            @test_throws MethodError ScatteringProcess(
                IN_PARTICLES,
                OUT_PARTICLES,
                _random_spin_pols(RNG, IN_PARTICLES),
                _random_spin_pols(RNG, OUT_PARTICLES)[1:1],
            )
        end

        @testset "2 particles, 3 spin/pols" begin
            @test_throws MethodError ScatteringProcess(
                IN_PARTICLES,
                OUT_PARTICLES,
                (_random_spin_pols(RNG, IN_PARTICLES)..., AllPol()),
                _random_spin_pols(RNG, OUT_PARTICLES),
            )
            @test_throws MethodError ScatteringProcess(
                IN_PARTICLES,
                OUT_PARTICLES,
                _random_spin_pols(RNG, IN_PARTICLES),
                (_random_spin_pols(RNG, OUT_PARTICLES)..., AllPol()),
            )
        end
    end
end

@testset "particle content" begin
    @testset "$n -> $m processes" for (n, m) in Base.product((2, 4), (3, 5))
        IN_PARTICLES = Tuple(rand(RNG, PARTICLES, n))
        OUT_PARTICLES = Tuple(rand(RNG, PARTICLES, m))
        proc = ScatteringProcess(IN_PARTICLES, OUT_PARTICLES)
        @testset "process $(proc)" begin
            @test incoming_particles(proc) == IN_PARTICLES
            @test outgoing_particles(proc) == OUT_PARTICLES
            @test number_incoming_particles(proc) == n
            @test number_outgoing_particles(proc) == m
            @test incoming_spin_pols(proc) == _groundtruth_spin_pols(IN_PARTICLES)
            @test outgoing_spin_pols(proc) == _groundtruth_spin_pols(OUT_PARTICLES)

            @test isphysical(proc, PerturbativeQED()) ==
                _groundtruth_is_physical(proc, PerturbativeQED())
        end

        IN_SPIN_POLS = _random_spin_pols(RNG, IN_PARTICLES)
        OUT_SPIN_POLS = _random_spin_pols(RNG, OUT_PARTICLES)
        proc = ScatteringProcess(IN_PARTICLES, OUT_PARTICLES, IN_SPIN_POLS, OUT_SPIN_POLS)
        @testset "process $(proc) with set spins/pols" begin
            @test incoming_particles(proc) == IN_PARTICLES
            @test outgoing_particles(proc) == OUT_PARTICLES
            @test number_incoming_particles(proc) == n
            @test number_outgoing_particles(proc) == m
            @test incoming_spin_pols(proc) == IN_SPIN_POLS
            @test outgoing_spin_pols(proc) == OUT_SPIN_POLS

            @test isphysical(proc, PerturbativeQED()) ==
                _groundtruth_is_physical(proc, PerturbativeQED())
        end
    end
end
