using QEDprocesses 
using Random
using QEDbase
using QEDprocesses

POLS = [PolX(), PolY(), AllPol()]
SPINS = [SpinUp(), SpinDown(), AllSpin()]
POL_AND_SPIN_COMBINATIONS = Iterators.product(SPINS, POLS, SPINS, POLS)
POL_COMBINATIONS = Iterators.product(POLS, POLS)

@testset "constructor" begin
    @testset "default" begin
        proc = Compton()
        @test QEDprocesses._spin_or_pol(proc, Photon(), Incoming()) == AllPol()
        @test QEDprocesses._spin_or_pol(proc, Electron(), Incoming()) == AllSpin()
        @test QEDprocesses._spin_or_pol(proc, Photon(), Outgoing()) == AllPol()
        @test QEDprocesses._spin_or_pol(proc, Electron(), Outgoing()) == AllSpin()
    end

    @testset "in_pol" begin
        @testset "$pol" for pol in POLS
            proc = Compton(pol)
            @test QEDprocesses._spin_or_pol(proc, Electron(), Incoming()) == AllSpin()
            @test QEDprocesses._spin_or_pol(proc, Photon(), Incoming()) == pol
            @test QEDprocesses._spin_or_pol(proc, Electron(), Outgoing()) == AllSpin()
            @test QEDprocesses._spin_or_pol(proc, Photon(), Outgoing()) == AllPol()
        end
    end

    @testset "in_pol+out_pol" begin
        @testset "$in_pol, $out_pol" for (in_pol, out_pol) in POL_COMBINATIONS
            proc = Compton(in_pol, out_pol)
            @test QEDprocesses._spin_or_pol(proc, Electron(), Incoming()) == AllSpin()
            @test QEDprocesses._spin_or_pol(proc, Photon(), Incoming()) == in_pol
            @test QEDprocesses._spin_or_pol(proc, Electron(), Outgoing()) == AllSpin()
            @test QEDprocesses._spin_or_pol(proc, Photon(), Outgoing()) == out_pol
        end
    end
    @testset "all spins+pols" begin
        @testset "$in_spin, $in_pol, $out_spin, $out_pol" for (
            in_spin, in_pol, out_spin, out_pol
        ) in POL_AND_SPIN_COMBINATIONS
            proc = Compton(in_spin, in_pol, out_spin, out_pol)
            @test QEDprocesses._spin_or_pol(proc, Electron(), Incoming()) == in_spin
            @test QEDprocesses._spin_or_pol(proc, Photon(), Incoming()) == in_pol
            @test QEDprocesses._spin_or_pol(proc, Electron(), Outgoing()) == out_spin
            @test QEDprocesses._spin_or_pol(proc, Photon(), Outgoing()) == out_pol
        end
    end
end
@testset "particle content" begin
    proc = Compton()
    @test incoming_particles(proc) == (Electron(), Photon())
    @test outgoing_particles(proc) == (Electron(), Photon())
    @test number_incoming_particles(proc) == 2
    @test number_outgoing_particles(proc) == 2
end
