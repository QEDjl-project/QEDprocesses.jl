using Random
using QEDbase
using QEDprocesses

POLS = [PolX(), PolY(), AllPol()]
SPINS = [SpinUp(), SpinDown(), AllSpin()]
POL_AND_SPIN_COMBINATIONS = Iterators.product(POLS, POLS, SPINS, SPINS)

@testset "Default Compton constructor" begin
    proc = Compton()
    @test QEDprocesses._spin_or_pol(proc, Photon(), Incoming()) == AllPol()
    @test QEDprocesses._spin_or_pol(proc, Electron(), Incoming()) == AllSpin()
    @test QEDprocesses._spin_or_pol(proc, Photon(), Outgoing()) == AllPol()
    @test QEDprocesses._spin_or_pol(proc, Electron(), Outgoing()) == AllSpin()
end

@testset "Compton with ($pol1, $pol2), ($spin1, $spin2)" for (pol1, pol2, spin1, spin2) in
                                                             POL_AND_SPIN_COMBINATIONS
    proc = nothing
    @testset "constructor" begin
        proc = Compton((pol1, pol2), (spin1, spin2))
        @test QEDprocesses._spin_or_pol(proc, Photon(), Incoming()) == pol1
        @test QEDprocesses._spin_or_pol(proc, Electron(), Incoming()) == spin1
        @test QEDprocesses._spin_or_pol(proc, Photon(), Outgoing()) == pol2
        @test QEDprocesses._spin_or_pol(proc, Electron(), Outgoing()) == spin2
    end

    @testset "scattering process interface" begin
        @test incoming_particles(proc) == (Photon(), Electron())
        @test outgoing_particles(proc) == (Photon(), Electron())
        @test number_incoming_particles(proc) == 2
        @test number_outgoing_particles(proc) == 2
    end

    model = PerturbativeQED()
    proc = Compton((pol1, pol2), (spin1, spin2))

    @testset "invalid inputs" begin
        momenta_2 = [zero(SFourMomentum) for _ in 1:2]
        momenta_3 = [zero(SFourMomentum) for _ in 1:3]

        momenta_2_2 = Matrix{SFourMomentum}(undef, 2, 2)
        momenta_3_2 = Matrix{SFourMomentum}(undef, 3, 2)

        for (valid, invalid) in
            Iterators.product([momenta_2, momenta_2_2], [momenta_3, momenta_3_2])
            # try compute single input with incorrect dimensions
            @test_throws DimensionMismatch differential_cross_section(
                proc, model, invalid, valid
            )
            @test_throws "incoming" differential_cross_section(proc, model, invalid, valid)
            @test_throws DimensionMismatch differential_cross_section(
                proc, model, valid, invalid
            )
            @test_throws "outgoing" differential_cross_section(proc, model, valid, invalid)
        end
    end

    #= See https://github.com/QEDjl-project/QEDbase.jl/issues/36
    @testset "valid inputs" begin
        momenta_in = [SFourMomentum(1.0, 0.0, 0.0, 0.0), SFourMomentum(1.0, 0.0, 0.0, 0.0)]
        momenta_out = [SFourMomentum(1.0, 0.0, 0.0, 0.0), SFourMomentum(1.0, 0.0, 0.0, 0.0)]

        differential_cross_section(proc, model, momenta_in, momenta_out)
    end
    =#
end
