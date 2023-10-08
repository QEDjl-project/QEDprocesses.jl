using Random
using QEDbase
using QEDprocesses

@testset "constructors" begin
    @testset "default" begin
        proc = Compton()
        @test QEDprocesses._spin_or_pol(proc, Photon(), Incoming()) == AllPol()
        @test QEDprocesses._spin_or_pol(proc, Electron(), Incoming()) == AllSpin()
        @test QEDprocesses._spin_or_pol(proc, Photon(), Outgoing()) == AllPol()
        @test QEDprocesses._spin_or_pol(proc, Electron(), Outgoing()) == AllSpin()
    end

    @testset "non-default" begin
        proc = Compton((PolX(), PolY()), (SpinUp(), SpinDown()))
        @test QEDprocesses._spin_or_pol(proc, Photon(), Incoming()) == PolX()
        @test QEDprocesses._spin_or_pol(proc, Electron(), Incoming()) == SpinUp()
        @test QEDprocesses._spin_or_pol(proc, Photon(), Outgoing()) == PolY()
        @test QEDprocesses._spin_or_pol(proc, Electron(), Outgoing()) == SpinDown()
    end
end

@testset "scattering process interface" begin
    proc = Compton()

    @test incoming_particles(proc) == (Photon(), Electron())
    @test outgoing_particles(proc) == (Photon(), Electron())

    @test number_incoming_particles(proc) == 2
    @test number_outgoing_particles(proc) == 2
end

@testset "invalid inputs" begin
    model = PerturbativeQED()
    proc = Compton()

    momenta_2 = [zero(SFourMomentum) for _ = 1:2]
    momenta_3 = [zero(SFourMomentum) for _ = 1:3]

    momenta_2_2 = Matrix{SFourMomentum}(undef, 2, 2)
    momenta_3_2 = Matrix{SFourMomentum}(undef, 3, 2)

    for (valid, invalid) in
        Iterators.product([momenta_2, momenta_2_2], [momenta_3, momenta_3_2])
        # try compute single input with incorrect dimensions
        @test_throws DimensionMismatch differential_cross_section(
            proc,
            model,
            invalid,
            valid,
        )
        @test_throws "incoming" differential_cross_section(proc, model, invalid, valid)
        @test_throws DimensionMismatch differential_cross_section(
            proc,
            model,
            valid,
            invalid,
        )
        @test_throws "outgoing" differential_cross_section(proc, model, valid, invalid)
    end
end
