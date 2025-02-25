
using QEDbase
using QEDcore
using QEDprocesses
using Random

const RNG = MersenneTwister(77697185)
const ATOL = eps()
const RTOL = sqrt(eps())

const PROC = Compton()
const MODEL = PerturbativeQED()

const OMEGAS = (1e-6 * rand(RNG), 1e-3 * rand(RNG), rand(RNG), 1e3 * rand(RNG))
const SQRT_S = (1.0, 1 + rand(RNG))

const COS_THETAS = [-1.0, 2 * rand(RNG) - 1, 0.0, 1.0]
const PHIS = [0, 2 * pi, rand(RNG) * 2 * pi]

@testset "in-phase-space layout" begin
    @testset "Compton rest frame" begin
        @testset "default" begin
            @test ComptonRestSystem() == ComptonRestSystem(Energy(2))
        end

        @testset "wrong rest system" begin
            @testset "$coord_type" for coord_type in (Energy, SpatialMagnitude)
                coord = coord_type(1)
                @test_throws ArgumentError ComptonRestSystem(coord)
            end

            @testset "no photon rapidity" begin
                coord = Rapidity(2)
                @test_throws ArgumentError ComptonRestSystem(coord)
            end
        end

        @testset "photon energy" begin
            in_psl = ComptonRestSystem(Energy(2))

            @testset "$om" for om in OMEGAS
                in_coords = (om,)
                in_psp = InPhaseSpacePoint(PROC, MODEL, in_psl, in_coords)

                test_P, test_K = momenta(in_psp, Incoming())

                @test isapprox(getMass2(test_K), mass(Photon())^2, atol=ATOL, rtol=RTOL)
                @test isapprox(getMass2(test_P), mass(Electron())^2, atol=ATOL, rtol=RTOL)
                @test isapprox(getE(test_K), om, atol=ATOL, rtol=RTOL)
                @test isapprox(getE(test_P), mass(Electron()), atol=ATOL, rtol=RTOL)
            end
        end

        @testset "photon magnitude" begin
            in_psl = ComptonRestSystem(SpatialMagnitude(2))

            # for photons: E === rho
            @testset "$om" for om in OMEGAS
                in_coords = (om,)
                in_psp = InPhaseSpacePoint(PROC, MODEL, in_psl, in_coords)

                test_P, test_K = momenta(in_psp, Incoming())

                @test isapprox(getMass2(test_K), mass(Photon())^2, atol=ATOL, rtol=RTOL)
                @test isapprox(getMass2(test_P), mass(Electron())^2, atol=ATOL, rtol=RTOL)
                @test isapprox(getE(test_K), om, atol=ATOL, rtol=RTOL)
                @test isapprox(getE(test_P), mass(Electron()), atol=ATOL, rtol=RTOL)
            end
        end

        @testset "cms energy" begin
            in_psl = ComptonRestSystem(CMSEnergy())

            @testset "$ss" for ss in SQRT_S
                in_coords = (ss,)
                in_psp = InPhaseSpacePoint(PROC, MODEL, in_psl, in_coords)

                test_P, test_K = momenta(in_psp, Incoming())

                @test isapprox(getMass2(test_K), mass(Photon())^2, atol=ATOL, rtol=RTOL)
                @test isapprox(getMass2(test_P), mass(Electron())^2, atol=ATOL, rtol=RTOL)
                @test isapprox(getMass(test_K + test_P), ss, atol=ATOL, rtol=RTOL)
                @test isapprox(getE(test_P), mass(Electron()), atol=ATOL, rtol=RTOL)
            end
        end
    end
end

@testset "out-phase-space-layout" begin
    @testset "Compton spherical system" begin
        in_psl = ComptonRestSystem(Energy(2))
        out_psl = ComptonSphericalLayout(in_psl)

        @testset "wrong in_psl" begin
            in_psl_fail = TwoBodyBeamSystem()
            @test_throws ArgumentError ComptonSphericalLayout(in_psl_fail)
        end

        @testset "$om, $cth, $phi" for (om, cth, phi) in
                                       Iterators.product(OMEGAS, COS_THETAS, PHIS)
            IN_COORDS = (om,)
            OUT_COORDS = (cth, phi)

            test_psp = PhaseSpacePoint(PROC, MODEL, out_psl, IN_COORDS, OUT_COORDS)
            @testset "mass shell" begin
                @testset "$dir $part" for (dir, part) in Iterators.product(
                    (Incoming(), Outgoing()), (Electron(), Photon())
                )
                    mom = momentum(test_psp, dir, part)
                    mom_square = mom * mom
                    mass_square = mass(part)^2

                    @test isapprox(mom_square, mass_square, atol=1e2 * ATOL, rtol=RTOL)
                end
            end

            @testset "energy-momentum conservation" begin
                in_moms = momenta(test_psp, Incoming())
                out_moms = momenta(test_psp, Outgoing())

                @test isapprox(sum(in_moms), sum(out_moms), atol=ATOL, rtol=RTOL)
            end
        end
    end
end
