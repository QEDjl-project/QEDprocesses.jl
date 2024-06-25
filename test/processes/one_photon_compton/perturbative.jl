
using QEDbase: QEDbase
using QEDcore
using QEDprocesses
using Random
using StaticArrays
using QuadGK

const RNG = MersenneTwister(77697185)
const ATOL = eps()
const RTOL = sqrt(eps())

include("groundtruths.jl")

const MODEL = PerturbativeQED()
const PS_DEF = PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame())
const OMEGAS = (1e-6 * rand(RNG), 1e-3 * rand(RNG), rand(RNG), 1e3 * rand(RNG))

const COS_THETAS = [-1.0, 2 * rand(RNG) - 1, 0.0, 1.0]
const PHIS = [0, 2 * pi, rand(RNG) * 2 * pi]

@testset "pretty-printing" begin
    buf = IOBuffer()
    print(buf, MODEL)
    @test String(take!(buf)) == "perturbative QED"

    show(buf, MIME"text/plain"(), MODEL)
    @test String(take!(buf)) == "perturbative QED"
end

@testset "perturbative kinematics" begin
    PROC = Compton()
    @testset "momentum generation" begin
        @testset "$om, $cth, $phi" for (om, cth, phi) in
                                       Iterators.product(OMEGAS, COS_THETAS, PHIS)
            IN_COORDS = (om,)
            OUT_COORDS = (cth, phi)
            IN_PS, OUT_PS = QEDprocesses._generate_momenta(
                PROC, MODEL, PS_DEF, IN_COORDS, OUT_COORDS
            )
            in_mom_square = getMass2.(IN_PS)
            out_mom_square = getMass2.(OUT_PS)
            in_masses = mass.(incoming_particles(PROC)) .^ 2
            out_masses = mass.(outgoing_particles(PROC)) .^ 2

            # we need a larger ATOL than eps() here because the error is accumulated over several additions
            @test all(isapprox.(in_mom_square, in_masses, atol=4 * ATOL, rtol=RTOL))
            @test all(isapprox.(out_mom_square, out_masses, atol=4 * ATOL, rtol=RTOL))
        end
    end
end

@testset "perturbative" begin
    @testset "$omega" for omega in OMEGAS
        @testset "differential cross section" begin
            @testset "spin and pol summed" begin
                PROC = Compton()

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = (omega,)
                    OUT_COORDS = (cos_theta, phi)
                    IN_PS, OUT_PS = QEDprocesses._generate_momenta(
                        PROC, MODEL, PS_DEF, IN_COORDS, OUT_COORDS
                    )

                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_polsum(
                        omega, cos_theta, 1.0
                    )

                    PSP = PhaseSpacePoint(PROC, MODEL, PS_DEF, IN_PS, OUT_PS)
                    test_val = unsafe_differential_cross_section(PSP)

                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "x-pol and spin summed" begin
                PROC = Compton(PolX())

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = (omega,)
                    OUT_COORDS = (cos_theta, phi)
                    IN_PS, OUT_PS = QEDprocesses._generate_momenta(
                        PROC, MODEL, PS_DEF, IN_COORDS, OUT_COORDS
                    )

                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_xpol(
                        omega, cos_theta, phi, 1.0
                    )

                    PSP = PhaseSpacePoint(PROC, MODEL, PS_DEF, IN_PS, OUT_PS)
                    test_val = unsafe_differential_cross_section(PSP)

                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "y-pol and spin summed" begin
                PROC = Compton(PolY())

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = (omega,)
                    OUT_COORDS = (cos_theta, phi)
                    IN_PS, OUT_PS = QEDprocesses._generate_momenta(
                        PROC, MODEL, PS_DEF, IN_COORDS, OUT_COORDS
                    )

                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_ypol(
                        omega, cos_theta, phi, 1.0
                    )

                    PSP = PhaseSpacePoint(PROC, MODEL, PS_DEF, IN_PS, OUT_PS)
                    test_val = unsafe_differential_cross_section(PSP)

                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end
        end
        @testset "total cross section" begin
            @testset "spin and pol summed" begin
                PROC = Compton()
                # Klein-Nishina: total cross section
                function klein_nishina_total_cross_section(in_ps)
                    function func(x)
                        return unsafe_differential_cross_section(
                            PhaseSpacePoint(
                                Compton(), PerturbativeQED(), PS_DEF, in_ps, (x, 0.0)
                            ),
                        )
                    end
                    res, err = quadgk(func, -1, 1)

                    # note: mul by 2pi instead of the phi-integration
                    return 2 * pi * res
                end

                IN_COORDS = (omega,)
                groundtruth = klein_nishina_total_cross_section(IN_COORDS)
                test_val = @inferred QEDprocesses.total_cross_section(
                    InPhaseSpacePoint(PROC, MODEL, PS_DEF, IN_COORDS)
                )
                @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    OUT_COORDS = (cos_theta, phi)

                    test_val = @inferred QEDprocesses.total_cross_section(
                        PhaseSpacePoint(PROC, MODEL, PS_DEF, IN_COORDS, OUT_COORDS)
                    )
                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)

                    out_moms = momenta(
                        PhaseSpacePoint(PROC, MODEL, PS_DEF, IN_COORDS, OUT_COORDS),
                        Outgoing(),
                    )
                    @test_throws MethodError QEDprocesses.total_cross_section(
                        OutPhaseSpacePoint(PROC, MODEL, PS_DEF, out_moms)
                    )
                end
            end
        end
    end
end
