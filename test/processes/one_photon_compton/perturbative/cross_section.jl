using QEDbase
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
const IN_PSL = ComptonRestSystem(Energy(2))
const OUT_PSL = ComptonSphericalLayout(IN_PSL)

const OMEGAS = (1e-6 * rand(RNG), 1e-3 * rand(RNG), rand(RNG), 1e3 * rand(RNG))
const COS_THETAS = [-1.0, 2 * rand(RNG) - 1, 0.0, 1.0]
const PHIS = [0, 2 * pi, rand(RNG) * 2 * pi]

@testset "perturbative" begin
    @testset "$omega" for omega in OMEGAS
        @testset "differential cross section" begin
            @testset "spin and pol summed" begin
                PROC = Compton()

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = (omega,)
                    OUT_COORDS = (cos_theta, phi)

                    PSP = PhaseSpacePoint(PROC, MODEL, OUT_PSL, IN_COORDS, OUT_COORDS)
                    test_val = unsafe_differential_cross_section(PSP)

                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_polsum_elab_sph(
                        omega, cos_theta, 1.0
                    )

                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "x-pol and spin summed" begin
                PROC = Compton(PolX())

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = (omega,)
                    OUT_COORDS = (cos_theta, phi)
                    PSP = PhaseSpacePoint(PROC, MODEL, OUT_PSL, IN_COORDS, OUT_COORDS)
                    test_val = unsafe_differential_cross_section(PSP)

                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_xpol_elab_sph(
                        omega, cos_theta, phi, 1.0
                    )

                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "y-pol and spin summed" begin
                PROC = Compton(PolY())

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = (omega,)
                    OUT_COORDS = (cos_theta, phi)
                    PSP = PhaseSpacePoint(PROC, MODEL, OUT_PSL, IN_COORDS, OUT_COORDS)
                    test_val = unsafe_differential_cross_section(PSP)

                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_ypol_elab_sph(
                        omega, cos_theta, phi, 1.0
                    )

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
                            PhaseSpacePoint(Compton(), MODEL, OUT_PSL, in_ps, (x, 0.0))
                        )
                    end
                    res, err = quadgk(func, -1, 1)

                    # note: mul by 2pi instead of the phi-integration
                    return 2 * pi * res
                end

                IN_COORDS = (omega,)
                groundtruth = klein_nishina_total_cross_section(IN_COORDS)
                test_val = @inferred QEDprocesses.total_cross_section(
                    InPhaseSpacePoint(PROC, MODEL, IN_PSL, IN_COORDS)
                )

                @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    OUT_COORDS = (cos_theta, phi)

                    test_val = @inferred QEDprocesses.total_cross_section(
                        PhaseSpacePoint(PROC, MODEL, OUT_PSL, IN_COORDS, OUT_COORDS)
                    )
                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)

                    out_moms = momenta(
                        PhaseSpacePoint(PROC, MODEL, OUT_PSL, IN_COORDS, OUT_COORDS),
                        Outgoing(),
                    )
                    @test_throws MethodError QEDprocesses.total_cross_section(
                        OutPhaseSpacePoint(PROC, MODEL, OUT_PSL, out_moms)
                    )
                end
            end
        end
    end
end
