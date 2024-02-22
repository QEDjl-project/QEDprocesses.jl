
using QEDbase
using QEDprocesses
using Random
using StaticArrays
using QuadGK

RNG = MersenneTwister(77697185)
ATOL = eps()
RTOL = sqrt(eps())

include("groundtruths.jl")

MODEL = PerturbativeQED()
PS_DEF = PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame())
OMEGAS = (1e-6 * rand(RNG), 1e-3 * rand(RNG), rand(RNG), 1e3 * rand(RNG))
#OMEGAS = (rand(RNG),)

COS_THETAS = [-1.0, 2 * rand(RNG) - 1, 0.0, 1.0]
PHIS = [0, 2 * pi, rand(RNG) * 2 * pi]

@testset "perturbative kinematics" begin
    PROC = Compton()
    @testset "momentum generation" begin
        @testset "$om, $cth, $phi" for (om, cth, phi) in
                                       Iterators.product(OMEGAS, COS_THETAS, PHIS)
            IN_COORDS = [om]
            OUT_COORDS = [cth, phi]
            IN_PS, OUT_PS = QEDprocesses._generate_momenta(
                PROC, MODEL, PS_DEF, IN_COORDS, PS_DEF, OUT_COORDS
            )
            in_mom_square = getMass2.(IN_PS)
            out_mom_square = getMass2.(OUT_PS)
            in_masses = mass.(incoming_particles(PROC)) .^ 2
            out_masses = mass.(outgoing_particles(PROC))
            @test isapprox(in_mom_square, SVector(in_masses))
            @test isapprox(out_mom_square, SVector(out_masses))
        end
    end
end

@testset "perturbative cross section" begin
    @testset "$omega" for omega in OMEGAS
        @testset "differential cross section" begin
            @testset "spin and pol summed" begin
                PROC = Compton()

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = [omega]
                    OUT_COORDS = [cos_theta, phi]
                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_polsum(
                        omega, cos_theta, 1.0
                    )
                    test_val = unsafe_differential_cross_section(
                        PROC, MODEL, PS_DEF, IN_COORDS, PS_DEF, OUT_COORDS
                    )
                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "x-pol and spin summed" begin
                PROC = Compton(PolX())

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = [omega]
                    OUT_COORDS = [cos_theta, phi]
                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_xpol(
                        omega, cos_theta, phi, 1.0
                    )
                    test_val = unsafe_differential_cross_section(
                        PROC, MODEL, PS_DEF, IN_COORDS, PS_DEF, OUT_COORDS
                    )
                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "y-pol and spin summed" begin
                PROC = Compton(PolY())

                @testset "$cos_theta $phi" for (cos_theta, phi) in
                                               Iterators.product(COS_THETAS, PHIS)
                    IN_COORDS = [omega]
                    OUT_COORDS = [cos_theta, phi]
                    groundtruth = _groundtruth_pert_compton_diffCS_spinsum_ypol(
                        omega, cos_theta, phi, 1.0
                    )
                    test_val = unsafe_differential_cross_section(
                        PROC, MODEL, PS_DEF, IN_COORDS, PS_DEF, OUT_COORDS
                    )
                    @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end
        end
    end
    # @testset "total cross section" begin
    #     @testset "spin and pol summed" begin
    #         PROC = Compton()
    #         # Klein-Nishina: total cross section
    #         function klein_nishina_total_cross_section(om, mass)
    #             function func(x)
    #                 return differential_cross_section_on_coord(
    #                     Compton(), PerturbativeQED(), om, [x, 0.0]
    #                 )
    #             end
    #             res, err = quadgk(func, -1, 1)
    #
    #             # note: mul by 2pi instead of the phi-integration
    #             return 2 * pi * res
    #         end
    #
    #         groundtruth = klein_nishina_total_cross_section(OMEGA, MASS)
    #         test_val = @inferred total_cross_section_on_coord(PROC, MODEL, OMEGA)
    #         @test isapprox(test_val, groundtruth, atol=ATOL, rtol=RTOL)
    #     end
    # end
end
