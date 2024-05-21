using Random
using QEDbase
using QEDprocesses

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

include("test_implementation/TestImplementation.jl")
TESTMODEL = TestImplementation.TestModel()
TESTPSDEF = TestImplementation.TestPhasespaceDef()

@testset "($N_INCOMING,$N_OUTGOING)" for (N_INCOMING, N_OUTGOING) in Iterators.product(
    (1, rand(RNG, 2:8)), (1, rand(RNG, 2:8))
)
    INCOMING_PARTICLES = rand(RNG, TestImplementation.PARTICLE_SET, N_INCOMING)
    OUTGOING_PARTICLES = rand(RNG, TestImplementation.PARTICLE_SET, N_OUTGOING)

    TESTPROC = TestImplementation.TestProcess(INCOMING_PARTICLES, OUTGOING_PARTICLES)

    # single ps points
    p_in_phys = TestImplementation._rand_momenta(RNG, N_INCOMING)
    p_in_phys_invalid = TestImplementation._rand_momenta(RNG, N_INCOMING + 1)
    p_in_unphys = TestImplementation._rand_in_momenta_failing(RNG, N_INCOMING)
    p_in_unphys_invalid = TestImplementation._rand_in_momenta_failing(RNG, N_INCOMING + 1)

    p_out_phys = TestImplementation._rand_momenta(RNG, N_OUTGOING)
    p_out_phys_invalid = TestImplementation._rand_momenta(RNG, N_OUTGOING + 1)
    p_out_unphys = TestImplementation._rand_out_momenta_failing(RNG, N_OUTGOING)
    p_out_unphys_invalid = TestImplementation._rand_out_momenta_failing(RNG, N_OUTGOING + 1)

    # sets of ps points
    p_in_set_phys = TestImplementation._rand_momenta(RNG, N_INCOMING, 2)
    p_in_set_unphys_mix = TestImplementation._rand_in_momenta_failing_mix(
        RNG, N_INCOMING, 2
    )
    p_in_set_unphys_all = TestImplementation._rand_in_momenta_failing_all(
        RNG, N_INCOMING, 2
    )
    p_in_set_phys_invalid = TestImplementation._rand_momenta(RNG, N_INCOMING + 1, 2)
    p_in_set_unphys_mix_invalid = TestImplementation._rand_in_momenta_failing_mix(
        RNG, N_INCOMING + 1, 2
    )
    p_in_set_unphys_all_invalid = TestImplementation._rand_in_momenta_failing_all(
        RNG, N_INCOMING + 1, 2
    )

    p_out_set_phys = TestImplementation._rand_momenta(RNG, N_OUTGOING, 2)
    p_out_set_unphys_mix = TestImplementation._rand_out_momenta_failing_mix(
        RNG, N_OUTGOING, 2
    )
    p_out_set_unphys_all = TestImplementation._rand_out_momenta_failing_all(
        RNG, N_OUTGOING, 2
    )
    p_out_set_phys_invalid = TestImplementation._rand_momenta(RNG, N_OUTGOING + 1, 2)
    p_out_set_unphys_mix_invalid = TestImplementation._rand_out_momenta_failing_mix(
        RNG, N_OUTGOING + 1, 2
    )
    p_out_set_unphys_all_invalid = TestImplementation._rand_out_momenta_failing_all(
        RNG, N_OUTGOING + 1, 2
    )

    p_in_all = (
        p_in_phys,
        p_in_unphys,
        p_in_phys_invalid,
        p_in_unphys_invalid,
        p_in_set_phys,
        p_in_set_unphys_mix,
        p_in_set_unphys_all,
        p_in_set_phys_invalid,
        p_in_set_unphys_mix_invalid,
        p_in_set_unphys_all_invalid,
    )

    p_out_all = (
        p_out_phys,
        p_out_phys_invalid,
        p_out_unphys,
        p_out_unphys_invalid,
        p_out_set_phys,
        p_out_set_unphys_mix,
        p_out_set_unphys_all,
        p_out_set_phys_invalid,
        p_out_set_unphys_mix_invalid,
        p_out_set_unphys_all_invalid,
    )
    # all combinations
    p_combs = Iterators.product(p_in_all, p_out_all)

    p_in_all_valid = (
        p_in_phys, p_in_unphys, p_in_set_phys, p_in_set_unphys_mix, p_in_set_unphys_all
    )

    p_out_all_valid = (
        p_out_phys, p_out_unphys, p_out_set_phys, p_out_set_unphys_mix, p_out_set_unphys_all
    )

    # all valid combinations
    p_combs_valid = Iterators.product(p_in_all_valid, p_out_all_valid)

    p_in_all_phys = (p_in_phys, p_in_set_phys)
    p_out_all_phys = (p_out_phys, p_out_set_phys)

    p_combs_phys = Iterators.product(p_in_all_phys, p_out_all_phys)

    @testset "cross section" begin
        @testset "unsafe" begin
            @testset "compute" begin
                for (P_IN, P_OUT) in p_combs_phys
                    diffCS = unsafe_differential_cross_section(
                        TESTPROC, TESTMODEL, TESTPSDEF, P_IN, P_OUT
                    )
                    groundtruth = TestImplementation._groundtruth_unsafe_diffCS(
                        TESTPROC, P_IN, P_OUT
                    )
                    @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "compute on phase space points" begin
                PS_POINT = generate_phase_space(
                    TESTPROC, TESTMODEL, TESTPSDEF, p_in_phys, p_out_phys
                )
                diffCS_on_psp = unsafe_differential_cross_section(PS_POINT)
                groundtruth = TestImplementation._groundtruth_unsafe_diffCS(
                    TESTPROC, p_in_phys, p_out_phys
                )
                @test isapprox(diffCS_on_psp, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "invalid input" begin
                for (P_IN, P_OUT) in p_combs

                    # filter out all valid combinations
                    if !((P_IN, P_OUT) in p_combs_valid)
                        @test_throws DimensionMismatch unsafe_differential_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN, P_OUT
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        COORDS_OUT = TestImplementation.flat_components(P_OUT)
                        @test_throws DimensionMismatch unsafe_differential_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, COORDS_OUT
                        )
                    end
                end
            end

            @testset "safe" begin
                @testset "compute" begin
                    for (P_IN, P_OUT) in p_combs_valid
                        diffCS_on_moms = differential_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN, P_OUT
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        COORDS_OUT = TestImplementation.flat_components(P_OUT)
                        diffCS_on_coords = differential_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, COORDS_OUT
                        )
                        groundtruth = TestImplementation._groundtruth_safe_diffCS(
                            TESTPROC, P_IN, P_OUT
                        )
                        @test isapprox(diffCS_on_moms, groundtruth, atol=ATOL, rtol=RTOL)
                        @test isapprox(diffCS_on_coords, groundtruth, atol=ATOL, rtol=RTOL)
                    end
                end

                @testset "compute on phase space points" begin
                    PS_POINT = generate_phase_space(
                        TESTPROC, TESTMODEL, TESTPSDEF, p_in_phys, p_out_phys
                    )
                    diffCS_on_psp = differential_cross_section(PS_POINT)
                    groundtruth = TestImplementation._groundtruth_safe_diffCS(
                        TESTPROC, p_in_phys, p_out_phys
                    )
                    @test isapprox(diffCS_on_psp, groundtruth, atol=ATOL, rtol=RTOL)
                end

                @testset "invalid input" begin
                    for (P_IN, P_OUT) in p_combs

                        # filter out all valid combinations
                        if !((P_IN, P_OUT) in p_combs_valid)
                            @test_throws DimensionMismatch differential_cross_section(
                                TESTPROC, TESTMODEL, TESTPSDEF, P_IN, P_OUT
                            )

                            COORDS_IN = TestImplementation.flat_components(P_IN)
                            COORDS_OUT = TestImplementation.flat_components(P_OUT)
                            @test_throws DimensionMismatch differential_cross_section(
                                TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, COORDS_OUT
                            )
                        end
                    end
                end
            end

            @testset "total cross section" begin
                @testset "compute" begin
                    for P_IN in (p_in_phys, p_in_set_phys)
                        groundtruth = TestImplementation._groundtruth_total_cross_section(
                            P_IN
                        )
                        totCS_on_moms = total_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        totCS_on_coords = total_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN
                        )

                        @test isapprox(totCS_on_moms, groundtruth, atol=ATOL, rtol=RTOL)
                        @test isapprox(totCS_on_coords, groundtruth, atol=ATOL, rtol=RTOL)
                    end
                end
                @testset "invalid input" begin
                    for P_IN in (p_in_phys_invalid, p_in_set_phys_invalid)
                        @test_throws DimensionMismatch total_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        @test_throws DimensionMismatch total_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN
                        )
                    end
                end
            end
        end

        @testset "differential probability" begin
            @testset "unsafe" begin
                @testset "compute" begin
                    for (P_IN, P_OUT) in p_combs_phys
                        prob_on_moms = unsafe_differential_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN, P_OUT
                        )
                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        COORDS_OUT = TestImplementation.flat_components(P_OUT)
                        prob_on_coords = unsafe_differential_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, COORDS_OUT
                        )
                        groundtruth = TestImplementation._groundtruth_unsafe_probability(
                            TESTPROC, P_IN, P_OUT
                        )
                        @test isapprox(prob_on_moms, groundtruth, atol=ATOL, rtol=RTOL)
                        @test isapprox(prob_on_coords, groundtruth, atol=ATOL, rtol=RTOL)
                    end
                end

                @testset "compute on phase space points" begin
                    PS_POINT = generate_phase_space(
                        TESTPROC, TESTMODEL, TESTPSDEF, p_in_phys, p_out_phys
                    )
                    prop_on_psp = unsafe_differential_probability(PS_POINT)
                    groundtruth = TestImplementation._groundtruth_unsafe_probability(
                        TESTPROC, p_in_phys, p_out_phys
                    )
                    @test isapprox(prop_on_psp, groundtruth, atol=ATOL, rtol=RTOL)
                end

                @testset "invalid input" begin
                    for (P_IN, P_OUT) in p_combs

                        # filter out all valid combinations
                        if !((P_IN, P_OUT) in p_combs_valid)
                            @test_throws DimensionMismatch unsafe_differential_probability(
                                TESTPROC, TESTMODEL, TESTPSDEF, P_IN, P_OUT
                            )

                            COORDS_IN = TestImplementation.flat_components(P_IN)
                            COORDS_OUT = TestImplementation.flat_components(P_OUT)
                            @test_throws DimensionMismatch unsafe_differential_probability(
                                TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, COORDS_OUT
                            )
                        end
                    end
                end
            end
            @testset "safe" begin
                @testset "compute" begin
                    for (P_IN, P_OUT) in p_combs_valid
                        prob_on_moms = differential_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN, P_OUT
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        COORDS_OUT = TestImplementation.flat_components(P_OUT)
                        prob_on_coords = differential_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, COORDS_OUT
                        )
                        groundtruth = TestImplementation._groundtruth_safe_probability(
                            TESTPROC, P_IN, P_OUT
                        )
                        @test isapprox(prob_on_moms, groundtruth, atol=ATOL, rtol=RTOL)
                        @test isapprox(prob_on_coords, groundtruth, atol=ATOL, rtol=RTOL)
                    end
                end

                @testset "compute on phase space points" begin
                    PS_POINT = generate_phase_space(
                        TESTPROC, TESTMODEL, TESTPSDEF, p_in_phys, p_out_phys
                    )
                    prop_on_psp = differential_probability(PS_POINT)
                    groundtruth = TestImplementation._groundtruth_safe_probability(
                        TESTPROC, p_in_phys, p_out_phys
                    )
                    @test isapprox(prop_on_psp, groundtruth, atol=ATOL, rtol=RTOL)
                end

                @testset "invalid input" begin
                    for (P_IN, P_OUT) in p_combs

                        # filter out all valid combinations
                        if !((P_IN, P_OUT) in p_combs_valid)
                            @test_throws DimensionMismatch differential_probability(
                                TESTPROC, TESTMODEL, TESTPSDEF, P_IN, P_OUT
                            )

                            COORDS_IN = TestImplementation.flat_components(P_IN)
                            COORDS_OUT = TestImplementation.flat_components(P_OUT)
                            @test_throws DimensionMismatch differential_probability(
                                TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, COORDS_OUT
                            )
                        end
                    end
                end
            end

            @testset "total probability" begin
                @testset "compute" begin
                    for P_IN in (p_in_phys, p_in_set_phys)
                        groundtruth = TestImplementation._groundtruth_total_probability(
                            P_IN
                        )
                        totCS_on_moms = total_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        totCS_on_coords = total_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN
                        )

                        @test isapprox(totCS_on_moms, groundtruth, atol=ATOL, rtol=RTOL)
                        @test isapprox(totCS_on_coords, groundtruth, atol=ATOL, rtol=RTOL)
                    end
                end
                @testset "invalid input" begin
                    for P_IN in (p_in_phys_invalid, p_in_set_phys_invalid)
                        @test_throws DimensionMismatch total_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        @test_throws DimensionMismatch total_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN
                        )
                    end
                end
            end
        end
    end
end
