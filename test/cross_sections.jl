using Random
using Suppressor
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

    # single ps points
    p_in_phys = _rand_momenta(RNG, N_INCOMING)
    p_in_phys_invalid = _rand_momenta(RNG, N_INCOMING + 1)
    p_in_unphys = _rand_momenta_failed(RNG, N_INCOMING)
    p_in_unphys_invalid = _rand_momenta_failed(RNG, N_INCOMING + 1)
    p_out_phys = _rand_momenta(RNG, N_OUTGOING)
    p_out_phys_invaid = _rand_momenta(RNG, N_OUTGOING + 1)

    # sets of ps points
    p_in_set_phys = _rand_momenta(RNG, N_INCOMING, 2)
    p_in_set_unphys_mix = _rand_momenta_failed_mix(RNG, N_INCOMING, 2)
    p_in_set_unphys_all = _rand_momenta_failed_all(RNG, N_INCOMING, 2)
    p_in_set_phys_invalid = _rand_momenta(RNG, N_INCOMING + 1, 2)
    p_in_set_unphys_mix_invalid = _rand_momenta_failed_mix(RNG, N_INCOMING + 1, 2)
    p_in_set_unphys_all_invalid = _rand_momenta_failed_all(RNG, N_INCOMING + 1, 2)
    p_out_set_phys = _rand_momenta(RNG, N_OUTGOING, 2)
    p_out_set_phys_invalid = _rand_momenta(RNG, N_OUTGOING + 1, 2)

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

    p_out_all = (p_out_phys, p_out_phys_invaid, p_out_set_phys, p_out_set_phys_invalid)
    # all combinations
    p_combs = Iterators.product(p_in_all, p_out_all)

    p_in_all_valid = (
        p_in_phys, p_in_unphys, p_in_set_phys, p_in_set_unphys_mix, p_in_set_unphys_all
    )

    p_out_all_valid = (p_out_phys, p_out_set_phys)

    # all valid combinations
    p_combs_valid = collect(Iterators.product(p_in_all_valid, p_out_all_valid))

    p_in_all_phys = (p_in_phys, p_in_set_phys)
    p_out_all_phys = (p_out_phys, p_out_set_phys)

    p_combs_phys = Iterators.product(p_in_all_phys, p_out_all_phys)

    @testset "cross section" begin
        @testset "unsafe" begin
            @testset "compute" begin
                for (P_IN, P_OUT) in p_combs_phys
                    diffCS = unsafe_differential_cross_section(
                        TESTPROC, TESTMODEL, TESTPSDEF, P_IN, TESTPSDEF, P_OUT
                    )
                    groundtruth = TestImplementation._groundtruth_unsafe_diffCS(
                        TESTPROC, P_IN, P_OUT
                    )
                    @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "invalid input" begin
                for (P_IN, P_OUT) in p_combs

                    # filter out all valid combinations
                    if !((P_IN, P_OUT) in p_combs_valid)
                        @test_throws DimensionMismatch unsafe_differential_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN, TESTPSDEF, P_OUT
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        COORDS_OUT = TestImplementation.flat_components(P_OUT)
                        @test_throws DimensionMismatch unsafe_differential_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, TESTPSDEF, COORDS_OUT
                        )
                    end
                end
            end

            @testset "safe" begin
                @testset "compute" begin
                    for (P_IN, P_OUT) in p_combs_valid
                        diffCS_on_moms = differential_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN, TESTPSDEF, P_OUT
                        )
                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        COORDS_OUT = TestImplementation.flat_components(P_OUT)
                        diffCS_on_coords = differential_cross_section(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, TESTPSDEF, COORDS_OUT
                        )
                        groundtruth = TestImplementation._groundtruth_safe_diffCS(
                            TESTPROC, P_IN, P_OUT
                        )
                        @test isapprox(diffCS_on_moms, groundtruth, atol=ATOL, rtol=RTOL)
                        @test isapprox(diffCS_on_coords, groundtruth, atol=ATOL, rtol=RTOL)
                    end
                end

                @testset "invalid input" begin
                    for (P_IN, P_OUT) in p_combs

                        # filter out all valid combinations
                        if !((P_IN, P_OUT) in p_combs_valid)
                            @test_throws DimensionMismatch differential_cross_section(
                                TESTPROC, TESTMODEL, TESTPSDEF, P_IN, TESTPSDEF, P_OUT
                            )

                            COORDS_IN = TestImplementation.flat_components(P_IN)
                            COORDS_OUT = TestImplementation.flat_components(P_OUT)
                            @test_throws DimensionMismatch differential_cross_section(
                                TESTPROC,
                                TESTMODEL,
                                TESTPSDEF,
                                COORDS_IN,
                                TESTPSDEF,
                                COORDS_OUT,
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
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN, TESTPSDEF, P_OUT
                        )
                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        COORDS_OUT = TestImplementation.flat_components(P_OUT)
                        prob_on_coords = unsafe_differential_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, TESTPSDEF, COORDS_OUT
                        )
                        groundtruth = TestImplementation._groundtruth_unsafe_probability(
                            TESTPROC, P_IN, P_OUT
                        )
                        @test isapprox(prob_on_moms, groundtruth, atol=ATOL, rtol=RTOL)
                        @test isapprox(prob_on_coords, groundtruth, atol=ATOL, rtol=RTOL)
                    end
                end

                @testset "invalid input" begin
                    for (P_IN, P_OUT) in p_combs

                        # filter out all valid combinations
                        if !((P_IN, P_OUT) in p_combs_valid)
                            @test_throws DimensionMismatch unsafe_differential_probability(
                                TESTPROC, TESTMODEL, TESTPSDEF, P_IN, TESTPSDEF, P_OUT
                            )

                            COORDS_IN = TestImplementation.flat_components(P_IN)
                            COORDS_OUT = TestImplementation.flat_components(P_OUT)
                            @test_throws DimensionMismatch unsafe_differential_probability(
                                TESTPROC,
                                TESTMODEL,
                                TESTPSDEF,
                                COORDS_IN,
                                TESTPSDEF,
                                COORDS_OUT,
                            )
                        end
                    end
                end
            end
            @testset "safe" begin
                @testset "compute" begin
                    for (P_IN, P_OUT) in p_combs_valid
                        prob_on_moms = differential_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, P_IN, TESTPSDEF, P_OUT
                        )

                        COORDS_IN = TestImplementation.flat_components(P_IN)
                        COORDS_OUT = TestImplementation.flat_components(P_OUT)
                        prob_on_coords = differential_probability(
                            TESTPROC, TESTMODEL, TESTPSDEF, COORDS_IN, TESTPSDEF, COORDS_OUT
                        )
                        groundtruth = TestImplementation._groundtruth_safe_probability(
                            TESTPROC, P_IN, P_OUT
                        )
                        @test isapprox(prob_on_moms, groundtruth, atol=ATOL, rtol=RTOL)
                        @test isapprox(prob_on_coords, groundtruth, atol=ATOL, rtol=RTOL)
                    end
                end

                @testset "invalid input" begin
                    for (P_IN, P_OUT) in p_combs

                        # filter out all valid combinations
                        if !((P_IN, P_OUT) in p_combs_valid)
                            @test_throws DimensionMismatch differential_probability(
                                TESTPROC, TESTMODEL, TESTPSDEF, P_IN, TESTPSDEF, P_OUT
                            )

                            COORDS_IN = TestImplementation.flat_components(P_IN)
                            COORDS_OUT = TestImplementation.flat_components(P_OUT)
                            @test_throws DimensionMismatch differential_probability(
                                TESTPROC,
                                TESTMODEL,
                                TESTPSDEF,
                                COORDS_IN,
                                TESTPSDEF,
                                COORDS_OUT,
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
                
         @testset "safe differential cross section" begin
            @testset "compute vector-vector" begin
                diffCS_physical = differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    p_out,
                )
                groundtruth_physical = _groundtruth_unsafe_diffCS(TestProcess(), p_in, p_out)
                @test isapprox(diffCS_physical, groundtruth_physical, atol=ATOL, rtol=RTOL)

                diffCS_unphysical= differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_unphysical,
                    TestPhasespaceDef(),
                    p_out,
                )
                groundtruth_unphysical = zero(Float64) 
                @test isapprox(diffCS_unphysical, groundtruth_unphysical, atol=ATOL, rtol=RTOL)

            end

            @testset "compute vector-matrix" begin
                diffCS_physical = differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    p_out_set,
                )

                groundtruth_physical = Vector{QEDprocesses._base_component_type(p_in)}(
                    undef, size(p_out_set, 2)
                )
                for i in 1:size(p_out_set, 2)
                    groundtruth_physical[i] = _groundtruth_unsafe_diffCS(
                        TestProcess(), p_in, view(p_out_set, :, i)
                    )
                end
                @test isapprox(diffCS_physical, groundtruth_physical, atol=ATOL, rtol=RTOL)

                diffCS_unphysical = differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_unphysical,
                    TestPhasespaceDef(),
                    p_out_set,
                )

                groundtruth_unphysical = zeros(size(p_out_set,2))
                @test isapprox(diffCS_unphysical, groundtruth_unphysical, atol=ATOL, rtol=RTOL)

            end

            @testset "compute matrix-vector" begin
                diffCS_physical = differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    p_out,
                )
                groundtruth_physical = Vector{QEDprocesses._base_component_type(p_in_set)}(
                    undef, size(p_in_set, 2)
                )
                for i in 1:size(p_in_set, 2)
                    groundtruth_physical[i] = _groundtruth_unsafe_diffCS(
                        TestProcess(), view(p_in_set, :, i), p_out
                    )
                end
                @test isapprox(diffCS_physical, groundtruth_physical, atol=ATOL, rtol=RTOL)


                
                diffCS_unphys_all= differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set_unphys_all,
                    TestPhasespaceDef(),
                    p_out,
                )
                groundtruth_unphys_all = zeros(size(p_in_set,2))
                @test isapprox(diffCS_unphys_all, groundtruth_unphys_all, atol=ATOL, rtol=RTOL)
                    
                diffCS_unphys_mix= differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set_unphys_mix,
                    TestPhasespaceDef(),
                    p_out,
                )
                groundtruth_unphys_mix= Vector{QEDprocesses._base_component_type(p_in_set)}(
                    undef, size(p_in_set, 2)
                )
                for i in 1:size(p_in_set, 2)
                        if i==1
                            groundtruth_unphys_mix[i] = zero(Float64)
                        else
                            groundtruth_unphys_mix[i] = _groundtruth_unsafe_diffCS(
                                TestProcess(), view(p_in_set_unphys_mix, :, i), p_out
                            )
                        end
                end
                @test isapprox(diffCS_unphys_mix, groundtruth_unphys_mix, atol=ATOL, rtol=RTOL)
            end

            @testset "compute matrix-matrix" begin
                diffCS_physical= differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    p_out_set,
                )
                groundtruth_physical= Matrix{QEDprocesses._base_component_type(p_in_set)}(
                    undef, size(p_in_set, 2), size(p_out_set, 2)
                )
                for i in 1:size(p_in_set, 2)
                    for j in 1:size(p_out_set, 2)
                        groundtruth_physical[i, j] = _groundtruth_unsafe_diffCS(
                            TestProcess(), view(p_in_set, :, i), view(p_out_set, :, j)
                        )
                    end
                end
                @test isapprox(diffCS_physical, groundtruth_physical, atol=ATOL, rtol=RTOL)
                    
                diffCS_unphys_all= differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set_unphys_all,
                    TestPhasespaceDef(),
                    p_out_set,
                )
                    groundtruth_unphys_all= zeros(size(p_in_set,2),size(p_out_set,2))
                @test isapprox(diffCS_unphys_all, groundtruth_unphys_all, atol=ATOL, rtol=RTOL)

                diffCS_unphys_mix = differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set_unphys_mix,
                    TestPhasespaceDef(),
                    p_out_set,
                )
                groundtruth_unphys_mix= Matrix{QEDprocesses._base_component_type(p_in_set)}(
                    undef, size(p_in_set, 2), size(p_out_set, 2)
                )
                for i in 1:size(p_in_set, 2)
                        for j in 1:size(p_out_set, 2)
                            if i==1
                                groundtruth_unphys_mix[i,j] = zero(Float64)
                            else
                                groundtruth_unphys_mix[i, j] = _groundtruth_unsafe_diffCS(
                                    TestProcess(), view(p_in_set, :, i), view(p_out_set, :, j)
                                )
                            end
                        end
                end
                @test isapprox(diffCS_unphys_mix, groundtruth_unphys_mix, atol=ATOL, rtol=RTOL)

            end

            @testset "fail vector-vector" begin
                    p_in_physical_invalid= _rand_momenta(RNG, N_INCOMING + 1)
                    p_in_unphysical_invalid = deepcopy(p_in_physical_invalid)
                    p_in_unphysical_invalid[1] = SFourMomentum(zeros(4))
                    for (P_IN,P_OUT) in Iterators.product(
                    (p_in,p_in_unphysical,p_in_physical_invalid,p_in_unphysical_invalid),
                    (p_out, _rand_momenta(RNG, N_OUTGOING + 1))
                    )
                        if (P_IN != p_in  || P_OUT != p_out) && (P_IN != p_in_unphysical || P_OUT != p_out)
                            @test_throws DimensionMismatch differential_cross_section(
                                TestProcess(),
                                TestModel(),
                                TestPhasespaceDef(),
                                P_IN,
                                TestPhasespaceDef(),
                                P_OUT,
                            )
                        end
                    end
            end

            @testset "fail vector-matrix" begin
                    p_in_physical_invalid= _rand_momenta(RNG, N_INCOMING + 1)
                    p_in_unphysical_invalid = deepcopy(p_in_physical_invalid)
                    p_in_unphysical_invalid[1] = SFourMomentum(zeros(4))
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_physical_invalid,
                    TestPhasespaceDef(),
                    p_out_set,
                )
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_unphysical_invalid,
                    TestPhasespaceDef(),
                    p_out_set,
                )
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )

                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_unphysical,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )
            end

            @testset "fail matrix-vector" begin
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                    TestPhasespaceDef(),
                    p_out,
                )
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1),
                )
            end

            @testset "fail matrix-matrix" begin
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                    TestPhasespaceDef(),
                    p_out_set,
                )
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )
            end
        end
        @testset "safe" begin
            @testset "compute" begin
                for (P_IN, P_OUT) in p_combs_valid
                    diffCS = differential_cross_section(
                        TestProcess(),
                        TestModel(),
                        TestPhasespaceDef(),
                        P_IN,
                        TestPhasespaceDef(),
                        P_OUT,
                    )
                    groundtruth = _groundtruth_safe_diffCS(TestProcess(), P_IN, P_OUT)
                    @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "invalid input" begin
                for (P_IN, P_OUT) in p_combs

                    # filter out all valid combinations
                    if !((P_IN, P_OUT) in p_combs_valid)
                        @test_throws DimensionMismatch differential_cross_section(
                            TestProcess(),
                            TestModel(),
                            TestPhasespaceDef(),
                            P_IN,
                            TestPhasespaceDef(),
                            P_OUT,
                        )
                    end
                end
            end
        end
    end

    @testset "probability" begin
        @testset "unsafe compute" begin
            for (P_IN, P_OUT) in p_combs_phys
                prob = unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    P_IN,
                    TestPhasespaceDef(),
                    P_OUT,
                )
                groundtruth = _groundtruth_unsafe_probability(TestProcess(), P_IN, P_OUT)
                @test isapprox(prob, groundtruth, atol=ATOL, rtol=RTOL)
            end
        end

        @testset "unsafe invalid input" begin
            for (P_IN, P_OUT) in p_combs

                # filter out all valid combinations
                if !((P_IN, P_OUT) in p_combs_valid)
                    @test_throws DimensionMismatch unsafe_probability(
                        TestProcess(),
                        TestModel(),
                        TestPhasespaceDef(),
                        P_IN,
                        TestPhasespaceDef(),
                        P_OUT,
                    )
                end
            end
        end
        @testset "safe compute" begin
            for (P_IN, P_OUT) in p_combs_valid
                prob = probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    P_IN,
                    TestPhasespaceDef(),
                    P_OUT,
                )
                groundtruth = _groundtruth_safe_probability(TestProcess(), P_IN, P_OUT)
                @test isapprox(prob, groundtruth, atol=ATOL, rtol=RTOL)
            end
        end

        @testset "safe invalid input" begin
            for (P_IN, P_OUT) in p_combs

                # filter out all valid combinations
                if !((P_IN, P_OUT) in p_combs_valid)
                    @test_throws DimensionMismatch probability(
                        TestProcess(),
                        TestModel(),
                        TestPhasespaceDef(),
                        P_IN,
                        TestPhasespaceDef(),
                        P_OUT,
                    )
                end
            end
        end
    end
end
