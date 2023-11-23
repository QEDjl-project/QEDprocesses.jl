using Random
using QEDbase
using QEDprocesses

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

include("utils/utils.jl")

@testset "($N_INCOMING,$N_OUTGOING)" for (N_INCOMING, N_OUTGOING) in Iterators.product(
(1, rand(RNG, 2:8)), (1, rand(RNG, 2:8))
)
    INCOMING_PARTICLES = rand(RNG, PARTICLE_SET, N_INCOMING)
    OUTGOING_PARTICLES = rand(RNG, PARTICLE_SET, N_OUTGOING)

    QEDprocesses.incoming_particles(::TestProcess) = INCOMING_PARTICLES
    QEDprocesses.outgoing_particles(::TestProcess) = OUTGOING_PARTICLES

        # single ps points
        p_in_phys= _rand_momenta(RNG, N_INCOMING)
        p_in_phys_invalid = _rand_momenta(RNG, N_INCOMING+1)
        p_in_unphys = _rand_momenta_failed(RNG,N_INCOMING)
        p_in_unphys_invalid = _rand_momenta_failed(RNG,N_INCOMING+1)
        p_out_phys= _rand_momenta(RNG, N_OUTGOING)
        p_out_phys_invaid = _rand_momenta(RNG, N_OUTGOING+1)

        # sets of ps points
        p_in_set_phys= _rand_momenta(RNG, N_INCOMING, 2)
        p_in_set_unphys_mix = _rand_momenta_failed_mix(RNG,N_INCOMING,2)
        p_in_set_unphys_all = _rand_momenta_failed_all(RNG,N_INCOMING,2)
        p_in_set_phys_invalid = _rand_momenta(RNG, N_INCOMING+1, 2)
        p_in_set_unphys_mix_invalid = _rand_momenta_failed_mix(RNG,N_INCOMING+1,2)
        p_in_set_unphys_all_invalid = _rand_momenta_failed_all(RNG,N_INCOMING+1,2)
        p_out_set_phys = _rand_momenta(RNG, N_OUTGOING, 2)
        p_out_set_phys_invalid = _rand_momenta(RNG, N_OUTGOING+1, 2)

        p_in_all = (p_in_phys, p_in_unphys,
            p_in_phys_invalid,p_in_unphys_invalid,
            p_in_set_phys, p_in_set_unphys_mix, p_in_set_unphys_all,
            p_in_set_phys_invalid, p_in_set_unphys_mix_invalid, p_in_set_unphys_all_invalid)

        p_out_all = (p_out_phys,
            p_out_phys_invaid,
            p_out_set_phys,
            p_out_set_phys_invalid)
        # all combinations
        p_combs = Iterators.product(p_in_all,p_out_all)

        p_in_all_valid = (p_in_phys, p_in_unphys,
            p_in_set_phys, p_in_set_unphys_mix, p_in_set_unphys_all)

        p_out_all_valid = (p_out_phys, p_out_set_phys)
        
        # all valid combinations
        p_combs_valid = collect(Iterators.product(p_in_all_valid,p_out_all_valid))

        p_in_all_phys = (p_in_phys,p_in_set_phys)
        p_out_all_phys = (p_out_phys,p_out_set_phys)

        p_combs_phys = Iterators.product(p_in_all_phys, p_out_all_phys)
        
    @testset "cross section" begin
        @testset "unsafe" begin
            @testset "compute" begin
                for (P_IN,P_OUT) in p_combs_phys
                    diffCS = unsafe_differential_cross_section(
                        TestProcess(),
                        TestModel(),
                        TestPhasespaceDef(),
                        P_IN,
                        TestPhasespaceDef(),
                        P_OUT
                    )
                    groundtruth = _groundtruth_unsafe_diffCS(TestProcess(), P_IN, P_OUT)
                    @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "invalid input" begin
                for (P_IN,P_OUT) in p_combs

                    # filter out all valid combinations
                    if !((P_IN,P_OUT) in p_combs_valid)
                        @test_throws DimensionMismatch unsafe_differential_cross_section(
                            TestProcess(),
                            TestModel(),
                            TestPhasespaceDef(),
                            P_IN,
                            TestPhasespaceDef(),
                            P_OUT
                        )

                    end

                end
            end
        end
        @testset "safe" begin
            @testset "compute" begin

                for (P_IN,P_OUT) in p_combs_valid 
                    diffCS = differential_cross_section(
                        TestProcess(),
                        TestModel(),
                        TestPhasespaceDef(),
                        P_IN,
                        TestPhasespaceDef(),
                        P_OUT
                    )
                    groundtruth = _groundtruth_safe_diffCS(TestProcess(), P_IN, P_OUT)
                    @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end


            @testset "invalid input" begin
                for (P_IN,P_OUT) in p_combs

                    # filter out all valid combinations
                    if !((P_IN,P_OUT) in p_combs_valid)
                        @test_throws DimensionMismatch differential_cross_section(
                            TestProcess(),
                            TestModel(),
                            TestPhasespaceDef(),
                            P_IN,
                            TestPhasespaceDef(),
                            P_OUT
                        )

                    end

                end
            end
        end
        # @testset "total cross section" begin
        #     @testset "compute vector" begin
        #         totCS = total_cross_section(TestProcess(), TestModel(), p_in)
        #         groundtruth = _groundtruth_totCS(p_in)
        #         @test isapprox(totCS, groundtruth, atol=ATOL, rtol=RTOL)
        #     end
        #
        #     @testset "compute matrix" begin
        #         totCS = total_cross_section(TestProcess(), TestModel(), p_in_set)
        #
        #         groundtruth = Vector{QEDprocesses._base_component_type(p_in)}(
        #             undef, size(p_in_set, 2)
        #         )
        #         for i in 1:size(p_in_set, 2)
        #             groundtruth[i] = _groundtruth_totCS(view(p_in_set, :, i))
        #         end
        #         @test isapprox(totCS, groundtruth, atol=ATOL, rtol=RTOL)
        #     end
        #
        #     @testset "fail vector" begin
        #         @test_throws DimensionMismatch total_cross_section(
        #             TestProcess(), TestModel(), _rand_momenta(RNG, N_INCOMING + 1)
        #         )
        #     end
        #
        #     @testset "fail matrix" begin
        #         @test_throws DimensionMismatch total_cross_section(
        #             TestProcess(), TestModel(), _rand_momenta(RNG, N_INCOMING + 1, 2)
        #         )
        #     end
        # end
    end
    
    @testset "probability" begin
            @testset "unsafe compute" begin
                for (P_IN,P_OUT) in p_combs_phys
                    prob = unsafe_probability(
                        TestProcess(),
                        TestModel(),
                        TestPhasespaceDef(),
                        P_IN,
                        TestPhasespaceDef(),
                        P_OUT
                    )
                    groundtruth = _groundtruth_unsafe_probability(TestProcess(), P_IN, P_OUT)
                    @test isapprox(prob, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end

            @testset "unsafe invalid input" begin
                for (P_IN,P_OUT) in p_combs

                    # filter out all valid combinations
                    if !((P_IN,P_OUT) in p_combs_valid)
                        @test_throws DimensionMismatch unsafe_probability(
                            TestProcess(),
                            TestModel(),
                            TestPhasespaceDef(),
                            P_IN,
                            TestPhasespaceDef(),
                            P_OUT
                        )

                    end

                end
            end
            @testset "safe compute" begin

                for (P_IN,P_OUT) in p_combs_valid 
                    prob = probability(
                        TestProcess(),
                        TestModel(),
                        TestPhasespaceDef(),
                        P_IN,
                        TestPhasespaceDef(),
                        P_OUT
                    )
                    groundtruth = _groundtruth_safe_probability(TestProcess(), P_IN, P_OUT)
                    @test isapprox(prob, groundtruth, atol=ATOL, rtol=RTOL)
                end
            end


            @testset "safe invalid input" begin
                for (P_IN,P_OUT) in p_combs

                    # filter out all valid combinations
                    if !((P_IN,P_OUT) in p_combs_valid)
                        @test_throws DimensionMismatch probability(
                            TestProcess(),
                            TestModel(),
                            TestPhasespaceDef(),
                            P_IN,
                            TestPhasespaceDef(),
                            P_OUT
                        )

                    end

                end
            end
    end
end
