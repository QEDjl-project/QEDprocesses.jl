using Random
using QEDbase
using QEDprocesses

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

include("test_implementation.jl")

@testset "($N_INCOMING,$N_OUTGOING)" for (N_INCOMING, N_OUTGOING) in Iterators.product(
    (1, rand(RNG, 2:8)), (1, rand(RNG, 2:8))
)
    INCOMING_PARTICLES = rand(RNG, PARTICLE_SET, N_INCOMING)
    OUTGOING_PARTICLES = rand(RNG, PARTICLE_SET, N_OUTGOING)

    QEDprocesses.incoming_particles(::TestProcess) = INCOMING_PARTICLES
    QEDprocesses.outgoing_particles(::TestProcess) = OUTGOING_PARTICLES

    @testset "cross section" begin
        p_in = _rand_momenta(RNG, N_INCOMING)
        p_in_unphysical = deepcopy(p_in)
        p_in_unphysical[1] = SFourMomentum(zeros(4))
        p_out = _rand_momenta(RNG, N_OUTGOING)
        p_in_set = _rand_momenta(RNG, N_INCOMING, 2)
        p_in_set_unphys_mix = deepcopy(p_in_set)
        p_in_set_unphys_mix[1,1] = SFourMomentum(zeros(4))
        p_in_set_unphys_all = deepcopy(p_in_set_unphys_mix)
        p_in_set_unphys_all[1,2] = SFourMomentum(zeros(4))
        p_out_set = _rand_momenta(RNG, N_OUTGOING, 2)

        @testset "unsafe differential cross section" begin
            @testset "compute vector-vector" begin
                diffCS = unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    p_out,
                )
                groundtruth = _groundtruth_unsafe_diffCS(TestProcess(), p_in, p_out)
                @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "compute vector-matrix" begin
                diffCS = unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    p_out_set,
                )

                groundtruth = Vector{QEDprocesses._base_component_type(p_in)}(
                    undef, size(p_out_set, 2)
                )
                for i in 1:size(p_out_set, 2)
                    groundtruth[i] = _groundtruth_unsafe_diffCS(
                        TestProcess(), p_in, view(p_out_set, :, i)
                    )
                end
                @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "compute matrix-vector" begin
                diffCS = unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    p_out,
                )
                groundtruth = Vector{QEDprocesses._base_component_type(p_in_set)}(
                    undef, size(p_in_set, 2)
                )
                for i in 1:size(p_in_set, 2)
                    groundtruth[i] = _groundtruth_unsafe_diffCS(
                        TestProcess(), view(p_in_set, :, i), p_out
                    )
                end
                @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "compute matrix-matrix" begin
                diffCS = unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    p_out_set,
                )
                groundtruth = Matrix{QEDprocesses._base_component_type(p_in_set)}(
                    undef, size(p_in_set, 2), size(p_out_set, 2)
                )
                for i in 1:size(p_in_set, 2)
                    for j in 1:size(p_out_set, 2)
                        groundtruth[i, j] = _groundtruth_unsafe_diffCS(
                            TestProcess(), view(p_in_set, :, i), view(p_out_set, :, j)
                        )
                    end
                end
                @test isapprox(diffCS, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "fail vector-vector" begin
                @test_throws DimensionMismatch unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1),
                    TestPhasespaceDef(),
                    p_out,
                )
                @test_throws DimensionMismatch unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1),
                )
            end

            @testset "fail vector-matrix" begin
                @test_throws DimensionMismatch unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1),
                    TestPhasespaceDef(),
                    p_out_set,
                )
                @test_throws DimensionMismatch unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )
            end

            @testset "fail matrix-vector" begin
                @test_throws DimensionMismatch unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                    TestPhasespaceDef(),
                    p_out,
                )
                @test_throws DimensionMismatch unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1),
                )
            end

            @testset "fail matrix-matrix" begin
                @test_throws DimensionMismatch unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                    TestPhasespaceDef(),
                    p_out_set,
                )
                @test_throws DimensionMismatch unsafe_differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )
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
end
