
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

    @testset "probabilities" begin
        p_in = _rand_momenta(RNG, N_INCOMING)
        p_out = _rand_momenta(RNG, N_OUTGOING)
        p_in_set = _rand_momenta(RNG, N_INCOMING, 2)
        p_out_set = _rand_momenta(RNG, N_OUTGOING, 2)

        @testset "unsafe differential probability" begin
            @testset "compute vector-vector" begin
                diffProb = unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    p_out,
                )
                groundtruth = _groundtruth_unsafe_probability(TestProcess(), p_in, p_out)
                @test isapprox(diffProb, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "compute vector-matrix" begin
                diffProb = unsafe_probability(
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
                    groundtruth[i] = _groundtruth_unsafe_probability(
                        TestProcess(), p_in, view(p_out_set, :, i)
                    )
                end
                @test isapprox(diffProb, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "compute matrix-vector" begin
                diffProb = unsafe_probability(
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
                    groundtruth[i] = _groundtruth_unsafe_probability(
                        TestProcess(), view(p_in_set, :, i), p_out
                    )
                end
                @test isapprox(diffProb, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "compute matrix-matrix" begin
                diffProb = unsafe_probability(
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
                        groundtruth[i, j] = _groundtruth_unsafe_probability(
                            TestProcess(), view(p_in_set, :, i), view(p_out_set, :, j)
                        )
                    end
                end
                @test isapprox(diffProb, groundtruth, atol=ATOL, rtol=RTOL)
            end

            @testset "fail vector-vector" begin
                @test_throws DimensionMismatch unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1),
                    TestPhasespaceDef(),
                    p_out,
                )
                @test_throws DimensionMismatch unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1),
                )
            end

            @testset "fail vector-matrix" begin
                @test_throws DimensionMismatch unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1),
                    TestPhasespaceDef(),
                    p_out_set,
                )
                @test_throws DimensionMismatch unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )
            end

            @testset "fail matrix-vector" begin
                @test_throws DimensionMismatch unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                    TestPhasespaceDef(),
                    p_out,
                )
                @test_throws DimensionMismatch unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1),
                )
            end

            @testset "fail matrix-matrix" begin
                @test_throws DimensionMismatch unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                    TestPhasespaceDef(),
                    p_out_set,
                )
                @test_throws DimensionMismatch unsafe_probability(
                    TestProcess(),
                    TestModel(),
                    TestPhasespaceDef(),
                    p_in_set,
                    TestPhasespaceDef(),
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )
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
