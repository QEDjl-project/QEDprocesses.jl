using Random
using QEDbase
using QEDprocesses

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

function _rand_momenta(rng::AbstractRNG, N)
    moms = Vector{SFourMomentum}(undef, N)
    for i = 1:N
        moms[i] = SFourMomentum(rand(rng, 4))
    end
    return moms
end

function _rand_momenta(rng::AbstractRNG, N1, N2)
    moms = Matrix{SFourMomentum}(undef, N1, N2)
    for i = 1:N1
        for j = 1:N2
            moms[i, j] = SFourMomentum(rand(rng, 4))
        end
    end
    return moms
end

struct TestParticle1 <: AbstractParticle end
struct TestParticle2 <: AbstractParticle end
struct TestParticle3 <: AbstractParticle end
struct TestParticle4 <: AbstractParticle end

PARTICLE_SET = [TestParticle1(), TestParticle2(), TestParticle3(), TestParticle4()]

struct TestProcess <: AbstractProcessDefinition end
struct TestProcess_FAIL <: AbstractProcessDefinition end

struct TestModel <: AbstractModelDefinition end
struct TestModel_FAIL <: AbstractModelDefinition end

_groundtruth_diffCS(initPS, finalPS) = sum(initPS) * sum(finalPS)
_groundtruth_totCS(initPS) = _groundtruth_diffCS(initPS, initPS)

@testset "interface fail" begin
    @test_throws MethodError incoming_particles(TestProcess_FAIL())
    @test_throws MethodError outgoing_particles(TestProcess_FAIL())
end
@testset "($N_INCOMING,$N_OUTGOING)" for (N_INCOMING, N_OUTGOING) in Iterators.product(
    (1, rand(RNG, 2:8)),
    (1, rand(RNG, 2:8)),
)

    INCOMING_PARTICLES = rand(RNG, PARTICLE_SET, N_INCOMING)
    OUTGOING_PARTICLES = rand(RNG, PARTICLE_SET, N_OUTGOING)

    QEDprocesses.incoming_particles(::TestProcess) = INCOMING_PARTICLES
    QEDprocesses.outgoing_particles(::TestProcess) = OUTGOING_PARTICLES

    function QEDprocesses._differential_cross_section(
        proc::TestProcess,
        model::TestModel,
        in_phasespace::AbstractVector{T},
        out_phasespace::AbstractVector{T},
    ) where {T<:QEDbase.AbstractFourMomentum}
        _groundtruth_diffCS(in_phasespace, out_phasespace)
    end

    function QEDprocesses._total_cross_section(
        proc::TestProcess,
        model::TestModel,
        in_phasespace::AbstractVector{T},
    ) where {T<:QEDbase.AbstractFourMomentum}
        _groundtruth_totCS(in_phasespace)
    end

    @testset "hard interface" begin
        @test incoming_particles(TestProcess()) == INCOMING_PARTICLES
        @test outgoing_particles(TestProcess()) == OUTGOING_PARTICLES
    end

    @testset "delegated functions" begin
        @test number_incoming_particles(TestProcess()) == N_INCOMING
        @test number_outgoing_particles(TestProcess()) == N_OUTGOING
    end


    @testset "cross section" begin

        p_in = _rand_momenta(RNG, N_INCOMING)
        p_out = _rand_momenta(RNG, N_OUTGOING)
        p_in_set = _rand_momenta(RNG, N_INCOMING, 2)
        p_out_set = _rand_momenta(RNG, N_OUTGOING, 2)

        @testset "interface fail" begin
            @test_throws MethodError differential_cross_section(
                TestProcess(),
                TestModel_FAIL(),
                p_in,
                p_out,
            )
            @test_throws MethodError total_cross_section(
                TestProcess(),
                TestModel_FAIL(),
                p_in,
                p_out,
            )
        end

        @testset "differential cross section" begin
            @testset "compute vector-vector" begin
                diffCS = differential_cross_section(TestProcess(), TestModel(), p_in, p_out)
                groundtruth = _groundtruth_diffCS(p_in, p_out)
                @test isapprox(diffCS, groundtruth, atol = ATOL, rtol = RTOL)
            end

            @testset "compute vector-matrix" begin
                diffCS =
                    differential_cross_section(TestProcess(), TestModel(), p_in, p_out_set)

                groundtruth = Vector{QEDprocesses._base_component_type(p_in)}(
                    undef,
                    size(p_out_set, 2),
                )
                for i = 1:size(p_out_set, 2)
                    groundtruth[i] = _groundtruth_diffCS(p_in, view(p_out_set, :, i))
                end
                @test isapprox(diffCS, groundtruth, atol = ATOL, rtol = RTOL)
            end

            @testset "compute matrix-vector" begin
                diffCS =
                    differential_cross_section(TestProcess(), TestModel(), p_in_set, p_out)
                groundtruth = Vector{QEDprocesses._base_component_type(p_in_set)}(
                    undef,
                    size(p_in_set, 2),
                )
                for i = 1:size(p_in_set, 2)
                    groundtruth[i] = _groundtruth_diffCS(view(p_in_set, :, i), p_out)
                end
                @test isapprox(diffCS, groundtruth, atol = ATOL, rtol = RTOL)
            end

            @testset "compute matrix-matrix" begin
                diffCS = differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    p_in_set,
                    p_out_set,
                )
                groundtruth = Matrix{QEDprocesses._base_component_type(p_in_set)}(
                    undef,
                    size(p_in_set, 2),
                    size(p_out_set, 2),
                )
                for i = 1:size(p_in_set, 2)
                    for j = 1:size(p_out_set, 2)
                        groundtruth[i, j] =
                            _groundtruth_diffCS(view(p_in_set, :, i), view(p_out_set, :, j))
                    end
                end
                @test isapprox(diffCS, groundtruth, atol = ATOL, rtol = RTOL)
            end

            @testset "fail vector-vector" begin
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    _rand_momenta(RNG, N_INCOMING + 1),
                    p_out,
                )
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    p_in,
                    _rand_momenta(RNG, N_OUTGOING + 1),
                )
            end

            @testset "fail vector-matrix" begin
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    _rand_momenta(RNG, N_INCOMING + 1),
                    p_out_set,
                )
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    p_in,
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )
            end

            @testset "fail matrix-vector" begin
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                    p_out,
                )
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    p_in_set,
                    _rand_momenta(RNG, N_OUTGOING + 1),
                )
            end

            @testset "fail matrix-matrix" begin
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                    p_out_set,
                )
                @test_throws DimensionMismatch differential_cross_section(
                    TestProcess(),
                    TestModel(),
                    p_in_set,
                    _rand_momenta(RNG, N_OUTGOING + 1, 2),
                )
            end
        end
        @testset "total cross section" begin
            @testset "compute vector" begin
                totCS = total_cross_section(TestProcess(), TestModel(), p_in)
                groundtruth = _groundtruth_totCS(p_in)
                @test isapprox(totCS, groundtruth, atol = ATOL, rtol = RTOL)
            end

            @testset "compute matrix" begin
                totCS = total_cross_section(TestProcess(), TestModel(), p_in_set)

                groundtruth = Vector{QEDprocesses._base_component_type(p_in)}(
                    undef,
                    size(p_in_set, 2),
                )
                for i = 1:size(p_in_set, 2)
                    groundtruth[i] = _groundtruth_totCS(view(p_in_set, :, i))
                end
                @test isapprox(totCS, groundtruth, atol = ATOL, rtol = RTOL)
            end

            @testset "fail vector" begin
                @test_throws DimensionMismatch total_cross_section(
                    TestProcess(),
                    TestModel(),
                    _rand_momenta(RNG, N_INCOMING + 1),
                )
            end

            @testset "fail matrix" begin
                @test_throws DimensionMismatch total_cross_section(
                    TestProcess(),
                    TestModel(),
                    _rand_momenta(RNG, N_INCOMING + 1, 2),
                )
            end
        end
    end
end
