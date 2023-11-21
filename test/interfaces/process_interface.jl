using Random
using QEDbase
using QEDprocesses

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

include("../test_implementation.jl")

@testset "interface fail" for (PROC, MODEL) in Iterators.product(
    (TestProcess(), TestProcess_FAIL()), (TestModel(), TestModel_FAIL())
)
    in_ps = _rand_momenta(RNG, 2)
    out_ps = _rand_momenta(RNG, 2)
    if _any_fail(PROC, MODEL)
        @test_throws MethodError incoming_particles(PROC)
        @test_throws MethodError outgoing_particles(PROC)
        @test_throws MethodError QEDprocesses._incident_flux(PROC, MODEL, in_ps)
        @test_throws MethodError QEDprocesses._matrix_element(PROC, MODEL, in_ps, out_ps)
    end

    for (IN_PS_DEF, OUT_PS_DEF) in Iterators.product(
        (TestPhasespaceDef(), TestPhasespaceDef_FAIL()),
        (TestPhasespaceDef(), TestPhasespaceDef_FAIL()),
    )
        if _any_fail(PROC, MODEL, IN_PS_DEF, OUT_PS_DEF)
            @test_throws MethodError QEDprocesses._phase_space_factor(
                PROC, MODEL, IN_PS_DEF, in_ps, OUT_PS_DEF, out_ps
            )
        end
    end
end

@testset "($N_INCOMING,$N_OUTGOING)" for (N_INCOMING, N_OUTGOING) in Iterators.product(
    (1, rand(RNG, 2:8)), (1, rand(RNG, 2:8))
)
    INCOMING_PARTICLES = rand(RNG, PARTICLE_SET, N_INCOMING)
    OUTGOING_PARTICLES = rand(RNG, PARTICLE_SET, N_OUTGOING)

    IN_PS = _rand_momenta(RNG, N_INCOMING)
    OUT_PS = _rand_momenta(RNG, N_OUTGOING)

    QEDprocesses.incoming_particles(::TestProcess) = INCOMING_PARTICLES
    QEDprocesses.outgoing_particles(::TestProcess) = OUTGOING_PARTICLES

    @testset "incoming/outgoing particles" begin
        @test incoming_particles(TestProcess()) == INCOMING_PARTICLES
        @test outgoing_particles(TestProcess()) == OUTGOING_PARTICLES
        @test number_incoming_particles(TestProcess()) == N_INCOMING
        @test number_outgoing_particles(TestProcess()) == N_OUTGOING
    end

    @testset "incident flux" begin
        test_incident_flux = QEDprocesses._incident_flux(TestProcess(), TestModel(), IN_PS)
        groundtruth = _groundtruth_incident_flux(IN_PS)
        @test isapprox(test_incident_flux, groundtruth, atol=ATOL, rtol=RTOL)
    end

    @testset "matrix element" begin
        test_avg_norm = QEDprocesses._averaging_norm(TestProcess())
        groundtruth = _groundtruth_averaging_norm(TestProcess())
        @test isapprox(test_avg_norm, groundtruth, atol=ATOL, rtol=RTOL)
    end

    @testset "matrix element" begin
        test_matrix_element = QEDprocesses._matrix_element(
            TestProcess(), TestModel(), IN_PS, OUT_PS
        )
        groundtruth = _groundtruth_matrix_element(IN_PS, OUT_PS)
        @test length(test_matrix_element) == length(groundtruth)
        for i in eachindex(test_matrix_element)
            @test isapprox(test_matrix_element[i], groundtruth[i], atol=ATOL, rtol=RTOL)
        end
    end

    @testset "phase space factor" begin
        test_phase_space_factor = QEDprocesses._phase_space_factor(
            TestProcess(),
            TestModel(),
            TestPhasespaceDef(),
            IN_PS,
            TestPhasespaceDef(),
            OUT_PS,
        )
        groundtruth = _groundtruth_phase_space_factor(IN_PS, OUT_PS)
        @test isapprox(test_phase_space_factor, groundtruth, atol=ATOL, rtol=RTOL)
    end
end
