using Random
using Suppressor
using QEDbase
using QEDprocesses

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

include("../test_implementation/TestImplementation.jl")

@testset "($N_INCOMING,$N_OUTGOING)" for (N_INCOMING, N_OUTGOING) in Iterators.product(
    (1, rand(RNG, 2:8)), (1, rand(RNG, 2:8))
)
    INCOMING_PARTICLES = rand(RNG, TestImplementation.PARTICLE_SET, N_INCOMING)
    OUTGOING_PARTICLES = rand(RNG, TestImplementation.PARTICLE_SET, N_OUTGOING)

    TESTPROC = TestImplementation.TestProcess(INCOMING_PARTICLES, OUTGOING_PARTICLES)
    TESTPROC_FAIL = TestImplementation.TestProcess_FAIL(
        INCOMING_PARTICLES, OUTGOING_PARTICLES
    )
    TESTMODEL = TestImplementation.TestModel()
    TESTMODEL_FAIL = TestImplementation.TestModel_FAIL()
    TESTPSDEF = TestImplementation.TestPhasespaceDef()
    TESTPSDEF_FAIL = TestImplementation.TestPhasespaceDef_FAIL()
    IN_PS = TestImplementation._rand_momenta(RNG, N_INCOMING)
    OUT_PS = TestImplementation._rand_momenta(RNG, N_OUTGOING)

    @testset "failed interface" begin
        @testset "failed process interface" begin
            @test_throws MethodError incoming_particles(TESTPROC_FAIL)
            @test_throws MethodError outgoing_particles(TESTPROC_FAIL)
        end
        @testset "$PROC $MODEL" for (PROC, MODEL) in Iterators.product(
            (TESTPROC, TESTPROC_FAIL), (TESTMODEL, TESTMODEL_FAIL)
        )
            in_ps = TestImplementation._rand_momenta(RNG, 2)
            out_ps = TestImplementation._rand_momenta(RNG, 2)
            if TestImplementation._any_fail(PROC, MODEL)
                @test_throws MethodError QEDprocesses._incident_flux(PROC, MODEL, in_ps)
                @test_throws MethodError QEDprocesses._averaging_norm(PROC, MODEL)
                @test_throws MethodError QEDprocesses._matrix_element(
                    PROC, MODEL, in_ps, out_ps
                )
            end

            for (IN_PS_DEF, OUT_PS_DEF) in
                Iterators.product((TESTPSDEF, TESTPSDEF_FAIL), (TESTPSDEF, TESTPSDEF_FAIL))
                if TestImplementation._any_fail(PROC, MODEL, IN_PS_DEF, OUT_PS_DEF)
                    @test_throws MethodError QEDprocesses._phase_space_factor(
                        PROC, MODEL, IN_PS_DEF, in_ps, OUT_PS_DEF, out_ps
                    )
                end
            end
        end
    end
    @testset "incoming/outgoing particles" begin
        @test incoming_particles(TESTPROC) == INCOMING_PARTICLES
        @test outgoing_particles(TESTPROC) == OUTGOING_PARTICLES
        @test number_incoming_particles(TESTPROC) == N_INCOMING
        @test number_outgoing_particles(TESTPROC) == N_OUTGOING
    end

    @testset "incident flux" begin
        test_incident_flux = QEDprocesses._incident_flux(TESTPROC, TESTMODEL, IN_PS)
        groundtruth = TestImplementation._groundtruth_incident_flux(IN_PS)
        @test isapprox(test_incident_flux, groundtruth, atol=ATOL, rtol=RTOL)
    end

    @testset "averaging norm" begin
        test_avg_norm = QEDprocesses._averaging_norm(TESTPROC)
        groundtruth = TestImplementation._groundtruth_averaging_norm(TESTPROC)
        @test isapprox(test_avg_norm, groundtruth, atol=ATOL, rtol=RTOL)
    end

    @testset "matrix element" begin
        test_matrix_element = QEDprocesses._matrix_element(
            TESTPROC, TESTMODEL, IN_PS, OUT_PS
        )
        groundtruth = TestImplementation._groundtruth_matrix_element(IN_PS, OUT_PS)
        @test length(test_matrix_element) == length(groundtruth)
        for i in eachindex(test_matrix_element)
            @test isapprox(test_matrix_element[i], groundtruth[i], atol=ATOL, rtol=RTOL)
        end
    end

    @testset "is in phasespace" begin
        @test QEDprocesses._is_in_phasespace(
            TESTPROC, TESTMODEL, TESTPSDEF, IN_PS, TESTPSDEF, OUT_PS
        )

        IN_PS_unphysical = deepcopy(IN_PS)
        IN_PS_unphysical[1] = SFourMomentum(zeros(4))

        @test !QEDprocesses._is_in_phasespace(
            TESTPROC, TESTMODEL, TESTPSDEF, IN_PS_unphysical, TESTPSDEF, OUT_PS
        )
    end

    @testset "is in phasespace" begin
        @test QEDprocesses._is_in_phasespace(
            TESTPROC, TESTMODEL, TESTPSDEF, IN_PS, TESTPSDEF, OUT_PS
        )

        IN_PS_unphysical = deepcopy(IN_PS)
        IN_PS_unphysical[1] = SFourMomentum(zeros(4))

        @test !QEDprocesses._is_in_phasespace(
            TESTPROC, TESTMODEL, TESTPSDEF, IN_PS_unphysical, TESTPSDEF, OUT_PS
        )

        OUT_PS_unphysical = deepcopy(OUT_PS)
        OUT_PS_unphysical[end] = ones(SFourMomentum)

        @test !QEDprocesses._is_in_phasespace(
            TESTPROC, TESTMODEL, TESTPSDEF, IN_PS, TESTPSDEF, OUT_PS_unphysical
        )
        @test !QEDprocesses._is_in_phasespace(
            TESTPROC, TESTMODEL, TESTPSDEF, IN_PS_unphysical, TESTPSDEF, OUT_PS_unphysical
        )
    end

    @testset "phase space factor" begin
        test_phase_space_factor = QEDprocesses._phase_space_factor(
            TESTPROC, TESTMODEL, TESTPSDEF, IN_PS, TESTPSDEF, OUT_PS
        )
        groundtruth = TestImplementation._groundtruth_phase_space_factor(IN_PS, OUT_PS)
        @test isapprox(test_phase_space_factor, groundtruth, atol=ATOL, rtol=RTOL)
=======
    IN_PS = _rand_momenta(RNG,N_INCOMING)
    OUT_PS = _rand_momenta(RNG,N_OUTGOING)

    QEDprocesses.incoming_particles(::TestProcess) = INCOMING_PARTICLES
    QEDprocesses.outgoing_particles(::TestProcess) = OUTGOING_PARTICLES
        

    @testset "incoming/outgoing particles" begin
        @test incoming_particles(TestProcess()) == INCOMING_PARTICLES
        @test outgoing_particles(TestProcess()) == OUTGOING_PARTICLES
        @test number_incoming_particles(TestProcess()) == N_INCOMING
        @test number_outgoing_particles(TestProcess()) == N_OUTGOING
    end

    @testset "incident flux" begin
        test_incident_flux = QEDprocesses._incident_flux(TestProcess(),TestModel(),IN_PS) 
        groundtruth = _groundtruth_incident_flux(IN_PS)
        @test isapprox(test_incident_flux,groundtruth,atol=ATOL,rtol=RTOL)
    end
        
    @testset "matrix element" begin
        test_avg_norm = QEDprocesses._averaging_norm(TestProcess())
        groundtruth = _groundtruth_averaging_norm(TestProcess())
        @test isapprox(test_avg_norm,groundtruth,atol=ATOL,rtol=RTOL)
    end
        
    @testset "matrix element" begin
        test_matrix_element = QEDprocesses._matrix_element(TestProcess(),TestModel(),IN_PS,OUT_PS) 
        groundtruth = _groundtruth_matrix_element(IN_PS,OUT_PS)
        @test length(test_matrix_element) == length(groundtruth)
        for i in eachindex(test_matrix_element)
            @test isapprox(test_matrix_element[i],groundtruth[i],atol=ATOL,rtol=RTOL)
        end
    end

    @testset "phase space factor" begin
        test_phase_space_factor = QEDprocesses._phase_space_factor(TestProcess(),TestModel(),TestPhasespaceDef(),IN_PS,TestPhasespaceDef(),OUT_PS) 
        groundtruth = _groundtruth_phase_space_factor(IN_PS,OUT_PS)
        @test isapprox(test_phase_space_factor,groundtruth,atol=ATOL,rtol=RTOL)
>>>>>>> a8e7d7d (updated differential cross section and tests)
    end
end
