
using Random
using QEDbase
using QEDprocesses

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

_groundthruth_compute(x) = x
_groundthruth_input_validation(x) = (x>0)
_transform_to_invalid(x) = -abs(x)
_groundtruth_post_computation(x,y) = x+y

# setups for which the interface is implemented
abstract type AbstractTestSetup <: AbstractComputationSetup end
QEDprocesses._compute(stp::AbstractTestSetup, x) = _groundthruth_compute(x)

# setup with default implementations
struct TestSetupDefault <: AbstractTestSetup end

# setup with custom input validation
struct TestSetupCustomValidation <: AbstractTestSetup end 
QEDprocesses._is_valid_input(::TestSetupCustomValidation, x) = _groundthruth_input_validation(x)

# setup with custom post computation
struct TestSetupCustomPostComputation <: AbstractTestSetup end
QEDprocesses._post_computation(::TestSetupCustomPostComputation,x,y) = _groundtruth_post_computation(x,y)

# setup with custom input validation and post computation
struct TestSetupCustom <: AbstractTestSetup end 
QEDprocesses._is_valid_input(::TestSetupCustom, x) = _groundthruth_input_validation(x)
QEDprocesses._post_computation(::TestSetupCustom,x,y) = _groundtruth_post_computation(x,y)

# setup which fail on computation with default implementations
struct TestSetupFAIL <: AbstractComputationSetup end

# setup which fail on computation with custom input validation
struct TestSetupCustomValidationFAIL <: AbstractComputationSetup end
QEDprocesses._is_valid_input(::TestSetupCustomValidationFAIL, x) = _groundthruth_input_validation(x)

# setup which fail on computation with custom post computation
struct TestSetupCustomPostComputationFAIL <: AbstractComputationSetup end
QEDprocesses._post_computation(::TestSetupCustomPostComputationFAIL,x,y) = _groundtruth_post_computation(x,y)
@testset "general computation setup interface" begin
    @testset "interface fail" begin
        rnd_input = rand(RNG)
        @test_throws MethodError QEDprocesses._compute(TestSetupFAIL(), rnd_input) 
        @test_throws MethodError compute(TestSetupFAIL(), rnd_input) 

        @test_throws MethodError QEDprocesses._compute(TestSetupCustomValidationFAIL(), rnd_input) 
        @test_throws MethodError compute(TestSetupCustomValidationFAIL(), rnd_input) 
        # invalid input should be caught without throwing a MethodError
        @test_throws ErrorException compute(TestSetupCustomValidationFAIL(), _transform_to_invalid(rnd_input)) 

        @test_throws MethodError QEDprocesses._compute(TestSetupCustomPostComputationFAIL(), rnd_input) 
        @test_throws MethodError compute(TestSetupCustomPostComputationFAIL(), rnd_input) 
    end

    @testset "default interface" begin
        stp = TestSetupDefault()

        rnd_input = rand(RNG)
        rnd_output = rand(RNG)
        @test QEDprocesses._is_valid_input(stp,rnd_input)
        @test QEDprocesses._post_computation(stp,rnd_input,rnd_output) == rnd_output
        @test isapprox(QEDprocesses._compute(stp, rnd_input), _groundthruth_compute(rnd_input), atol=ATOL,rtol=RTOL)
        @test isapprox(compute(stp, rnd_input), _groundthruth_compute(rnd_input), atol=ATOL,rtol=RTOL)
    end

    @testset "custom input validation" begin
        stp = TestSetupCustomValidation()
        rnd_input = rand(RNG)
        @test QEDprocesses._is_valid_input(stp, _groundthruth_input_validation(rnd_input))
        @test !QEDprocesses._is_valid_input(stp, !_groundthruth_input_validation(rnd_input))
        @test isapprox(compute(stp, rnd_input), _groundthruth_compute(rnd_input), atol=ATOL,rtol=RTOL)
        @test_throws ErrorException compute(stp, _transform_to_invalid(rnd_input))
    end

    @testset "custom post computation" begin
        stp = TestSetupCustomPostComputation()
        rnd_input = rand(RNG)
        rnd_output = rand(RNG)
        @test isapprox(QEDprocesses._post_computation(stp,rnd_input,rnd_output), _groundtruth_post_computation(rnd_input,rnd_output))
        @test isapprox(compute(stp,rnd_input), _groundtruth_post_computation(rnd_input,_groundthruth_compute(rnd_input)))
    end

    @testset "custom input validation and post computation" begin
        stp = TestSetupCustom()
        rnd_input = rand(RNG)
        rnd_output = rand(RNG)

        @test QEDprocesses._is_valid_input(stp, _groundthruth_input_validation(rnd_input))
        @test !QEDprocesses._is_valid_input(stp, !_groundthruth_input_validation(rnd_input))
        @test_throws ErrorException compute(stp, _transform_to_invalid(rnd_input))
        @test isapprox(QEDprocesses._post_computation(stp,rnd_input,rnd_output), _groundtruth_post_computation(rnd_input,rnd_output))
        @test isapprox(compute(stp,rnd_input), _groundtruth_post_computation(rnd_input,_groundthruth_compute(rnd_input)))
    end
end
# process setup 

struct TestParticle1 <: AbstractParticle end
struct TestParticle2 <: AbstractParticle end
struct TestParticle3 <: AbstractParticle end
struct TestParticle4 <: AbstractParticle end

PARTICLE_SET = [TestParticle1(), TestParticle2(), TestParticle3(),TestParticle4()]

struct TestProcess <: AbstractProcessDefinition end
struct TestModel <: AbstractModelDefinition end

struct TestProcessSetup <: AbstractProcessSetup end
QEDprocesses.scattering_process(::TestProcessSetup) = TestProcess()
QEDprocesses.compute_model(::TestProcessSetup) = TestModel()

struct TestProcessSetupFAIL <: AbstractProcessSetup end

@testset "process setup interface" begin
    @testset "interface fail" begin
        @test_throws MethodError scattering_process(TestProcessSetupFAIL)
        @test_throws MethodError compute_model(TestProcessSetupFAIL)
    end

    @testset "hard interface" begin
        stp = TestProcessSetup()
        
        @test QEDprocesses._is_computation_setup(stp) 
        @test scattering_process(stp) == TestProcess()
        @test compute_model(stp) == TestModel()
    end

    @testset "($N_INCOMING,$N_OUTGOING)" for (N_INCOMING,N_OUTGOING) in Iterators.product(
                                                                                          (1, rand(RNG, 2:8)), (1, rand(RNG, 2:8))
                                                                                         )   
        INCOMING_PARTICLES = rand(RNG,PARTICLE_SET,N_INCOMING)
        OUTGOING_PARTICLES = rand(RNG,PARTICLE_SET,N_OUTGOING)

        QEDprocesses.incoming_particles(::TestProcess) = INCOMING_PARTICLES
        QEDprocesses.outgoing_particles(::TestProcess) = OUTGOING_PARTICLES 

        @testset "delegated functions" begin
            stp = TestProcessSetup()
            @test number_incoming_particles(stp) == N_INCOMING
            @test number_outgoing_particles(stp) == N_OUTGOING
        end

    end
end
