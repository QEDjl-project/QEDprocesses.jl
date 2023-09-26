
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

abstract type AbstractTestSetup <: AbstractComputationSetup end
QEDprocesses._compute(stp::AbstractTestSetup, x) = _groundthruth_compute(x)

struct TestSetupDefault <: AbstractTestSetup end

struct TestSetupCustomValidation <: AbstractTestSetup end 
QEDprocesses._input_validation(::TestSetupCustomValidation, x) = _groundthruth_input_validation(x)

struct TestSetupCustomPostComputation <: AbstractTestSetup end
QEDprocesses._post_computation(::TestSetupCustomPostComputation,x,y) = _groundtruth_post_computation(x,y)

struct TestSetupCustom <: AbstractTestSetup end 
QEDprocesses._input_validation(::TestSetupCustom, x) = _groundthruth_input_validation(x)
QEDprocesses._post_computation(::TestSetupCustom,x,y) = _groundtruth_post_computation(x,y)

# setups which fail on computation
struct TestSetupFAIL <: AbstractComputationSetup end

struct TestSetupCustomValidationFAIL <: AbstractComputationSetup end
QEDprocesses._input_validation(::TestSetupCustomValidationFAIL, x) = _groundthruth_input_validation(x)

struct TestSetupCustomPostComputationFAIL <: AbstractComputationSetup end
QEDprocesses._post_computation(::TestSetupCustomPostComputationFAIL,x,y) = _groundtruth_post_computation(x,y)

@testset "interface fail" begin
    @test_throws MethodError QEDprocesses._compute(TestSetupFAIL()) 
    @test_throws MethodError compute(TestSetupFAIL()) 

    @test_throws MethodError QEDprocesses._compute(TestSetupCustomValidationFAIL()) 
    @test_throws MethodError compute(TestSetupCustomValidationFAIL()) 

    @test_throws MethodError QEDprocesses._compute(TestSetupCustomPostComputationFAIL()) 
    @test_throws MethodError compute(TestSetupCustomPostComputationFAIL()) 
end

@testset "default interface" begin
    stp = TestSetupDefault()

    rnd_input = rand(RNG)
    rnd_output = rand(RNG)
    @test QEDprocesses._input_validation(stp,rnd_input)
    @test QEDprocesses._post_computation(stp,rnd_input,rnd_output) == rnd_output
    @test isapprox(QEDprocesses._compute(stp, rnd_input), _groundthruth_compute(rnd_input), atol=ATOL,rtol=RTOL)
    @test isapprox(compute(stp, rnd_input), _groundthruth_compute(rnd_input), atol=ATOL,rtol=RTOL)
end

@testset "custom input validation" begin
    stp = TestSetupCustomValidation()
    rnd_input = rand(RNG)
    @test QEDprocesses._input_validation(stp, _groundthruth_input_validation(rnd_input))
    @test !QEDprocesses._input_validation(stp, !_groundthruth_input_validation(rnd_input))
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

    @test QEDprocesses._input_validation(stp, _groundthruth_input_validation(rnd_input))
    @test !QEDprocesses._input_validation(stp, !_groundthruth_input_validation(rnd_input))
    @test_throws ErrorException compute(stp, _transform_to_invalid(rnd_input))
    @test isapprox(QEDprocesses._post_computation(stp,rnd_input,rnd_output), _groundtruth_post_computation(rnd_input,rnd_output))
    @test isapprox(compute(stp,rnd_input), _groundtruth_post_computation(rnd_input,_groundthruth_compute(rnd_input)))
end
