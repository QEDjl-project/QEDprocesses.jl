using QEDprocesses
using QEDbase
using AMDGPU
using CUDA
using Random
using SafeTestsets

include("../test_implementation/TestImplementation.jl")

GPU_VECTOR_TYPES = Vector{Type}()
if AMDGPU.functional()
    println("Testing with AMDGPU.jl")
    push!(GPU_VECTOR_TYPES, ROCVector)
end
if CUDA.functional()
    println("Testing with CUDA.jl")
    push!(GPU_VECTOR_TYPES, CuVector)
end

if isempty(GPU_VECTOR_TYPES)
    println("No functional GPUs found!")
    return nothing
end

PROC_DEF_TUPLES = [
    (
        Compton(),
        PerturbativeQED(),
        PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame()),
    ),
    (
        Compton(SpinUp(), PolX(), SpinDown(), PolY()),
        PerturbativeQED(),
        PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame()),
    ),
]

RNG = Random.MersenneTwister(573)

@testset "Testing with $VECTOR_TYPE" for VECTOR_TYPE in GPU_VECTOR_TYPES
    @testset "$proc $model $ps_def" for (proc, model, ps_def) in PROC_DEF_TUPLES
        N = 10_000

        psps = [
            PhaseSpacePoint(
                proc,
                model,
                ps_def,
                TestImplementation._rand_coordinates(RNG, proc, model, ps_def)...,
            ) for _ in 1:N
        ]
        procs = [proc for _ in 1:N]

        gpupsps = VECTOR_TYPE(psps)
        gpuprocs = VECTOR_TYPE(procs)

        @testset "PSP interface" begin
            in_moms_gpu = Vector(momenta.(gpupsps, Incoming()))
            out_moms_gpu = Vector(momenta.(gpupsps, Outgoing()))
            in_moms = momenta.(psps, Incoming())
            out_moms = momenta.(psps, Outgoing())

            @test getindex.(in_moms_gpu, Ref(1)) == getindex.(in_moms, Ref(1))
            @test getindex.(in_moms_gpu, Ref(2)) == getindex.(in_moms, Ref(2))
            @test getindex.(out_moms_gpu, Ref(1)) == getindex.(out_moms, Ref(1))
            @test getindex.(out_moms_gpu, Ref(2)) == getindex.(out_moms, Ref(2))
        end

        @testset "Private Process Functions" begin
            @test all(
                isapprox.(
                    Vector(QEDprocesses._averaging_norm.(gpuprocs)),
                    QEDprocesses._averaging_norm.(procs),
                ),
            )
        end

        @testset "Public Process Functions" begin
            @test Vector(incoming_particles.(gpuprocs)) == incoming_particles.(procs)
            @test Vector(outgoing_particles.(gpuprocs)) == outgoing_particles.(procs)

            @test Vector(particles.(gpuprocs, Incoming())) == particles.(procs, Incoming())
            @test Vector(particles.(gpuprocs, Outgoing())) == particles.(procs, Outgoing())

            @test Vector(number_incoming_particles.(gpuprocs)) ==
                number_incoming_particles.(procs)
            @test Vector(number_outgoing_particles.(gpuprocs)) ==
                number_outgoing_particles.(procs)

            @test Vector(number_particles.(gpuprocs, Incoming())) ==
                number_particles.(procs, Incoming())
            @test Vector(number_particles.(gpuprocs, Outgoing())) ==
                number_particles.(procs, Outgoing())

            @test Vector(in_phase_space_dimension.(gpuprocs, model)) ==
                in_phase_space_dimension.(procs, model)
            @test Vector(out_phase_space_dimension.(gpuprocs, model)) ==
                out_phase_space_dimension.(procs, model)
        end

        @testset "Private PhaseSpacePoint Interface" begin
            @test all(
                isapprox.(
                    Vector(QEDprocesses._incident_flux.(gpupsps)),
                    QEDprocesses._incident_flux.(psps),
                ),
            )

            @test all(
                TestImplementation.tuple_isapprox.(
                    Vector(QEDprocesses._matrix_element.(gpupsps)),
                    QEDprocesses._matrix_element.(psps),
                ),
            )

            @test Vector(QEDprocesses._is_in_phasespace.(gpupsps)) ==
                QEDprocesses._is_in_phasespace.(psps)

            @test all(
                isapprox.(
                    Vector(QEDprocesses._phase_space_factor.(gpupsps)),
                    QEDprocesses._phase_space_factor.(psps),
                ),
            )

            #= TODO: this is currently broken
            @test all(
                isapprox.(
                    Vector(QEDprocesses._total_probability.(gpupsps)),
                    QEDprocesses._total_probability.(psps),
                ),
            )=#
        end

        @testset "Public PhaseSpacePoint Interface" begin
            @test all(
                isapprox.(
                    Vector(differential_probability.(gpupsps)),
                    differential_probability.(psps),
                ),
            )

            @test all(
                isapprox.(
                    Vector(QEDprocesses._is_in_phasespace.(gpupsps)),
                    QEDprocesses._is_in_phasespace.(psps),
                ),
            )

            @test all(
                isapprox.(
                    Vector(differential_cross_section.(gpupsps)),
                    differential_cross_section.(psps),
                ),
            )

            #= TODO: this is currently broken
            @test all(
                isapprox.(Vector(total_cross_section.(gpupsps)), total_cross_section.(psps))
            )=#
        end
    end
end
