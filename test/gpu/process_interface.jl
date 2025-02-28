if isempty(GPUS)
    @info """No GPU tests are enabled, skipping tests...
    To test GPU functionality, please use 'TEST_<GPU> = 1 julia ...' for one of GPU=[CUDA, AMDGPU, METAL, ONEAPI]"""
    return nothing
end

using QEDprocesses
using QEDbase
using QEDcore

using Random
using SafeTestsets

DEF_POLS = (PolX(), PolY())
DEF_SPINS = (SpinUp(), SpinDown())

const MODEL = PerturbativeQED()
const IN_PSL = ComptonRestSystem(Energy(2))
const OUT_PSL = ComptonSphericalLayout(IN_PSL)

PROC_DEF_TUPLES = [
    (Compton(), MODEL, OUT_PSL),
    [
        (Compton(s1, p1, s2, p2), MODEL, OUT_PSL) for
        (s1, p1, s2, p2) in Iterators.product(DEF_SPINS, DEF_POLS, DEF_SPINS, DEF_POLS)
    ]...,
]

RNG = Random.MersenneTwister(573)

@testset "Testing with $GPU_MODULE" for (GPU_MODULE, VECTOR_TYPE) in GPUS
    @testset "$proc $model $ps_def" for (proc, model, ps_def) in PROC_DEF_TUPLES
        N = 100

        @info "Testing $proc $model $ps_def"
        flush(stdout)

        psps = [
            PhaseSpacePoint(
                proc, model, ps_def, _rand_coordinates(RNG, proc, model, ps_def)...
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
                    Vector(QEDbase._averaging_norm.(gpuprocs)),
                    QEDbase._averaging_norm.(procs),
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

            @test Vector(QEDbase.in_phase_space_dimension.(gpuprocs, model)) ==
                QEDbase.in_phase_space_dimension.(procs, model)
            @test Vector(QEDbase.out_phase_space_dimension.(gpuprocs, model)) ==
                QEDbase.out_phase_space_dimension.(procs, model)
        end

        @testset "Private PSP/Process Interface" begin
            @test all(
                isapprox.(
                    Vector(QEDbase._incident_flux.(gpupsps)), QEDbase._incident_flux.(psps)
                ),
            )

            @test all(
                tuple_isapprox.(
                    Vector(QEDbase._matrix_element.(gpupsps)),
                    QEDbase._matrix_element.(psps);
                    rtol=sqrt(eps(Float64)),
                ),
            )

            @test Vector(QEDbase._is_in_phasespace.(gpupsps)) ==
                QEDbase._is_in_phasespace.(psps)

            @test all(
                isapprox.(
                    Vector(QEDbase._phase_space_factor.(gpupsps)),
                    QEDbase._phase_space_factor.(psps),
                ),
            )

            # this currently throws an exception because QuadGK does not work on the GPU
            @test all(
                isapprox.(
                    Vector(QEDprocesses._total_probability.(gpupsps)),
                    QEDprocesses._total_probability.(psps),
                ),
            ) broken = true
        end

        @testset "Public PSP/Process Interface" begin
            @test all(
                isapprox.(
                    Vector(differential_probability.(gpupsps)),
                    differential_probability.(psps),
                ),
            )

            @test all(
                isapprox.(
                    Vector(QEDbase._is_in_phasespace.(gpupsps)),
                    QEDbase._is_in_phasespace.(psps),
                ),
            )

            @test all(
                isapprox.(
                    Vector(differential_cross_section.(gpupsps)),
                    differential_cross_section.(psps),
                ),
            )

            # as above, this currently throws an exception because QuadGK does not work on the GPU
            @test all(
                isapprox.(Vector(total_cross_section.(gpupsps)), total_cross_section.(psps))
            ) broken = true
        end
    end
end
