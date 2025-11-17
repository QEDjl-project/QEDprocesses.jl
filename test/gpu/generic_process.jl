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

const IN_PSL_GENERIC = TwoBodyRestSystem()
const OUT_PSL_GENERIC = FlatPhaseSpaceLayout(IN_PSL_GENERIC)

PROC_DEF_TUPLES = [
    (ScatteringProcess((Electron(), Photon()), (Electron(), Photon())), MODEL, OUT_PSL_GENERIC),
    (ScatteringProcess((Electron(), Positron()), (Electron(), Positron())), MODEL, OUT_PSL_GENERIC),
]

RNG = Random.MersenneTwister(573)

if !LARGE_TESTS()
    @info "Skipping large tests...\nEnable them explicitly with an environment variable LARGE_TESTS=1"
    PROC_DEF_TUPLES = PROC_DEF_TUPLES[1:5]
end

@testset "Testing with $GPU_MODULE" for (GPU_MODULE, VECTOR_TYPE) in GPUS
    @testset "Float type $FLOAT_T" for FLOAT_T in GPU_FLOAT_TYPES[GPU_MODULE]
        @testset "$proc ($(incoming_spin_pols(proc)), $(outgoing_spin_pols(proc)))" for (proc, model, psl) in PROC_DEF_TUPLES
            N = 128

            @info "Testing $proc ($(incoming_spin_pols(proc)), $(outgoing_spin_pols(proc))) ($FLOAT_T)"
            flush(stdout)

            coords = [_rand_coordinates(RNG, proc, model, psl, FLOAT_T) for _ in 1:N]

            psps = PhaseSpacePoint.(
                proc, model, Ref(psl), getindex.(coords, 1), getindex.(coords, 2)
            )

            procs = [proc for _ in 1:N]

            gpupsps = VECTOR_TYPE(psps)
            gpuprocs = VECTOR_TYPE(procs)

            @testset "PSP interface" begin
                in_moms_gpu = Vector(momenta.(gpupsps, Incoming()))
                out_moms_gpu = Vector(momenta.(gpupsps, Outgoing()))
                in_moms = momenta.(psps, Incoming())
                out_moms = momenta.(psps, Outgoing())

                @test eltype(eltype(eltype(in_moms_gpu))) == FLOAT_T
                @test eltype(eltype(eltype(out_moms_gpu))) == FLOAT_T
                @test eltype(eltype(eltype(in_moms))) == FLOAT_T
                @test eltype(eltype(eltype(out_moms))) == FLOAT_T

                @test getindex.(in_moms_gpu, Ref(1)) == getindex.(in_moms, Ref(1))
                @test getindex.(in_moms_gpu, Ref(2)) == getindex.(in_moms, Ref(2))
                @test getindex.(out_moms_gpu, Ref(1)) == getindex.(out_moms, Ref(1))
                @test getindex.(out_moms_gpu, Ref(2)) == getindex.(out_moms, Ref(2))
            end

            @testset "Private Process Functions" begin
                # TODO: this isn't very pretty but necessary so the return type is stable
                let FLOAT_T = FLOAT_T
                    function wrap(::Type{T}) where {T}
                        function _stable_norm(proc)
                            return QEDbase._averaging_norm(T, proc)
                        end

                        @test sum(
                            isapprox.(Vector(_stable_norm.(gpuprocs)), _stable_norm.(procs))
                        ) == N
                    end
                    wrap(FLOAT_T)
                end
            end

            @testset "Public Process Functions" begin
                @test sum(Vector(incoming_particles.(gpuprocs)) .== incoming_particles.(procs)) == N
                @test sum(Vector(outgoing_particles.(gpuprocs)) .== outgoing_particles.(procs)) == N

                @test sum(
                    Vector(particles.(gpuprocs, Incoming())) .==
                        particles.(procs, Incoming())
                ) == N
                @test sum(
                    Vector(particles.(gpuprocs, Outgoing())) .==
                        particles.(procs, Outgoing())
                ) == N

                @test sum(
                    Vector(number_incoming_particles.(gpuprocs)) .==
                        number_incoming_particles.(procs)
                ) == N
                @test sum(
                    Vector(number_outgoing_particles.(gpuprocs)) .==
                        number_outgoing_particles.(procs)
                ) == N

                @test sum(
                    Vector(number_particles.(gpuprocs, Incoming())) .==
                        number_particles.(procs, Incoming())
                ) == N
                @test sum(
                    Vector(number_particles.(gpuprocs, Outgoing())) .==
                        number_particles.(procs, Outgoing())
                ) == N

                @test sum(
                    Vector(QEDbase.in_phase_space_dimension.(gpuprocs, model)) .==
                        QEDbase.in_phase_space_dimension.(procs, model)
                ) == N
                @test sum(
                    Vector(QEDbase.out_phase_space_dimension.(gpuprocs, model)) .==
                        QEDbase.out_phase_space_dimension.(procs, model)
                ) == N
            end

            @testset "Private PSP/Process Interface" begin
                gpu = Vector(QEDbase._incident_flux.(gpupsps))
                gt = QEDbase._incident_flux.(psps)
                @test eltype(gpu) == FLOAT_T
                @test sum(isapprox.(gpu, gt)) == N

                gpu = Vector(QEDbase._is_in_phasespace.(gpupsps))
                gt = QEDbase._is_in_phasespace.(psps)
                @test eltype(gpu) == Bool
                @test sum(gpu .== gt) == N

                gpu = Vector(QEDbase._phase_space_factor.(gpupsps))
                gt = QEDbase._phase_space_factor.(psps)
                @test eltype(gpu) == FLOAT_T
                @test sum(isapprox.(gpu, gt)) == N
            end

            @testset "KernelAbstractions Probability" begin
                dest = similar(gpupsps, FLOAT_T)
                gt = unsafe_differential_probability.(psps)
                unsafe_differential_probability!(dest, gpupsps)
                @test sum(isapprox.(Vector(dest), gt)) == N

                #=
                fill!(dest, zero(FLOAT_T))
                gt = differential_probability.(psps)
                differential_probability!(dest, gpupsps)
                @test sum(isapprox.(Vector(dest), gt)) == N
                =#
            end

            @testset "KernelAbstractions Cross Section" begin
                dest = similar(gpupsps, FLOAT_T)
                gt = unsafe_differential_cross_section.(psps)
                unsafe_differential_cross_section!(dest, gpupsps)
                @test sum(isapprox.(Vector(dest), gt)) == N

                #=
                fill!(dest, zero(FLOAT_T))
                gt = differential_cross_section.(psps)
                differential_cross_section!(dest, gpupsps)
                @test sum(isapprox.(Vector(dest), gt)) == N
                =#
            end
        end
    end
end
