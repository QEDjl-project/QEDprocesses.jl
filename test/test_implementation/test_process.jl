
# dummy particles
struct TestParticle1 <: AbstractParticle end
struct TestParticle2 <: AbstractParticle end
struct TestParticle3 <: AbstractParticle end
struct TestParticle4 <: AbstractParticle end

const PARTICLE_SET = [TestParticle1(), TestParticle2(), TestParticle3(), TestParticle4()]

"""

    TestProcess(rng,incoming_particles,outgoing_particles)

"""
struct TestProcess{IP<:AbstractVector,OP<:AbstractVector} <: AbstractProcessDefinition
    incoming_particles::IP
    outgoing_particles::OP
end

function TestProcess(rng::AbstractRNG, N_in::Int, N_out::Int)
    in_particles = rand(rng, PARTICLE_SET, N_in)
    out_particles = rand(rng, PARTICLE_SET, N_out)
    return TestProcess(in_particles, out_particles)
end

QEDprocesses.incoming_particles(proc::TestProcess) = proc.incoming_particles
QEDprocesses.outgoing_particles(proc::TestProcess) = proc.outgoing_particles

struct TestProcess_FAIL{IP<:AbstractVector,OP<:AbstractVector} <: AbstractProcessDefinition
    incoming_particles::IP
    outgoing_particles::OP
end

function TestProcess_FAIL(rng::AbstractRNG, N_in::Int, N_out::Int)
    in_particles = rand(rng, PARTICLE_SET, N_in)
    out_particles = rand(rng, PARTICLE_SET, N_out)
    return TestProcess_FAIL(in_particles, out_particles)
end

# dummy phase space definition + failing phase space definition
struct TestPhasespaceDef <: AbstractPhasespaceDefinition end
struct TestPhasespaceDef_FAIL <: AbstractPhasespaceDefinition end

# dummy implementation of the process interface

function QEDprocesses._incident_flux(
    ::TestProcess, ::TestModel, in_ps::AbstractVector{T}
) where {T<:QEDbase.AbstractFourMomentum}
    return _groundtruth_incident_flux(in_ps)
end

function QEDprocesses._averaging_norm(proc::TestProcess)
    return _groundtruth_averaging_norm(proc)
end

function QEDprocesses._matrix_element(
    ::TestProcess, ::TestModel, in_ps::AbstractVector{T}, out_ps::AbstractVector{T}
) where {T<:QEDbase.AbstractFourMomentum}
    return _groundtruth_matrix_element(in_ps, out_ps)
end

function QEDprocesses._is_in_phasespace(
    ::TestProcess,
    ::TestModel,
    in_ps_def::TestPhasespaceDef,
    in_ps::AbstractVector{T},
    out_ps_def::TestPhasespaceDef,
    out_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return _groundtruth_is_in_phasespace(in_ps, out_ps)
end

function QEDprocesses._phase_space_factor(
    ::TestProcess,
    ::TestModel,
    in_ps_def::TestPhasespaceDef,
    in_ps::AbstractVector{T},
    out_ps_def::TestPhasespaceDef,
    out_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return _groundtruth_phase_space_factor(in_ps, out_ps)
end

function QEDprocesses._total_probability(
    proc::TestProcess,
    model::TestModel,
    in_ps_def::TestPhasespaceDef,
    in_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return _groundtruth_total_probability(in_ps)
end
