
# dummy particles
struct TestParticleFermion <: FermionLike end
struct TestParticleBoson <: BosonLike end

const PARTICLE_SET = [TestParticleFermion(), TestParticleBoson()]

"""

    TestProcess(rng,incoming_particles,outgoing_particles)

"""
struct TestProcess{IP<:Tuple,OP<:Tuple} <: AbstractProcessDefinition
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

struct TestProcess_FAIL{IP<:Tuple,OP<:Tuple} <: AbstractProcessDefinition
    incoming_particles::IP
    outgoing_particles::OP
end

function TestProcess_FAIL(rng::AbstractRNG, N_in::Int, N_out::Int)
    in_particles = Tuple(rand(rng, PARTICLE_SET, N_in))
    out_particles = Tuple(rand(rng, PARTICLE_SET, N_out))
    return TestProcess_FAIL(in_particles, out_particles)
end

function QEDprocesses.in_phase_space_dimension(proc::TestProcess, ::TestModel)
    return number_incoming_particles(proc) * 4
end
function QEDprocesses.out_phase_space_dimension(proc::TestProcess, ::TestModel)
    return number_outgoing_particles(proc) * 4
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
    ps_def::TestPhasespaceDef,
    in_ps::AbstractVector{T},
    out_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return _groundtruth_is_in_phasespace(in_ps, out_ps)
end

function QEDprocesses._phase_space_factor(
    ::TestProcess,
    ::TestModel,
    ps_def::TestPhasespaceDef,
    in_ps::AbstractVector{T},
    out_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return _groundtruth_phase_space_factor(in_ps, out_ps)
end

function QEDprocesses._generate_incoming_momenta(
    proc::TestProcess,
    model::TestModel,
    phase_space_def::TestPhasespaceDef,
    in_phase_space::AbstractVector{T},
) where {T<:Real}
    return _groundtruth_generate_momenta(in_phase_space)
end

function QEDprocesses._generate_outgoing_momenta(
    proc::TestProcess,
    model::TestModel,
    phase_space_def::TestPhasespaceDef,
    out_phase_space::AbstractVector{T},
) where {T<:Real}
    return _groundtruth_generate_momenta(out_phase_space)
end

function QEDprocesses._total_probability(
    proc::TestProcess, model::TestModel, ps_def::TestPhasespaceDef, in_ps::AbstractVector{T}
) where {T<:QEDbase.AbstractFourMomentum}
    return _groundtruth_total_probability(in_ps)
end
