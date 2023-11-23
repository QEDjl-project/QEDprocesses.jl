
function _rand_momenta(rng::AbstractRNG, N)
    moms = Vector{SFourMomentum}(undef, N)
    for i in 1:N
        moms[i] = SFourMomentum(rand(rng, 4))
    end
    return moms
end

function _rand_momenta(rng::AbstractRNG, N1, N2)
    moms = Matrix{SFourMomentum}(undef, N1, N2)
    for i in 1:N1
        for j in 1:N2
            moms[i, j] = SFourMomentum(rand(rng, 4))
        end
    end
    return moms
end

function _groundtruth_incident_flux(in_ps)
    s = sum(in_ps)
    return s * s
end

function _groundtruth_matrix_element(in_ps, out_ps)
    s_in = sum(in_ps)
    s_out = sum(out_ps)
    res = s_in * s_in + 1im * (s_out * s_out)
    return (res, 2 * res, 3 * res)
end

function _groundtruth_averaging_norm(proc)
    return 1.0 / (number_incoming_particles(proc) + number_outgoing_particles(proc))
end
function _groundtruth_is_in_phasespace(in_ps,out_ps)
    if in_ps[1] == SFourMomentum(zeros(4))
        return false
    end
        return true
end
    
function _groundtruth_phase_space_factor(in_ps, out_ps)
    en_in = getE.(in_ps)
    en_out = getE.(out_ps)
    return 1 / (prod(en_in) * prod(en_out))
end

function _groundtruth_unsafe_probability(proc, in_ps, out_ps)
    mat_el = _groundtruth_matrix_element(in_ps, out_ps)
    mat_el_sq = abs2.(mat_el)
    normalization = _groundtruth_averaging_norm(proc)
    ps_fac = _groundtruth_phase_space_factor(in_ps, out_ps)
    return sum(mat_el_sq) * ps_fac * normalization
end

function _groundtruth_unsafe_diffCS(proc, in_ps, out_ps)
    init_flux = _groundtruth_incident_flux(in_ps)
    return _groundtruth_unsafe_probability(proc, in_ps, out_ps) / (4 * init_flux)
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

struct TestPhasespaceDef <: AbstractPhasespaceDefinition end
struct TestPhasespaceDef_FAIL <: AbstractPhasespaceDefinition end

_any_fail(x...) = true
_any_fail(::TestProcess, ::TestModel) = false
_any_fail(::TestProcess, ::TestModel, ::TestPhasespaceDef, ::TestPhasespaceDef) = false

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
