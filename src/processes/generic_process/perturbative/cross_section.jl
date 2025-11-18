@inline function _all_onshell(psp::PhaseSpacePoint{<:ScatteringProcess})
    return _all_onshell(particles(psp, Incoming())) && _all_onshell(particles(psp, Outgoing()))
end

function _scattering_proc_from_type(::Type{ScatteringProcess{IN_T, OUT_T, IN_SP, OUT_SP}}) where {
        I, O,
        IN_T <: NTuple{I, AbstractParticleType},
        OUT_T <: NTuple{O, AbstractParticleType},
        IN_SP <: NTuple{I, AbstractSpinOrPolarization},
        OUT_SP <: NTuple{O, AbstractSpinOrPolarization},
    }
    _ctor(::Type{T}) where {T <: Tuple} = tuple((Base.issingletontype(Ti) ? Ti() : Ti() for Ti in T.parameters)...)
    return ScatteringProcess(
        _ctor(IN_T),
        _ctor(OUT_T),
        _ctor(IN_SP),
        _ctor(OUT_SP),
    )
end

function QEDbase._matrix_element(psp::PhaseSpacePoint{PROC, PerturbativeQED}) where {PROC <: ScatteringProcess}
    mat_el_func = _mat_el_sq_sum_func(process(psp), typeof(psp))
    return sqrt(mat_el_func(psp))
end

function QEDbase._matrix_element_square_sum(psp::PhaseSpacePoint{PROC, PerturbativeQED}) where {PROC <: ScatteringProcess}
    mat_el_func = _mat_el_sq_sum_func(process(psp), typeof(psp))
    return mat_el_func(psp)
end

function QEDbase._averaging_norm(::Type{T}, proc::ScatteringProcess) where {T <: Number}
    return one(T) / incoming_multiplicity(proc)
end

function QEDbase._is_in_phasespace(psp::PhaseSpacePoint{<:ScatteringProcess, PerturbativeQED})
    if (
            !isapprox(
                sum(momenta(psp, Incoming())),
                sum(momenta(psp, Outgoing()));
                rtol = 100 * eps(momentum_eltype(psp)),
                atol = 100 * eps(momentum_eltype(psp))
            )
        )
        return false
    end
    return _all_onshell(psp)
end

function QEDbase._incident_flux(psp::InPhaseSpacePoint{PROC, PerturbativeQED}) where {PROC <: ScatteringProcess}
    proc = process(psp)
    if length(incoming_particles(proc)) > 2
        throw("_incident_flux is unimplemented for general scattering processes with more than 2 incoming particles")
    end
    if length(incoming_particles(proc)) < 1
        throw("_incident_flux is not defined for scattering processes with less than 2 incoming particles")
    end

    p1_mom = momentum(psp, Incoming(), Val(1))
    p2_mom = momentum(psp, Incoming(), Val(2))

    EL_TYPE = eltype(p1_mom)

    p1_mass = mass(EL_TYPE, incoming_particles(proc)[1])
    p2_mass = mass(EL_TYPE, incoming_particles(proc)[2])

    return QEDcore.sq_diff_sqrt(p1_mom * p2_mom, p1_mass * p2_mass)
end

function QEDbase.unsafe_differential_cross_section!(dest::AbstractVector, in_psps::AbstractVector{PSP}) where {
        MODEL <: PerturbativeQED,
        PROC <: ScatteringProcess,
        PSP <: AbstractPhaseSpacePoint{PROC, MODEL},
    }
    @assert length(in_psps) == length(dest)
    k = _diff_cs_kernel(_scattering_proc_from_type(PROC), PSP)
    return k(get_backend(dest))(dest, in_psps; ndrange = length(dest))
end

function QEDbase.unsafe_differential_probability!(dest::AbstractVector, in_psps::AbstractVector{PSP}) where {
        MODEL <: PerturbativeQED,
        PROC <: ScatteringProcess,
        PSP <: AbstractPhaseSpacePoint{PROC, MODEL},
    }
    @assert length(in_psps) == length(dest)
    k = _diff_prob_kernel(_scattering_proc_from_type(PROC), PSP)
    return k(get_backend(dest))(dest, in_psps; ndrange = length(dest))
end
