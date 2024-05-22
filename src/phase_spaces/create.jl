# PSP constructors from particle statefuls

"""
    InPhaseSpacePoint(
        proc::AbstractProcessDefinition, 
        model::AbstractModelDefinition, 
        ps_def::AbstractPhasespaceDefinition, 
        in_ps::Tuple{ParticleStateful},
    )

    Construct a [`PhaseSpacePoint`](@ref) with only input particles from [`ParticleStateful`](@ref)s. The result will be `<: InPhaseSpacePoint` but **not** `<: OutPhaseSpacePoint`.
"""
function InPhaseSpacePoint(
    proc::PROC, model::MODEL, ps_def::PSDEF, in_ps::IN_PARTICLES
) where {
    PROC<:AbstractProcessDefinition,
    MODEL<:AbstractModelDefinition,
    PSDEF<:AbstractPhasespaceDefinition,
    IN_PARTICLES<:Tuple{Vararg{ParticleStateful}},
}
    return PhaseSpacePoint(proc, model, ps_def, in_ps, ())
end

"""
    OutPhaseSpacePoint(
        proc::AbstractProcessDefinition, 
        model::AbstractModelDefinition, 
        ps_def::AbstractPhasespaceDefinition, 
        out_ps::Tuple{ParticleStateful},
    )

Construct a [`PhaseSpacePoint`](@ref) with only output particles from [`ParticleStateful`](@ref)s. The result will be `<: OutPhaseSpacePoint` but **not** `<: InPhaseSpacePoint`.
"""
function OutPhaseSpacePoint(
    proc::PROC, model::MODEL, ps_def::PSDEF, out_ps::OUT_PARTICLES
) where {
    PROC<:AbstractProcessDefinition,
    MODEL<:AbstractModelDefinition,
    PSDEF<:AbstractPhasespaceDefinition,
    OUT_PARTICLES<:Tuple{Vararg{ParticleStateful}},
}
    return PhaseSpacePoint(proc, model, ps_def, (), out_ps)
end

# PSP constructors from momenta

"""
    PhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        in_ps::NTuple{N,AbstractFourMomentum},
        out_ps::NTuple{M,AbstractFourMomentum},
    )

Construct the phase space point from given momenta of incoming and outgoing particles regarding a given process.
"""
function PhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    in_ps::NTuple{N,ELEMENT},
    out_ps::NTuple{M,ELEMENT},
) where {N,M,ELEMENT<:AbstractFourMomentum}
    in_particles = _build_particle_statefuls(proc, in_ps, Incoming())
    out_particles = _build_particle_statefuls(proc, out_ps, Outgoing())

    # no need to check anything since it is constructed correctly above

    return PhaseSpacePoint(proc, model, ps_def, in_particles, out_particles)
end

"""
    InPhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        in_ps::NTuple{N,AbstractFourMomentum},
    )

Construct a [`PhaseSpacePoint`](@ref) with only input particles from given momenta. The result will be `<: InPhaseSpacePoint` but **not** `<: OutPhaseSpacePoint`.
"""
function InPhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    in_ps::NTuple{N,ELEMENT},
) where {N,ELEMENT<:AbstractFourMomentum}
    in_particles = _build_particle_statefuls(proc, in_ps, Incoming())

    # no need to check anything since it is constructed correctly above

    return PhaseSpacePoint(proc, model, ps_def, in_particles, ())
end

"""
    OutPhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        out_ps::NTuple{N,AbstractFourMomentum},
    )

Construct a [`PhaseSpacePoint`](@ref) with only output particles from given momenta. The result will be `<: OutPhaseSpacePoint` but **not** `<: InPhaseSpacePoint`.
"""
function OutPhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    out_ps::NTuple{N,ELEMENT},
) where {N,ELEMENT<:AbstractFourMomentum}
    out_particles = _build_particle_statefuls(proc, out_ps, Outgoing())
    # no need to check anything since it is constructed correctly above

    return PhaseSpacePoint(proc, model, ps_def, (), out_particles)
end

# PSP constructors from coordinates

"""
    PhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        in_ps::NTuple{N,Real},
        out_ps::NTuple{M,Real},
    )

Construct a [`PhaseSpacePoint`](@ref) from given coordinates by using the [`_generate_momenta`](@ref) interface.
"""
function PhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    in_ps::NTuple{N,Real},
    out_ps::NTuple{M,Real},
) where {N,M}
    in_ps, out_ps = _generate_momenta(proc, model, ps_def, in_ps, out_ps)
    return PhaseSpacePoint(proc, model, ps_def, in_ps, out_ps)
end

"""
    InPhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        in_ps::NTuple{N,Real},
    )

Construct a [`PhaseSpacePoint`](@ref) from given coordinates by using the [`_generate_momenta`](@ref) interface. The result will be `<: InPhaseSpacePoint` but **not** `<: OutPhaseSpacePoint`.

!!! note
    A similar function for [`OutPhaseSpacePoint`](@ref) does not exist from coordinates, only a full [`PhaseSpacePoint`](@ref).
"""
function InPhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    in_ps::NTuple{N,Real},
) where {N}
    in_ps = _generate_incoming_momenta(proc, model, ps_def, in_ps)
    return InPhaseSpacePoint(proc, model, ps_def, in_ps)
end
