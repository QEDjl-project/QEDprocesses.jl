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
        in_momenta::NTuple{N,QEDbase.AbstractFourMomentum},
        out_momenta::NTuple{M,QEDbase.AbstractFourMomentum},
    )

Construct the phase space point from given momenta of incoming and outgoing particles regarding a given process.
"""
function PhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    in_momenta::NTuple{N,ELEMENT},
    out_momenta::NTuple{M,ELEMENT},
) where {N,M,ELEMENT<:QEDbase.AbstractFourMomentum}
    in_particles = _build_particle_statefuls(proc, in_momenta, Incoming())
    out_particles = _build_particle_statefuls(proc, out_momenta, Outgoing())

    return PhaseSpacePoint(proc, model, ps_def, in_particles, out_particles)
end

"""
    InPhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        in_momenta::NTuple{N,QEDbase.AbstractFourMomentum},
    )

Construct a [`PhaseSpacePoint`](@ref) with only input particles from given momenta. The result will be `<: InPhaseSpacePoint` but **not** `<: OutPhaseSpacePoint`.
"""
function InPhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    in_momenta::NTuple{N,ELEMENT},
) where {N,ELEMENT<:QEDbase.AbstractFourMomentum}
    in_particles = _build_particle_statefuls(proc, in_momenta, Incoming())

    return PhaseSpacePoint(proc, model, ps_def, in_particles, ())
end

"""
    OutPhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        out_momenta::NTuple{N,QEDbase.AbstractFourMomentum},
    )

Construct a [`PhaseSpacePoint`](@ref) with only output particles from given momenta. The result will be `<: OutPhaseSpacePoint` but **not** `<: InPhaseSpacePoint`.
"""
function OutPhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    out_momenta::NTuple{N,ELEMENT},
) where {N,ELEMENT<:QEDbase.AbstractFourMomentum}
    out_particles = _build_particle_statefuls(proc, out_momenta, Outgoing())

    return PhaseSpacePoint(proc, model, ps_def, (), out_particles)
end

# PSP constructors from coordinates

"""
    PhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        in_coords::NTuple{N,Real},
        out_coords::NTuple{M,Real},
    )

Construct a [`PhaseSpacePoint`](@ref) from given coordinates by using the [`_generate_momenta`](@ref) interface.
"""
function PhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    in_coords::NTuple{N,Real},
    out_coords::NTuple{M,Real},
) where {N,M}
    in_ps, out_ps = _generate_momenta(proc, model, ps_def, in_coords, out_coords)
    return PhaseSpacePoint(proc, model, ps_def, in_ps, out_ps)
end

"""
    InPhaseSpacePoint(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        ps_def::AbstractPhasespaceDefinition,
        in_coords::NTuple{N,Real},
    )

Construct a [`PhaseSpacePoint`](@ref) from given coordinates by using the [`_generate_momenta`](@ref) interface. The result will be `<: InPhaseSpacePoint` but **not** `<: OutPhaseSpacePoint`.

!!! note
    A similar function for [`OutPhaseSpacePoint`](@ref) does not exist from coordinates, only a full [`PhaseSpacePoint`](@ref).
"""
function InPhaseSpacePoint(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    ps_def::AbstractPhasespaceDefinition,
    in_coords::NTuple{N,Real},
) where {N}
    in_ps = _generate_incoming_momenta(proc, model, ps_def, in_coords)
    return InPhaseSpacePoint(proc, model, ps_def, in_ps)
end
