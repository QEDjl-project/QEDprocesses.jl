@inline function _check_in_phase_space_dimension(
    proc::AbstractProcessDefinition,
    ::AbstractModelDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:AbstractFourMomentum}
    return size(in_phase_space, 1) == number_incoming_particles(proc) || throw(
        DimensionMismatch(
            "the number of incoming particles <$(number_incoming_particles(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )
end

@inline function _check_out_phase_space_dimension(
    proc::AbstractProcessDefinition,
    ::AbstractModelDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:AbstractFourMomentum}
    return size(out_phase_space, 1) == number_outgoing_particles(proc) || throw(
        DimensionMismatch(
            "the number of outgoing particles <$(number_outgoing_particles(proc))> is inconsistent with input size <$(size(out_phase_space,1))>",
        ),
    )
end

@inline function _check_in_phase_space_dimension(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    return size(in_phase_space, 1) == in_phase_space_dimension(proc, model) || throw(
        DimensionMismatch(
            "the dimension of the in-phase-space <$(in_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )
end

@inline function _check_out_phase_space_dimension(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    return size(out_phase_space, 1) == out_phase_space_dimension(proc, model) || throw(
        DimensionMismatch(
            "the dimension of the out-phase-space <$(out_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(out_phase_space,1))>",
        ),
    )
end

# recursion termination: base case
@inline _assemble_tuple_type(::Tuple{}, ::ParticleDirection, ::Type) = ()

# function assembling the correct type information for the tuple of ParticleStatefuls in a phasespace point constructed from momenta
@inline function _assemble_tuple_type(
    particle_types::Tuple{SPECIES_T,Vararg{AbstractParticleType}}, dir::DIR_T, ELTYPE::Type
) where {SPECIES_T<:AbstractParticleType,DIR_T<:ParticleDirection}
    return (
        ParticleStateful{DIR_T,SPECIES_T,ELTYPE},
        _assemble_tuple_type(particle_types[2:end], dir, ELTYPE)...,
    )
end

# recursion termination: success
@inline _recursive_type_check(::Tuple{}, ::Tuple{}, ::ParticleDirection) = nothing

# recursion termination: overload for unequal number of particles
@inline function _recursive_type_check(
    ::Tuple{Vararg{ParticleStateful,N}},
    ::Tuple{Vararg{AbstractParticleType,M}},
    dir::ParticleDirection,
) where {N,M}
    return throw(
        InvalidInputError("expected $(M) $(dir) particles for the process but got $(N)")
    )
end

# recursion termination: overload for invalid types
@inline function _recursive_type_check(
    ::Tuple{ParticleStateful{DIR_IN_T,SPECIES_IN_T},Vararg{ParticleStateful,N}},
    ::Tuple{SPECIES_T,Vararg{AbstractParticleType,N}},
    dir::DIR_T,
) where {
    N,
    DIR_IN_T<:ParticleDirection,
    DIR_T<:ParticleDirection,
    SPECIES_IN_T<:AbstractParticleType,
    SPECIES_T<:AbstractParticleType,
}
    return throw(
        InvalidInputError(
            "expected $(dir) $(SPECIES_T()) but got $(DIR_IN_T()) $(SPECIES_IN_T())"
        ),
    )
end

@inline function _recursive_type_check(
    t::Tuple{ParticleStateful{DIR_T,SPECIES_T},Vararg{ParticleStateful,N}},
    p::Tuple{SPECIES_T,Vararg{AbstractParticleType,N}},
    dir::DIR_T,
) where {N,DIR_T<:ParticleDirection,SPECIES_T<:AbstractParticleType}
    return _recursive_type_check(t[2:end], p[2:end], dir)
end

@inline function _check_psp(
    in_proc::P_IN_Ts, out_proc::P_OUT_Ts, in_p::IN_Ts, out_p::OUT_Ts
) where {
    P_IN_Ts<:Tuple{Vararg{AbstractParticleType}},
    P_OUT_Ts<:Tuple{Vararg{AbstractParticleType}},
    IN_Ts<:Tuple{Vararg{ParticleStateful}},
    OUT_Ts<:Tuple{},
}
    # specific overload for IncomingPhaseSpacePoint
    _recursive_type_check(in_p, in_proc, Incoming())
    return nothing
end

@inline function _check_psp(
    in_proc::P_IN_Ts, out_proc::P_OUT_Ts, in_p::IN_Ts, out_p::OUT_Ts
) where {
    P_IN_Ts<:Tuple{Vararg{AbstractParticleType}},
    P_OUT_Ts<:Tuple{Vararg{AbstractParticleType}},
    IN_Ts<:Tuple{},
    OUT_Ts<:Tuple{Vararg{ParticleStateful}},
}
    # specific overload for OutgoingPhaseSpacePoint
    _recursive_type_check(out_p, out_proc, Outgoing())
    return nothing
end

@inline function _check_psp(
    in_proc::P_IN_Ts, out_proc::P_OUT_Ts, in_p::IN_Ts, out_p::OUT_Ts
) where {
    P_IN_Ts<:Tuple{Vararg{AbstractParticleType}},
    P_OUT_Ts<:Tuple{Vararg{AbstractParticleType}},
    IN_Ts<:Tuple{Vararg{ParticleStateful}},
    OUT_Ts<:Tuple{Vararg{ParticleStateful}},
}
    # in_proc/out_proc contain only species types
    # in_p/out_p contain full ParticleStateful types

    _recursive_type_check(in_p, in_proc, Incoming())
    _recursive_type_check(out_p, out_proc, Outgoing())
    return nothing
end
