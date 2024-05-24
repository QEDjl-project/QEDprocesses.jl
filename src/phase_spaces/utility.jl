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
    throw(InvalidInputError("expected $(M) $(dir) particles for the process but got $(N)"))
    return nothing
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
    throw(
        InvalidInputError(
            "expected $(dir) $(SPECIES_T()) but got $(DIR_IN_T()) $(SPECIES_IN_T())"
        ),
    )
    return nothing
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
    # specific overload for InPhaseSpacePoint
    _recursive_type_check(in_p, in_proc, Incoming())

    return typeof(in_p[1].mom)
end

@inline function _check_psp(
    in_proc::P_IN_Ts, out_proc::P_OUT_Ts, in_p::IN_Ts, out_p::OUT_Ts
) where {
    P_IN_Ts<:Tuple{Vararg{AbstractParticleType}},
    P_OUT_Ts<:Tuple{Vararg{AbstractParticleType}},
    IN_Ts<:Tuple{},
    OUT_Ts<:Tuple{Vararg{ParticleStateful}},
}
    # specific overload for OutPhaseSpacePoint
    _recursive_type_check(out_p, out_proc, Outgoing())

    return typeof(out_p[1].mom)
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

    return typeof(out_p[1].mom)
end

"""
    _momentum_type(psp::PhaseSpacePoint)
    _momentum_type(type::Type{PhaseSpacePoint})

Returns the element type of the [`PhaseSpacePoint`](@ref) object or type, e.g. `SFourMomentum`.

```jldoctest
julia> using QEDprocesses; using QEDbase;

julia> psp = PhaseSpacePoint(Compton(), PerturbativeQED(), PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame()), Tuple(rand(SFourMomentum) for _ in 1:2), Tuple(rand(SFourMomentum) for _ in 1:2));

julia> QEDprocesses._momentum_type(psp)
SFourMomentum

julia> QEDprocesses._momentum_type(typeof(psp))
SFourMomentum
```
"""
@inline function _momentum_type(
    ::Type{T}
) where {P,M,D,I,O,E,T<:PhaseSpacePoint{P,M,D,I,O,E}}
    return E
end

@inline _momentum_type(::T) where {T<:PhaseSpacePoint} = _momentum_type(T)

# convenience function building a type stable tuple of ParticleStatefuls from the given process, momenta, and direction
function _build_particle_statefuls(
    proc::AbstractProcessDefinition, moms::NTuple{N,ELEMENT}, dir::ParticleDirection
) where {N,ELEMENT<:QEDbase.AbstractFourMomentum}
    N == number_particles(proc, dir) || throw(
        InvalidInputError(
            "expected $(number_particles(proc, dir)) $(dir) particles for the process but got $(N)",
        ),
    )
    res = Tuple{_assemble_tuple_type(particles(proc, dir), dir, ELEMENT)...}(
        ParticleStateful(dir, particle, mom) for
        (particle, mom) in zip(particles(proc, dir), moms)
    )

    return res
end
