using QEDbase

import QEDbase: AbstractFourMomentum

abstract type AbstractCoordinateSystem end
struct SphericalCoordinateSystem <: AbstractCoordinateSystem end

abstract type AbstractFrameOfReference end
struct CenterOfMomentumFrame <: AbstractFrameOfReference end
struct ElectronRestFrame <: AbstractFrameOfReference end

abstract type AbstractPhasespaceDefinition end

"""
    PhasespaceDefinition(coord_sys::AbstractCoordinateSystem, frame::AbstractFrameOfReference)

Convenient type to dispatch on coordiante systems and frames of reference.
"""
struct PhasespaceDefinition{CS<:AbstractCoordinateSystem,F<:AbstractFrameOfReference} <:
       AbstractPhasespaceDefinition
    coord_sys::CS
    frame::F
end

Broadcast.broadcastable(ps_def::AbstractPhasespaceDefinition) = Ref(ps_def)

# abstract type for generic phase spaces
#
# Currently, elements can be either four-momenta, or real numbers,
# i.e. coordinates.
AbstractPhasespaceElement = Union{AbstractFourMomentum,Real}

"""
    ParticleStateful <: AbstractParticle

Representation of a particle with a state. It has four fields:
- `dir::ParticleDirection`: The direction of the particle, `QEDbase.Incoming()` or `QEDbase.Outgoing()`.
- `species::AbstractParticleType`: The species of the particle, `QEDbase.Electron()`, `QEDbase.Positron()` etc.
- `mom::AbstractFourMomentum`: The momentum of the particle.

Overloads for `QEDbase.is_fermion`, `QEDbase.is_boson`, `QEDbase.is_particle`, `QEDbase.is_anti_particle`, `QEDbase.is_incoming`, `QEDbase.is_outgoing`, `QEDbase.mass`, and `QEDbase.charge` are provided, delegating the call to the correct field and thus implementing the `QEDbase.AbstractParticle` interface.

```jldoctest
julia> using QEDbase; using QEDprocesses

julia> ParticleStateful(Incoming(), Electron(), SFourMomentum(1, 0, 0, 0))
ParticleStateful: incoming electron
    momentum: [1.0, 0.0, 0.0, 0.0]

julia> ParticleStateful(Outgoing(), Photon(), SFourMomentum(1, 0, 0, 0))
ParticleStateful: outgoing photon
    momentum: [1.0, 0.0, 0.0, 0.0]
```
"""
struct ParticleStateful{
    DIR<:ParticleDirection,SPECIES<:AbstractParticleType,ELEMENT<:AbstractFourMomentum
} <: AbstractParticle
    dir::DIR
    species::SPECIES
    mom::ELEMENT

    function ParticleStateful(
        dir::DIR, species::SPECIES, mom::ELEMENT
    ) where {
        DIR<:ParticleDirection,SPECIES<:AbstractParticleType,ELEMENT<:AbstractFourMomentum
    }
        return new{DIR,SPECIES,ELEMENT}(dir, species, mom)
    end
end

"""
    PhaseSpacePoint

Representation of a point in the phase space of a process. Contains the process ([`AbstractProcessDefinition`](@ref)), the model ([`AbstractModelDefinition`](@ref)), the phase space definition ([`AbstractPhasespaceDefinition`]), and stateful incoming and outgoing particles ([`ParticleStateful`](@ref)).

The legality of the combination of the given process and the incoming and outgoing particles is checked on construction. If the numbers of particles mismatch, the types of particles mismatch (note that order is important), or incoming particles have an `Outgoing` direction, an error is thrown.

```jldoctest
julia> using QEDprocesses; using QEDbase

julia> PhaseSpacePoint(
            Compton(), 
            PerturbativeQED(), 
            PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame()), 
            (
                ParticleStateful(Incoming(), Electron(), SFourMomentum(1, 0, 0, 0)), 
                ParticleStateful(Incoming(), Photon(), SFourMomentum(1, 0, 0, 0))
            ), 
            (
                ParticleStateful(Outgoing(), Electron(), SFourMomentum(1, 0, 0, 0)), 
                ParticleStateful(Outgoing(), Photon(), SFourMomentum(1, 0, 0, 0))
            )
        )
PhaseSpacePoint:
    process: one-photon Compton scattering
    model: perturbative QED
    phasespace definition: spherical coordinates in electron rest frame
    incoming particles:
     -> incoming electron: [1.0, 0.0, 0.0, 0.0]
     -> incoming photon: [1.0, 0.0, 0.0, 0.0]
    outgoing particles:
     -> outgoing electron: [1.0, 0.0, 0.0, 0.0]
     -> outgoing photon: [1.0, 0.0, 0.0, 0.0]
```
"""
struct PhaseSpacePoint{
    PROC<:AbstractProcessDefinition,
    MODEL<:AbstractModelDefinition,
    PSDEF<:AbstractPhasespaceDefinition,
    IN_PARTICLES<:Tuple{Vararg{ParticleStateful}},
    OUT_PARTICLES<:Tuple{Vararg{ParticleStateful}},
}
    proc::PROC
    model::MODEL
    ps_def::PSDEF

    in_particles::IN_PARTICLES
    out_particles::OUT_PARTICLES

    function PhaseSpacePoint(
        proc::PROC, model::MODEL, ps_def::PSDEF, in_p::IN_PARTICLES, out_p::OUT_PARTICLES
    ) where {
        PROC<:AbstractProcessDefinition,
        MODEL<:AbstractModelDefinition,
        PSDEF<:AbstractPhasespaceDefinition,
        IN_PARTICLES<:Tuple{Vararg{ParticleStateful}},
        OUT_PARTICLES<:Tuple{Vararg{ParticleStateful}},
    }
        _check_psp(incoming_particles(proc), outgoing_particles(proc), in_p, out_p)

        return new{PROC,MODEL,PSDEF,IN_PARTICLES,OUT_PARTICLES}(
            proc, model, ps_def, in_p, out_p
        )
    end

    """
        PhaseSpacePoint(
            proc::AbstractProcessDefinition, 
            model::AbstractModelDefinition, 
            ps_def::AbstractPhasespaceDefinition, 
            in_ps::AbstractVector{ELEMENT}, 
            out_ps::AbstractVector{ELEMENT},
        ) where {ELEMENT<:QEDbase.AbstractFourMomentum}

    Construct the phase space point from given momenta of incoming and outgoing particles regarding a given process.
    """
    function PhaseSpacePoint(
        proc::PROC,
        model::MODEL,
        ps_def::PSDEF,
        in_ps::AbstractVector{ELEMENT},
        out_ps::AbstractVector{ELEMENT},
    ) where {
        PROC<:AbstractProcessDefinition,
        MODEL<:AbstractModelDefinition,
        PSDEF<:AbstractPhasespaceDefinition,
        ELEMENT<:AbstractFourMomentum,
    }
        length(in_ps) == number_incoming_particles(proc) || throw(
            InvalidInputError(
                "expected $(number_incoming_particles(proc)) incoming particles for the process but got $(length(in_ps))",
            ),
        )
        length(out_ps) == number_outgoing_particles(proc) || throw(
            InvalidInputError(
                "expected $(number_outgoing_particles(proc)) outgoing particles for the process but got $(length(out_ps))",
            ),
        )

        in_particles = Union{
            Tuple{_assemble_tuple_type(incoming_particles(proc), Incoming(), ELEMENT)...}
        }(
            ParticleStateful(Incoming(), particle, mom) for
            (particle, mom) in zip(incoming_particles(proc), in_ps)
        )

        out_particles = Union{
            Tuple{_assemble_tuple_type(outgoing_particles(proc), Outgoing(), ELEMENT)...}
        }(
            ParticleStateful(Outgoing(), particle, mom) for
            (particle, mom) in zip(outgoing_particles(proc), out_ps)
        )

        # no need to check anything since it is constructed correctly above

        return new{PROC,MODEL,PSDEF,typeof(in_particles),typeof(out_particles)}(
            proc, model, ps_def, in_particles, out_particles
        )
    end

    """
        PhaseSpacePoint(
            proc::AbstractProcessDefinition,
            model::AbstractModelDefinition,
            ps_def::AbstractPhasespaceDefinition,
            in_ps::AbstractVector{ELEMENT},
            out_ps::Tuple{},
        ) where {ELEMENT<:AbstractFourMomentum}

    Construct a PhaseSpacePoint with only input particles. The result will be `<: IncomingPhaseSpacePoint` but **not** `<: OutgoingPhaseSpacePoint`. Call this by simply passing an empty `Tuple` as the `out_phasespace`.
    """
    function PhaseSpacePoint(
        proc::PROC, model::MODEL, ps_def::PSDEF, in_ps::AbstractVector{ELEMENT}, ::Tuple{}
    ) where {
        PROC<:AbstractProcessDefinition,
        MODEL<:AbstractModelDefinition,
        PSDEF<:AbstractPhasespaceDefinition,
        ELEMENT<:AbstractFourMomentum,
    }
        length(in_ps) == number_incoming_particles(proc) || throw(
            InvalidInputError(
                "expected $(number_incoming_particles(proc)) incoming particles for the process but got $(length(in_ps))",
            ),
        )
        in_particles = Union{
            Tuple{_assemble_tuple_type(incoming_particles(proc), Incoming(), ELEMENT)...}
        }(
            ParticleStateful(Incoming(), particle, mom) for
            (particle, mom) in zip(incoming_particles(proc), in_ps)
        )

        # no need to check anything since it is constructed correctly above

        return new{PROC,MODEL,PSDEF,typeof(in_particles),Tuple{}}(
            proc, model, ps_def, in_particles, ()
        )
    end

    """
        PhaseSpacePoint(
            proc::AbstractProcessDefinition,
            model::AbstractModelDefinition,
            ps_def::AbstractPhasespaceDefinition,
            in_ps::Tuple{},
            out_ps::AbstractVector{ELEMENT},
        ) where {ELEMENT<:AbstractFourMomentum}

    Construct a PhaseSpacePoint with only output particles. The result will be `<: OutgoingPhaseSpacePoint` but **not** `<: IncomingPhaseSpacePoint`. Call this by simply passing an empty `Tuple` as the `in_phasespace`.
    """
    function PhaseSpacePoint(
        proc::PROC, model::MODEL, ps_def::PSDEF, ::Tuple{}, out_ps::AbstractVector{ELEMENT}
    ) where {
        PROC<:AbstractProcessDefinition,
        MODEL<:AbstractModelDefinition,
        PSDEF<:AbstractPhasespaceDefinition,
        ELEMENT<:AbstractFourMomentum,
    }
        length(out_ps) == number_outgoing_particles(proc) || throw(
            InvalidInputError(
                "expected $(number_outgoing_particles(proc)) outgoing particles for the process but got $(length(out_ps))",
            ),
        )
        out_particles = Union{
            Tuple{_assemble_tuple_type(outgoing_particles(proc), Outgoing(), ELEMENT)...}
        }(
            ParticleStateful(Outgoing(), particle, mom) for
            (particle, mom) in zip(outgoing_particles(proc), out_ps)
        )

        # no need to check anything since it is constructed correctly above

        return new{PROC,MODEL,PSDEF,Tuple{},typeof(out_particles)}(
            proc, model, ps_def, (), out_particles
        )
    end
end

IncomingPhaseSpacePoint{P,M,D,IN,OUT} = PhaseSpacePoint{
    P,M,D,IN,OUT
} where {IN<:Tuple{ParticleStateful,Vararg},OUT<:Tuple{Vararg}}
OutgoingPhaseSpacePoint{P,M,D} = PhaseSpacePoint{
    P,M,D,IN,OUT
} where {IN<:Tuple{Vararg},OUT<:Tuple{ParticleStateful,Vararg}}
