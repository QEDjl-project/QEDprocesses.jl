#####################
# phase spaces  
#
# This file contains a collection of types and functions to handle phase spaces
# for scattering processes.
#
# TODO: ship to interfaces!
#####################

using StaticArrays

import QEDbase.is_particle
import QEDbase.is_anti_particle
import QEDbase.is_fermion
import QEDbase.is_boson
import QEDbase.is_incoming
import QEDbase.is_outgoing
import QEDbase.mass
import QEDbase.charge

import QEDbase:
    is_particle,
    is_anti_particle,
    is_fermion,
    is_boson,
    is_incoming,
    is_outgoing,
    mass,
    charge,
    AbstractFourMomentum,
    AbstractParticleType,
    AbstractParticle

abstract type AbstractCoordinateSystem end
struct SphericalCoordinateSystem <: AbstractCoordinateSystem end

abstract type AbstractFrameOfReference end
struct CenterOfMomentumFrame <: AbstractFrameOfReference end
struct ElectronRestFrame <: AbstractFrameOfReference end

abstract type AbstractPhasespaceDefinition end

Broadcast.broadcastable(ps_def::AbstractPhasespaceDefinition) = Ref(ps_def)

"""
    PhasespaceDefinition(coord_sys::AbstractCoordinateSystem, frame::AbstractFrameOfReference)

Convenient type to dispatch on coordiante systems and frames of reference.
"""
struct PhasespaceDefinition{CS<:AbstractCoordinateSystem,F<:AbstractFrameOfReference} <:
       AbstractPhasespaceDefinition
    coord_sys::CS
    frame::F
end

function Base.show(io::IO, ::SphericalCoordinateSystem)
    print(io, "spherical coordinates")
    return nothing
end
function Base.show(io::IO, ::CenterOfMomentumFrame)
    print(io, "center-of-momentum frame")
    return nothing
end
function Base.show(io::IO, ::ElectronRestFrame)
    print(io, "electron rest frame")
    return nothing
end
function Base.show(io::IO, m::MIME"text/plain", ps_def::PhasespaceDefinition)
    println(io, "PhasespaceDefinition")
    println(io, "    coordinate system: $(ps_def.coord_sys)")
    println(io, "    frame: $(ps_def.frame)")
    return nothing
end
function Base.show(io::IO, ps_def::PhasespaceDefinition)
    print(io, "$(ps_def.coord_sys) in $(ps_def.frame)")
    return nothing
end

# abstract type for generic phase spaces
#
# Currently, elements can be either four-momenta, or real numbers,
# i.e. coordinates.
AbstractPhasespaceElement = Union{AbstractFourMomentum,Real}

# utility functions 

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
    spin: all spins
    momentum: [1.0, 0.0, 0.0, 0.0]

julia> ParticleStateful(Outgoing(), Photon(), SFourMomentum(1, 0, 0, 0), PolX())
ParticleStateful: outgoing photon
    polarization: x-polarized
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

# particle interface
@inline is_incoming(particle::ParticleStateful) = is_incoming(particle.dir)
@inline is_outgoing(particle::ParticleStateful) = is_outgoing(particle.dir)
@inline is_fermion(particle::ParticleStateful) = is_fermion(particle.species)
@inline is_boson(particle::ParticleStateful) = is_boson(particle.species)
@inline is_particle(particle::ParticleStateful) = is_particle(particle.species)
@inline is_anti_particle(particle::ParticleStateful) = is_anti_particle(particle.species)
@inline mass(particle::ParticleStateful) = mass(particle.species)
@inline charge(particle::ParticleStateful) = charge(particle.species)

# accessors
particle_direction(part::ParticleStateful) = part.dir
particle_species(part::ParticleStateful) = part.species
momentum(part::ParticleStateful) = part.mom

# recursion termination: base case
@inline _assemble_tuple_type(proc::Tuple{}, dir::ParticleDirection, ELTYPE::Type) = ()

# function assembling the correct type information for the tuple of ParticleStatefuls in a phasespace point constructed from momenta
@inline function _assemble_tuple_type(
    proc::Tuple{SPECIES_T,Vararg{AbstractParticleType}}, dir::DIR_T, ELTYPE::Type
) where {SPECIES_T<:AbstractParticleType,DIR_T<:ParticleDirection}
    return (
        ParticleStateful{DIR_T,SPECIES_T,ELTYPE},
        _assemble_tuple_type(proc[2:end], dir, ELTYPE)...,
    )
end

function Base.show(io::IO, particle::ParticleStateful)
    print(
        io, "$(particle.dir) $(particle.species) ($(particle.spin_or_pol)): $(particle.mom)"
    )
    return nothing
end

function Base.show(io::IO, m::MIME"text/plain", particle::ParticleStateful)
    println(io, "ParticleStateful: $(particle.dir) $(particle.species)")
    println(io, "    momentum: $(particle.mom)")
    return nothing
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
            [
                ParticleStateful(Incoming(), Electron(), SFourMomentum(1, 0, 0, 0)), 
                ParticleStateful(Incoming(), Photon(), SFourMomentum(1, 0, 0, 0))
            ], 
            [
                ParticleStateful(Outgoing(), Electron(), SFourMomentum(1, 0, 0, 0)), 
                ParticleStateful(Outgoing(), Photon(), SFourMomentum(1, 0, 0, 0))
            ]
        )
PhaseSpacePoint:
    process: one-photon Compton scattering
    model: perturbative QED
    phasespace definition: spherical coordinates in electron rest frame
    incoming particles:
     -> incoming electron (all spins): [1.0, 0.0, 0.0, 0.0]
     -> incoming photon (all polarizations): [1.0, 0.0, 0.0, 0.0]
    outgoing particles:
     -> outgoing electron (all spins): [1.0, 0.0, 0.0, 0.0]
     -> outgoing photon (all polarizations): [1.0, 0.0, 0.0, 0.0]
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
        in_particles = Tuple{
            _assemble_tuple_type(incoming_particles(proc), Incoming(), ELEMENT)...
        }(
            ParticleStateful(Incoming(), particle, mom) for
            (particle, mom) in zip(incoming_particles(proc), in_ps)
        )

        out_particles = Tuple{
            _assemble_tuple_type(outgoing_particles(proc), Outgoing(), ELEMENT)...
        }(
            ParticleStateful(Outgoing(), particle, mom) for
            (particle, mom) in zip(outgoing_particles(proc), out_ps)
        )

        # no need to check anything since it is constructed correctly above

        return new{PROC,MODEL,PSDEF,typeof(in_particles),typeof(out_particles)}(
            proc, model, ps_def, in_particles, out_particles
        )
    end
end

# recursion termination: success
@inline _recursive_type_check(t::Tuple{}, p::Tuple{}, dir::ParticleDirection) = nothing

# recursion termination: overload for unequal number of particles
@inline function _recursive_type_check(
    ::Tuple{Vararg{ParticleStateful,N}},
    ::Tuple{Vararg{AbstractParticleType,M}},
    dir::ParticleDirection,
) where {N,M}
    return throw(
        InvalidInputError(
            "the number of $(dir) particles in the process $(M) does not match number of particles in the input $(N)",
        ),
    )
end

# recursion termination: overload for invalid types
@inline function _recursive_type_check(
    ::Tuple{ParticleStateful{DIR_IN_T,SPECIES_IN_T},Vararg{ParticleStateful,N}},
    ::Tuple{SPECIES_T,Vararg{AbstractParticleType,M}},
    dir::DIR_T,
) where {
    N,
    M,
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
    OUT_Ts<:Tuple{Vararg{ParticleStateful}},
}
    # in_proc/out_proc contain only species types
    # in_p/out_p contain full ParticleStateful types

    _recursive_type_check(in_p, in_proc, Incoming())
    _recursive_type_check(out_p, out_proc, Outgoing())
    return nothing
end

"""
    Base.getindex(psp::PhaseSpacePoint, dir::Incoming, n::Int)

Overload for the array indexing operator `[]`. Returns the nth incoming particle in this phase space point.
"""
function Base.getindex(psp::PhaseSpacePoint, ::Incoming, n::Int)
    return psp.in_particles[n]
end

"""
    Base.getindex(psp::PhaseSpacePoint, dir::Outgoing, n::Int)

Overload for the array indexing operator `[]`. Returns the nth outgoing particle in this phase space point.
"""
function Base.getindex(psp::PhaseSpacePoint, ::Outgoing, n::Int)
    return psp.out_particles[n]
end

"""
    momentum(psp::PhaseSpacePoint, dir::ParticleDirection, n::Int)

Returns the momentum of the `n`th particle in the given [`PhaseSpacePoint`](@ref) which has direction `dir`. If `n` is outside the valid range for this phase space point, a `BoundsError` is thrown.
"""
function momentum(psp::PhaseSpacePoint, dir::ParticleDirection, n::Int)
    return psp[dir, n].mom
end

function Base.show(io::IO, psp::PhaseSpacePoint)
    print(io, "PhaseSpacePoint of $(psp.proc)")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", psp::PhaseSpacePoint)
    println(io, "PhaseSpacePoint:")
    println(io, "    process: $(psp.proc)")
    println(io, "    model: $(psp.model)")
    println(io, "    phasespace definition: $(psp.ps_def)")
    println(io, "    incoming particles:")
    for p in psp.in_particles
        println(io, "     -> $(p)")
    end
    println(io, "    outgoing particles:")
    for p in psp.out_particles
        println(io, "     -> $(p)")
    end
    return nothing
end
