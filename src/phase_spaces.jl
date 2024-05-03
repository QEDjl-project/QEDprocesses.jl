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

# abstract type for generic phase spaces
#
# Currently, elements can be either four-momenta, or real numbers,
# i.e. coordinates.
AbstractPhasespaceElement = Union{QEDbase.AbstractFourMomentum,Real}

# utility functions 

@inline function _check_in_phase_space_dimension(
    proc::AbstractProcessDefinition,
    ::AbstractModelDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return size(in_phase_space, 1) == number_incoming_particles(proc) || throw(
        DimensionMismatch(
            "The number of incoming particles <$(number_incoming_particles(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )
end

@inline function _check_out_phase_space_dimension(
    proc::AbstractProcessDefinition,
    ::AbstractModelDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return size(out_phase_space, 1) == number_outgoing_particles(proc) || throw(
        DimensionMismatch(
            "The number of outgoing particles <$(number_outgoing_particles(proc))> is inconsistent with input size <$(size(out_phase_space,1))>",
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
            "The dimension of the in-phase-space <$(in_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(in_phase_space,1))>",
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
            "The dimension of the out-phase-space <$(out_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(out_phase_space,1))>",
        ),
    )
end

"""
    ParticleStateful

Representation of a particle with a state. It has four fields:
- `mom::AbstractPhasespaceElement`: The momentum of the particle.
- `type::AbstractParticle`: The type of the particle, `QEDbase.Electron()`, `QEDbase.Positron()` etc..
- `spinorpol::AbstractSpinOrPolarization`: The spin or polarization of the particle, `QEDbase.SpinUp()`, `QEDbase.PolX() etc.. Can only use spins with fermions and polarizations with bosons.
- `dir::ParticleDirection`: The direction of the particle, `QEDbase.Incoming()` or `QEDbase.Outgoing()`.

Overloads for `QEDbase.is_fermion`, `QEDbase.is_boson`, `QEDbase.is_particle`, `QEDbase.is_anti_particle`, `QEDbase.is_incoming`, and `QEDbase.is_outgoing` are provided, delegating the call to the correct field.

The legality of the combination of `type` and `spinorpol` is checked on construction. If, for example, the construction of an `Electron()` with a polarization is attempted, an [`InvalidInputError`](@ref) is thrown.
"""
struct ParticleStateful{ElType<:AbstractPhasespaceElement}
    mom::ElType
    type::QEDbase.AbstractParticle
    spinorpol::AbstractSpinOrPolarization
    dir::ParticleDirection

    function ParticleStateful(
        mom::ElType,
        type::QEDbase.AbstractParticle,
        spinorpol::AbstractSpinOrPolarization,
        dir::ParticleDirection,
    ) where {ElType<:AbstractPhasespaceElement}
        if is_fermion(type) && !(spinorpol isa AbstractSpin)
            throw(
                InvalidInputError(
                    "Tried to construct a stateful fermion with a $(spinorpol) instead of a spin",
                ),
            )
        elseif is_boson(type) && !(spinorpol isa AbstractPolarization)
            throw(
                InvalidInputError(
                    "Tried to construct a stateful boson with a $(spinorpol) instead of a polarization",
                ),
            )
        end

        return new{ElType}(mom, type, spinorpol, dir)
    end
end

"""
    PhaseSpacePoint

Representation of a point in the phase space of a process. Contains the process ([`AbstractProcessDefinition`](@ref)), the model ([`AbstractModelDefinition`](@ref)), the phase space definition ([`AbstractPhasespaceDefinition`]), and stateful incoming and outgoing particles ([`ParticleStateful`](@ref)).

The legality of the combination of the given process and the incoming and outgoing particles is checked on construction. If the numbers of particles mismatch, the types of particles mismatch (note that order is important), or incoming particles have a `Outgoing` direction, an error is thrown.
"""
struct PhaseSpacePoint{
    PROC<:AbstractProcessDefinition,
    MODEL<:AbstractModelDefinition,
    PSDEF<:AbstractPhasespaceDefinition,
    PhaseSpaceElementType<:AbstractPhasespaceElement,
    N_IN_PARTICLES,
    N_OUT_PARTICLES,
}
    proc::PROC
    model::MODEL
    ps_def::PSDEF

    in_particles::SVector{N_IN_PARTICLES,ParticleStateful{PhaseSpaceElementType}}
    out_particles::SVector{N_OUT_PARTICLES,ParticleStateful{PhaseSpaceElementType}}

    function PhaseSpacePoint(
        proc::PROC, model::MODEL, ps_def::PSDEF, in_p::IN_P, out_p::OUT_P
    ) where {
        PROC<:AbstractProcessDefinition,
        MODEL<:AbstractModelDefinition,
        PSDEF<:AbstractPhasespaceDefinition,
        PhaseSpaceElementType<:AbstractPhasespaceElement,
        IN_P<:AbstractVector{ParticleStateful{PhaseSpaceElementType}},
        OUT_P<:AbstractVector{ParticleStateful{PhaseSpaceElementType}},
    }
        length(incoming_particles(proc)) == length(in_p) || throw(
            InvalidInputError(
                "The number of incoming particles given by the process ($(incoming_particles(proc))) mismatches the number of given stateful incoming particles ($(length(in_p)))",
            ),
        )
        length(outgoing_particles(proc)) == length(out_p) || throw(
            InvalidInputError(
                "The number of outgoing particles given by the process ($(outgoing_particles(proc))) mismatches the number of given stateful outgoing particles ($(length(out_p)))",
            ),
        )

        for (proc_p, p) in zip(incoming_particles(proc), in_p)
            proc_p == p.type || throw(
                InvalidInputError(
                    "Process given particle type ($(proc_p)) does not match stateful particle type ($(p.type))",
                ),
            )
            is_incoming(p) || throw(
                InvalidInputError(
                    "Stateful particle $(p) is given as an incoming particle but is outgoing",
                ),
            )
        end
        for (proc_p, p) in zip(outgoing_particles(proc), out_p)
            proc_p == p.type || throw(
                InvalidInputError(
                    "Process given particle type ($(proc_p)) does not match stateful particle type ($(p.type))",
                ),
            )
            is_outgoing(p) || throw(
                InvalidInputError(
                    "Stateful particle $(p) is given as an outgoing particle but is incoming",
                ),
            )
        end

        return new{PROC,MODEL,PSDEF,PhaseSpaceElementType,length(in_p),length(out_p)}(
            proc, model, ps_def, in_p, out_p
        )
    end
end

@inline is_incoming(particle::ParticleStateful) = is_incoming(particle.dir)
@inline is_outgoing(particle::ParticleStateful) = is_outgoing(particle.dir)
@inline is_fermion(particle::ParticleStateful) = is_fermion(particle.type)
@inline is_boson(particle::ParticleStateful) = is_boson(particle.type)
@inline is_particle(particle::ParticleStateful) = is_particle(particle.type)
@inline is_anti_particle(particle::ParticleStateful) = is_anti_particle(particle.type)
