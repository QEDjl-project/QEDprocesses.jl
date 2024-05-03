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
- `type::AbstractParticle`: The type of the particle, [`Electron`](@ref)(), [`Positron`](@ref)() etc..
- `spinorpol::AbstractSpinOrPolarization`: The spin or polarization of the particle, [`SpinUp`](@ref)(), [`PolX`](@ref)() etc.. Can only use spins with fermions and polarizations with bosons.
- `dir::ParticleDirection`: The direction of the particle, [`Incoming`](@ref)() or [`Outgoing`](@ref)().

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
end

@inline is_incoming(particle::ParticleStateful) = is_incoming(particle.dir)
@inline is_outgoing(particle::ParticleStateful) = is_outgoing(particle.dir)
@inline is_fermion(particle::ParticleStateful) = is_fermion(particle.type)
@inline is_boson(particle::ParticleStateful) = is_boson(particle.type)
@inline is_particle(particle::ParticleStateful) = is_particle(particle.type)
@inline is_anti_particle(particle::ParticleStateful) = is_anti_particle(particle.type)
