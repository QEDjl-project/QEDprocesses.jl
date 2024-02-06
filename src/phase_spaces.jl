#####################
# phase spaces  
#
# This file contains a colection of types and functions to handle phase spaces
# for scattering processes.
#
# TODO: ship to interfaces!
#####################

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
