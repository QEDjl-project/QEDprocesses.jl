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
struct PhasespaceDefinition{CS<:AbstractCoordinateSystem,F<:AbstractFrameOfReference} <:
       AbstractPhasespaceDefinition
    coord_sys::CS
    frame::F
end

"""
    _generate_momenta(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVector{T},
    ) where {T<:Real}
"""
function _generate_incoming_momenta end

"""
    _generate_momenta(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        out_phase_space_def::AbstractPhasespaceDefinition,
        out_phase_space::AbstractVector{T},
    ) where {T<:Real}
"""
function _generate_outgoing_momenta end


# abstract type for generic phase spaces
#
# Currently, elements can be either four-momenta, or real numbers,
# i.e. coordinates.
AbstractPhasespaceElement = Union{QEDbase.AbstractFourMomentum,Real}
