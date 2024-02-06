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
