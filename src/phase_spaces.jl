#####################
# phase spaces  
#
# This file contains a colection of types and functions to handle phase spaces
# for scattering processes.
#####################

abstract type AbstractCoordinateSystem end
struct SphericalCoordinateSystem <: AbstractCoordinateSystem end

abstract type AbstractFrameOfReference end
struct CenterOfMomentumFrame <: AbstractFrameOfReference end
struct ElectronRestFrame <: AbstractFrameOfReference end

abstract type AbstractPhasespaceDefinition end
struct PhasespaceDefinition{CS<:AbstractCoordinateSystem,F<:AbstractFrameOfReference} <: AbstractPhasespaceDefinition
    coord_sys::CS
    frame::F
end 
