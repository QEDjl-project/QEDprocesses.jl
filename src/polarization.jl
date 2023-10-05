###############
# Polarization
#
# In this file, we define the Polarization type and 
# implement some subtypes of it for use with processes.
###############

abstract type AbstractPolarization end

struct NoPolarization <: AbstractPolarization end

struct XPolarization <: AbstractPolarization end

struct YPolarization <: AbstractPolarization end
