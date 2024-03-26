#################
# The one-photon Compton process
#
# This file contains the implementation of the abstract process setup for
# Compton
##################

export Compton
export omega_prime

include("process.jl")
include("perturbative/kinematics.jl")
include("perturbative/cross_section.jl")
