#################
# The one-photon Compton process
#
# This file contains the implementation of the abstract process setup for
# Compton
##################

export Compton

include("process.jl")
include("perturbative/kinematics.jl")
include("perturbative/cross_section.jl")
