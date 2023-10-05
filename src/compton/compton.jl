###############
# Compton
#
# In this directory, we define the Compton type and 
# implement the AbstractScatteringProcess interface to provide its functionality.
###############

export Compton
export polarization

include("scattering_process.jl")
include("differential_cross_section.jl")
