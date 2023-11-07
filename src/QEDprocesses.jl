module QEDprocesses

# Abstract model interface
export AbstractModelDefinition, fundamental_interaction_type

# Abstract process interface
export AbstractProcessDefinition, incoming_particles, outgoing_particles
export in_phasespace_dimension, out_phasespace_dimension
export number_incoming_particles, number_outgoing_particles
export differential_cross_section, total_cross_section

# Abstract setup interface
export AbstractComputationSetup, InvalidInputError, compute
export AbstractProcessSetup, scattering_process, physical_model

# propagator
export propagator

using DocStringExtensions
using QEDbase

include("utils.jl")
include("interfaces/model_interface.jl")
include("interfaces/process_interface.jl")
include("interfaces/setup_interface.jl")
include("propagators.jl")
end
