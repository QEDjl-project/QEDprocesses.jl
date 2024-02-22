module QEDprocesses

# constants
export 
    ALPHA,
    ALPHA_SQUARE,
    ELEMENTARY_CHARGE,
    ELEMENTARY_CHARGE_SQUARE,
    ELECTRONMASS,
    ONE_OVER_FOURPI

# Abstract model interface
export AbstractModelDefinition, fundamental_interaction_type

# Abstract process interface
export AbstractProcessDefinition, incoming_particles, outgoing_particles
export number_incoming_particles, number_outgoing_particles

# probabilities
export differential_probability, unsafe_differential_probability
export total_probability

# probabilities
export differential_probability, unsafe_differential_probability
export total_probability

# probabilities
export differential_probability, unsafe_differential_probability
export total_probability

# differential cross section
export differential_cross_section, unsafe_differential_cross_section
export total_cross_section

# Abstract setup interface
export AbstractComputationSetup, InvalidInputError, compute
export AbstractProcessSetup, scattering_process, physical_model

# propagator
export propagator

# phase space 
export AbstractCoordinateSystem, SphericalCoordinateSystem
export AbstractFrameOfReference, CenterOfMomentumFrame, ElectronRestFrame
export AbstractPhasespaceDefinition, PhasespaceDefinition

using QEDbase

include("utils.jl")
include("interfaces/model_interface.jl")
include("interfaces/process_interface.jl")
include("interfaces/setup_interface.jl")
include("phase_spaces.jl")
include("momentum_generation.jl")
include("propagators.jl")
include("probabilities.jl")
include("cross_sections.jl")

include("processes/one_photon_compton/one_photon_compton.jl")

include("utils.jl")

include("patch_QEDbase.jl")
end
