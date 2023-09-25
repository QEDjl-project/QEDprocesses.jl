module QEDprocesses

export AbstractParticle
export is_fermion, is_boson, is_particle, is_anti_particle
export mass, charge

# Abstract model interface
export AbstractModelDefinition, fundamental_interaction_type

# Abstract process interface
export AbstractProcessDefinition, incoming_particles, outgoing_particles
export initial_phasespace_dimension, final_phasespace_dimension
export number_incoming_particles, number_outgoing_particles
export differential_cross_section, total_cross_section

# particle types
export AbstractParticleType
export FermionLike, Fermion, AntiFermion, MajoranaFermion
export BosonLike, Boson, AntiBoson, MajoranaBoson
export Electron, Positron, Photon

using DocStringExtensions
using QEDbase

include("utils.jl")
include("interfaces/particle_interface.jl")
include("interfaces/model_interface.jl")
include("interfaces/process_interface.jl")
include("particle_types.jl")
end
