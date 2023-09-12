module QEDprocesses

# Write your package code here.
export AbstractParticle
export is_fermion, is_boson, is_particle, is_anti_particle
export mass, charge

using DocStringExtensions

include("interfaces/particle_interface.jl")
end
