module QEDprocesses

export AbstractParticle
export is_fermion, is_boson, is_particle, is_anti_particle
export mass, charge

<<<<<<< HEAD
# particle types
export AbstractParticleType
export FermionLike, Fermion, AntiFermion, MajoranaFermion
export BosonLike, Boson, AntiBoson, MajoranaBoson
export Electron, Positron, Photon

using DocStringExtensions

include("interfaces/particle_interface.jl")
include("particle_types.jl")
end
