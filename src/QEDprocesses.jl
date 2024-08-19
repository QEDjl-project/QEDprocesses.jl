module QEDprocesses

# constants
export ALPHA,
    ALPHA_SQUARE, ELEMENTARY_CHARGE, ELEMENTARY_CHARGE_SQUARE, ELECTRONMASS, ONE_OVER_FOURPI

# specific compute models
export PerturbativeQED

# specific scattering processes
export Compton, omega_prime
export GenericQEDProcess, isphysical

using QEDbase
using QEDcore
using StaticArrays
using QuadGK

include("constants.jl")
include("utils.jl")

include("models/models.jl")

# generic qed process
include("processes/generic_process/utility.jl")
include("processes/generic_process/process.jl")
include("processes/generic_process/perturbative/cross_section.jl")

# one photon compton
include("processes/one_photon_compton/process.jl")
include("processes/one_photon_compton/perturbative/kinematics.jl")
include("processes/one_photon_compton/perturbative/cross_section.jl")
include("processes/one_photon_compton/perturbative/total_probability.jl")

include("patch_QEDbase.jl")
end
