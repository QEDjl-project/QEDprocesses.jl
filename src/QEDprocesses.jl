module QEDprocesses

# constants
export ALPHA,
    ALPHA_SQUARE, ELEMENTARY_CHARGE, ELEMENTARY_CHARGE_SQUARE, ELECTRONMASS, ONE_OVER_FOURPI

# propagator
export propagator

# specific compute models
export PerturbativeQED

# specific scattering processes
export Compton, omega_prime

# probabilities
export differential_probability, unsafe_differential_probability
export total_probability

# differential cross sections
export differential_cross_section, unsafe_differential_cross_section
export total_cross_section

using QEDbase
using QEDcore
using StaticArrays
using QuadGK

include("constants.jl")
include("utils.jl")

include("cross_section/diff_probability.jl")
include("cross_section/diff_cross_section.jl")
include("cross_section/total_probability.jl")
include("cross_section/total_cross_section.jl")

include("models/models.jl")

# one photon compton
include("processes/one_photon_compton/process.jl")
include("processes/one_photon_compton/perturbative/kinematics.jl")
include("processes/one_photon_compton/perturbative/cross_section.jl")
include("processes/one_photon_compton/perturbative/total_probability.jl")

include("patch_QEDbase.jl")
end
