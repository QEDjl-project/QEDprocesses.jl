module QEDprocesses

# specific compute models
export PerturbativeQED

# specific scattering processes
export Compton, omega_prime
export ComptonRestSystem
export ComptonSphericalLayout

# generic scattering process
export ScatteringProcess, isphysical

using Reexport
using QEDbase
using QEDcore
using StaticArrays
using QuadGK

include("utils.jl")

include("models/models.jl")

include("processes/utils.jl")

# generic qed process
include("processes/generic_process/utility.jl")
include("processes/generic_process/process.jl")
include("processes/generic_process/perturbative/cross_section.jl")

# one photon compton
include("processes/one_photon_compton/process.jl")
include("processes/one_photon_compton/perturbative/kinematics.jl")
include("processes/one_photon_compton/perturbative/cross_section.jl")
include("processes/one_photon_compton/perturbative/total_probability.jl")

end
