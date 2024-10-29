# recursion success, all particles are compatible with their spin/pol
_assert_spin_pol_particle_compatability(::Tuple{}, ::Tuple{}) = nothing

# recursion base case: check first particle against first spin/pol, then recurse
# note: the length of the tuples is expected to be the same, the constructor ensures this by using NTuples
function _assert_spin_pol_particle_compatability(
    particles::Tuple{AbstractParticleType,Vararg},
    spin_pols::Tuple{AbstractSpinOrPolarization,Vararg},
)
    if is_fermion(particles[1]) && !(spin_pols[1] isa AbstractSpin)
        throw(
            InvalidInputError(
                "particle \"$(particles[1])\" is a fermion and should have a spin, but has \"$(spin_pols[1])\"",
            ),
        )
    end
    if is_boson(particles[1]) && !(spin_pols[1] isa AbstractPolarization)
        throw(
            InvalidInputError(
                "particle \"$(particles[1])\" is a boson and should have a polarization, but has \"$(spin_pols[1])\"",
            ),
        )
    end
    return _assert_spin_pol_particle_compatability(particles[2:end], spin_pols[2:end])
end

# this should move to QEDbase as part of the interface, see https://github.com/QEDjl-project/QEDbase.jl/issues/114
_particle_to_letter(::Electron) = "e"
_particle_to_letter(::Positron) = "p"
_particle_to_letter(::Photon) = "k"
