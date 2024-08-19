_assert_particle_type_tuple(::Tuple{}) = nothing
function _assert_particle_type_tuple(t::Tuple{AbstractParticleType,Vararg})
    return _assert_particle_type_tuple(t[2:end])
end
function _assert_particle_type_tuple(t::Any)
    throw(
        InvalidInputError(
            "invalid input, provide a tuple of AbstractParticleTypes to construct a GenericQEDProcess",
        ),
    )
end

_assert_spin_pol_particle_compatability(::Tuple{}, ::Tuple{}) = nothing
function _assert_spin_pol_particle_compatability(::Tuple{}, ::Tuple{Vararg})
    throw(InvalidInputError("more spins/pols than particles given"))
end
function _assert_spin_pol_particle_compatability(::Tuple{Vararg}, ::Tuple{})
    throw(InvalidInputError("more particles than spins/pols given"))
end

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
