_assert_particle_type_tuple(::Tuple{}) = nothing
function _assert_particle_type_tuple(t::Tuple{AbstractParticleType,Vararg})
    return _assert_particle_type_tuple(t[2:end])
end
function _assert_particle_type_tuple(t::Any)
    throw(
        InvalidInputError(
            "invalid input, provide a tuple of AbstractParticleTypes to construct a QEDProcess",
        ),
    )
end
