############
# scattering probabilities
#
# This file contains implementations of the scattering probability based on the
# process interface with and without input validation and/or phase space
# constraint.
############

# convenience function
# can be overloaded if an analytical version is known
function _matrix_element_square(psp::PhaseSpacePoint)
    mat_el = _matrix_element(psp)
    return abs2.(mat_el)
end

"""
    unsafe_differential_probability(phase_space_point::PhaseSpacePoint)

Return differential probability evaluated on a phase space point without checking if the given phase space(s) are physical.
"""
function unsafe_differential_probability(psp::PhaseSpacePoint)
    matrix_elements_sq = _matrix_element_square(psp)

    normalization = _averaging_norm(psp.proc)

    ps_fac = _phase_space_factor(psp)

    return normalization * sum(matrix_elements_sq) * ps_fac
end

"""
    differential_probability(phase_space_point::PhaseSpacePoint)

If the given phase spaces are physical, return differential probability evaluated on a phase space point. Zero otherwise.
"""
function differential_probability(phase_space_point::PhaseSpacePoint)
    if !_is_in_phasespace(phase_space_point)
        return zero(eltype(momentum(phase_space_point,Incoming(),1)))
    end

    return unsafe_differential_probability(phase_space_point)
end


