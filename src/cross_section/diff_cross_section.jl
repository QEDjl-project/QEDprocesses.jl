########################
# differential and total cross sections.
#
# This file contains default implementations for differential 
# cross sections based on the scattering process interface
########################

function _incident_flux(psp::PhaseSpacePoint)
    return _incident_flux(psp.proc, psp.model, momentum.(psp.in_particles))
end

"""
    unsafe_differential_cross_section(phase_space_point::PhaseSpacePoint)

Return the differential cross section evaluated on a phase space point without checking if the given phase space is physical.
"""
function unsafe_differential_cross_section(phase_space_point::PhaseSpacePoint)
    I = 1 / (4 * _incident_flux(phase_space_point))

    return I * unsafe_differential_probability(phase_space_point)
end

"""
    differential_cross_section(phase_space_point::PhaseSpacePoint)

If the given phase spaces are physical, return differential cross section evaluated on a phase space point. Zero otherwise.
"""
function differential_cross_section(phase_space_point::PhaseSpacePoint)
    if !_is_in_phasespace(phase_space_point)
        return zero(eltype(momentum(phase_space_point, Incoming(), 1)))
    end

    return unsafe_differential_cross_section(phase_space_point)
end
