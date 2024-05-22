
############
# Total cross sections
############

"""
    total_cross_section(in_psp::IncomingPhaseSpacePoint)

Return the total cross section for a given [`IncomingPhaseSpacePoint`](@ref).
"""
function total_cross_section(in_psp::IncomingPhaseSpacePoint)
    I = 1 / (4 * _incident_flux(in_psp))
    return I * _total_probability(in_psp)
end

"""
    total_cross_section(in_psps::AbstractVector{IncomingPhaseSpacePoint})

Return the total cross section for a given vector of [`IncomingPhaseSpacePoint`](@ref).
"""
function _total_cross_section(in_psps::AbstractVector{<:IncomingPhaseSpacePoint})
    return _total_cross_section.(in_psps)
end
