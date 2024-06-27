
############
# Total cross sections
############

"""
    total_cross_section(in_psp::InPhaseSpacePoint)

Return the total cross section for a given `QEDcore.InPhaseSpacePoint`.
"""
function total_cross_section(in_psp::InPhaseSpacePoint)
    I = 1 / (4 * QEDbase._incident_flux(in_psp))
    return I * QEDbase._total_probability(in_psp)
end
