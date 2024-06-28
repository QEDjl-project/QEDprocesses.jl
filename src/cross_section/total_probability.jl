###########
# Total probability
###########

"""
    total_probability(in_psp::InPhaseSpacePoint)

Return the total probability of a given `QEDcore.InPhaseSpacePoint`.
"""
function total_probability(in_psp::InPhaseSpacePoint)
    return QEDbase._total_probability(in_psp)
end
