###########
# Total probability
###########

"""
    total_probability(in_psp::InPhaseSpacePoint)

Return the total probability of a given [`InPhaseSpacePoint`](@ref).
"""
function total_probability(in_psp::InPhaseSpacePoint)
    return _total_probability(in_psp)
end

"""
    total_probability(in_psps::AbstractVector{InPhaseSpacePoint})

Return the total probability of a given vector of [`InPhaseSpacePoint`](@ref)s.
"""
function total_probability(in_psps::AbstractVector{<:InPhaseSpacePoint})
    return _total_probability.(in_psps)
end
