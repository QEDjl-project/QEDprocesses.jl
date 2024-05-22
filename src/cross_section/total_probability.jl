###########
# Total probability
###########

"""
    total_probability(in_psp::IncomingPhaseSpacePoint)

Return the total probability of a given [`IncomingPhaseSpacePoint`](@ref).
"""
function total_probability(in_psp::IncomingPhaseSpacePoint)
    return _total_probability(in_psp)
end

"""
    total_probability(in_psps::AbstractVector{IncomingPhaseSpacePoint})

Return the total probability of a given vector of [`IncomingPhaseSpacePoint`](@ref)s.
"""
function total_probability(in_psps::AbstractVector{<:IncomingPhaseSpacePoint})
    return _total_probability.(in_psps)
end
