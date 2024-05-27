"""
Return a tuple of tuples of incoming and outgoing coordinates for a given process, model and ps_def that make up a physical phase space point.
"""
function _rand_coordinates(
    rng::AbstractRNG, ::PROCESS, ::MODEL, ::PS_DEF
) where {PROCESS<:Compton,MODEL<:PerturbativeQED,PS_DEF<:PhasespaceDefinition}
    return ((rand(rng, Float64),), (rand(rng, Float64), rand(rng, Float64)))
end
