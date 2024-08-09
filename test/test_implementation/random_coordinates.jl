"""
Return a tuple of tuples of incoming and outgoing coordinates for a given process, model and ps_def that make up a physical phase space point.
"""
function _rand_coordinates(
    rng::AbstractRNG, ::PROCESS, ::MODEL, ::PS_DEF
) where {PROCESS<:Compton,MODEL<:PerturbativeQED,PS_DEF<:PhasespaceDefinition}
    return ((rand(rng, Float64),), (rand(rng, Float64), rand(rng, Float64)))
end

tuple_isapprox(::Tuple{}, ::Tuple{}; atol=0.0, rtol=eps()) = true
function tuple_isapprox(
    a::Tuple{<:Number,Vararg}, b::Tuple{<:Number,Vararg}; atol=0.0, rtol=eps()
)
    return isapprox(a[1], b[1]; atol=atol, rtol=rtol) &&
           tuple_isapprox(a[2:end], b[2:end]; atol=atol, rtol=rtol)
end
