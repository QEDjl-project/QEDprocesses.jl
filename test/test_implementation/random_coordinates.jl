using QEDcore, QEDprocesses
using Random

"""
Return a tuple of tuples of incoming and outgoing coordinates for a given process, model and ps_def that make up a physical phase space point.
"""
function _rand_coordinates(
    rng::AbstractRNG, ::PROCESS, ::MODEL, ::PSL, FLOAT_T=Float64
) where {PROCESS<:Compton,MODEL<:PerturbativeQED,PSL<:AbstractPhaseSpaceLayout}
    return ((rand(rng, FLOAT_T),), (rand(rng, FLOAT_T), rand(rng, FLOAT_T)))
end

tuple_iaspprox(::Tuple{}, ::Tuple{Vararg}; atol=0.0, rtol=eps()) = false
tuple_iaspprox(::Tuple{Vararg}, ::Tuple{}; atol=0.0, rtol=eps()) = false
tuple_isapprox(::Tuple{}, ::Tuple{}; atol=0.0, rtol=eps()) = true
function tuple_isapprox(
    a::Tuple{<:Number,Vararg}, b::Tuple{<:Number,Vararg}; atol=0.0, rtol=eps()
)
    return isapprox(a[1], b[1]; atol=atol, rtol=rtol) &&
           tuple_isapprox(a[2:end], b[2:end]; atol=atol, rtol=rtol)
end
