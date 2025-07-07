@inline function _is_onshell(::Photon, mom::AbstractFourMomentum{T}) where {T <: Number}
    # photons are massless, so use an atol here
    return isapprox(getMass2(mom), mass(T, Photon())^2; atol = eps(T))
end
@inline function _is_onshell(
        p::P, mom::AbstractFourMomentum{T}
    ) where {P <: AbstractParticle, T <: Number}
    return isapprox(getMass2(mom), mass(T, p)^2; rtol = sqrt(eps(T)))
end

@inline function _all_onshell(
        species::Tuple{P, Vararg}, momenta::Tuple{AbstractFourMomentum{T}, Vararg}
    ) where {P <: AbstractParticle, T <: Number}
    return _is_onshell(species[1], momenta[1]) && _all_onshell(species[2:end], momenta[2:end])
end
@inline _all_onshell(::Tuple{}, ::Tuple{}) = true
