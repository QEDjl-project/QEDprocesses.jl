@inline function _is_onshell(::Photon, mom::AbstractFourMomentum{T}) where {T <: Number}
    # photons are massless, so use an atol here
    return isapprox(getMass2(mom), mass(T, Photon())^2; atol = sqrt(eps(T)))
end
@inline function _is_onshell(
        p::P, mom::AbstractFourMomentum{T}
    ) where {P <: AbstractParticleType, T <: Number}
    return isapprox(getMass2(mom), mass(T, p)^2; rtol = sqrt(eps(T)))
end

@inline function _is_onshell(p::AbstractParticleStateful)
    return _is_onshell(particle_species(p), momentum(p))
end

@inline function _all_onshell(
        particles::Tuple{P, Vararg}
    ) where {P <: AbstractParticleStateful}
    return _is_onshell(particles[1]) && _all_onshell(particles[2:end])
end
@inline _all_onshell(::Tuple{}) = true
