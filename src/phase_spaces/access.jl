import QEDbase:
    is_particle,
    is_anti_particle,
    is_fermion,
    is_boson,
    is_incoming,
    is_outgoing,
    mass,
    charge

# particle interface
@inline is_incoming(particle::ParticleStateful) = is_incoming(particle.dir)
@inline is_outgoing(particle::ParticleStateful) = is_outgoing(particle.dir)
@inline is_fermion(particle::ParticleStateful) = is_fermion(particle.species)
@inline is_boson(particle::ParticleStateful) = is_boson(particle.species)
@inline is_particle(particle::ParticleStateful) = is_particle(particle.species)
@inline is_anti_particle(particle::ParticleStateful) = is_anti_particle(particle.species)
@inline mass(particle::ParticleStateful) = mass(particle.species)
@inline charge(particle::ParticleStateful) = charge(particle.species)

# accessors
particle_direction(part::ParticleStateful) = part.dir
particle_species(part::ParticleStateful) = part.species
momentum(part::ParticleStateful) = part.mom

"""
    momenta(psp::PhaseSpacePoint, ::ParticleDirection)

Return a `Tuple` of all the particles' momenta for the given `ParticleDirection`.
"""
momenta(psp::PhaseSpacePoint, ::Incoming) = momentum.(psp.in_particles)
momenta(psp::PhaseSpacePoint, ::Outgoing) = momentum.(psp.out_particles)

"""
    Base.getindex(psp::PhaseSpacePoint, dir::Incoming, n::Int)

Overload for the array indexing operator `[]`. Returns the nth incoming particle in this phase space point.
"""
function Base.getindex(psp::PhaseSpacePoint, ::Incoming, n::Int)
    return psp.in_particles[n]
end

"""
    Base.getindex(psp::PhaseSpacePoint, dir::Outgoing, n::Int)

Overload for the array indexing operator `[]`. Returns the nth outgoing particle in this phase space point.
"""
function Base.getindex(psp::PhaseSpacePoint, ::Outgoing, n::Int)
    return psp.out_particles[n]
end

"""
    momentum(psp::PhaseSpacePoint, dir::ParticleDirection, n::Int)

Returns the momentum of the `n`th particle in the given [`PhaseSpacePoint`](@ref) which has direction `dir`. If `n` is outside the valid range for this phase space point, a `BoundsError` is thrown.
"""
function momentum(psp::PhaseSpacePoint, dir::ParticleDirection, n::Int)
    return psp[dir, n].mom
end
