###############
# The paricle interface
#
# In this file, we define the interface of working with particles in a general
# sense. 
# 
# This file is part of `QEDprocesses.jl` which is by itself part of the `QED.jl`
# ecosystem.
#
###############


"""
Abstract base type for every type which might be considered as a `particle` in the context of `QED.jl`. For every (concrete) subtype of `AbstractParticle`, there are two kinds of functions implemented: static functions and property functions. 
The static functions say something, what kind of particle it is (defaults are written in square brackets)

```julia
    is_fermion(::AbstractParticle)::Bool [= false]
    is_boson(::AbstractParticle)::Bool [= false]
    is_particle(::AbstractParticle)::Bool [= true]
    is_anti_particle(::AbstractParticle)::Bool [= false]
``` 
If the output of those functions differ from the defaults for a subtype of `AbstractParticle`, these functions need to be overwritten.
The second type of functions define a soft interface `AbstractParticle`:

```julia
    mass(::AbstractParticle)::Real
    charge(::AbstractParticle)::Real
```
These functions need to be implemented in order to have the subtype of `AbstractParticle` work with the functionalities of `QEDprocesses.jl`.
"""
abstract type AbstractParticle end

"""
    $(TYPEDSIGNATURES)

Interface function for particles. Return `true` if the passed subtype of `AbstractParticle` can be considered as a `fermion` in the sense of particle statistics, and `false` otherwise.

The default implementation of `is_fermion` for every subtype of `AbstractParticle` will always return `false`.
"""
Base.@pure is_fermion(::AbstractParticle) = false

"""
    $(TYPEDSIGNATURES)

Interface function for particles. Return `true` if the passed subtype of `AbstractParticle` can be considered as a `boson` in the sense of particle statistics, and `false` otherwise.
The default implementation of `is_boson` for every subtype of `AbstractParticle` will always return `false`.
"""
Base.@pure is_boson(::AbstractParticle) = false

"""
    $(TYPEDSIGNATURES)
    
Interface function for particles. Return `true` if the passed subtype of `AbstractParticle` can be considered as a *particle* as distinct from anti-particles, and `false` otherwise.
The default implementation of `is_particle` for every subtype of `AbstractParticle` will always return `true`.
"""
Base.@pure is_particle(::AbstractParticle) = true

"""
    $(TYPEDSIGNATURES)

Interface function for particles. Return true if the passed subtype of `AbstractParticle` can be considered as a *anti particle* as distinct from their particle counterpart, and `false` otherwise.
The default implementation of `is_anti_particle` for every subtype of `AbstractParticle` will always return `false`.
"""
Base.@pure is_anti_particle(::AbstractParticle) = false

"""
    $(TYPEDSIGNATURES)

Interface function for particles. Return the rest mass of a particle (in units of the electron mass).

This needs to be implemented for each concrete subtype of `AbstractParticle` and will throw an error otherwise.
"""
function mass(particle::AbstractParticle)::Real
    return error(
        "The function mass($(typeof(particle))) is not implemented. You need to implement it to use the particle interface.",
    )
end

"""
    $(TYPEDSIGNATURES)

Interface function for particles. Return the electric charge of a particle (in units of the elementary electric charge).

This needs to be implemented for each concrete subtype of `AbstractParticle` and will throw an error otherwise.
"""
function charge(::AbstractParticle)::Real
    return error(
        "The function mass($(typeof(particle))) is not implemented. You need to implement it to use the particle interface.",
    )
end
