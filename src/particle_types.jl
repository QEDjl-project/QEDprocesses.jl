###############
# The particle types 
#
# In this file, we define the types of particles used in `QEDprocesses.jl` and
# implement the abstact particle interface accordingly. 
# 
# This file is part of `QEDprocesses.jl` which is by itself part of the `QED.jl`
# ecosystem.
###############

"""
    AbstractParticleType <: AbstractParticle

This is the abstract base type for every species of particles. All functionalities defined on subtypes of `AbstractParticleType` should be static, i.e. known at compile time. 
For adding runtime information, e.g. four-momenta or particle states, to a particle, consider implementing concrete subtype of `AbstractParticle` instead, which may has a type parameter `P<:AbstractParticleType`. See the concrete type `Particle{P,ST,MT}`

Concrete built-in subtypes of `AbstractParticleType` are 

```julia
    Electron
    Positron
    Photon
```

"""
abstract type AbstractParticleType <: AbstractParticle end

"""
Abstract base types for particle species that act like fermions in the sense of particle statistics. 
    
!!! note "particle interface"
    Every concrete subtype of `FermionLike` has `is_fermion(::FermionLike) = true`.
"""
abstract type FermionLike <: AbstractParticleType end

Base.@pure is_fermion(::FermionLike) = true

"""
Abstract base type for fermions as distinct from anti-fermions. 
    
!!! note "particle interface"
    All subtypes of `Fermion` have
    ```julia 
    is_fermion(::Fermion) = true
    is_particle(::Fermion) = true
    is_anti_particle(::Fermion) = false
    ```

"""
abstract type Fermion <: FermionLike end

Base.@pure is_particle(::Fermion) = true

Base.@pure is_anti_particle(::Fermion) = false

"""
Abstract base type for anti-fermions as distinct from its particle counterpart `Fermion`.

!!! note "particle interface"
    All subtypes of `AntiFermion` have 
    ```julia 
    is_fermion(::AntiFermion) = true
    is_particle(::AntiFermion) = false
    is_anti_particle(::AntiFermion) = true
    ```
    
"""
abstract type AntiFermion <: FermionLike end

Base.@pure is_particle(::AntiFermion) = false

Base.@pure is_anti_particle(::AntiFermion) = true

"""
Abstract base type for majorana-fermions, i.e. fermions which are their own anti-particles .

!!! note "particle interface"
    All subtypes of `MajoranaFermion` have 
    ```julia 
    is_fermion(::Fermion) = true
    is_particle(::MajoranaFermion) = true
    is_anti_particle(::MajoranaFermion) = true
    ```
    
"""
abstract type MajoranaFermion <: FermionLike end

Base.@pure is_particle(::MajoranaFermion) = true

Base.@pure is_anti_particle(::MajoranaFermion) = true

"""
Concrete type for the *electrons* as a particle species. Mostly used for dispatch. 

!!! note "particle interface"
    Besides being a subtype of `Fermion`, `Electrons` have

    ```julia
    mass(::Electron) = 1.0
    charge(::Electron) = -1.0
    ```
"""
struct Electron <: Fermion end
mass(::Electron) = 1.0
charge(::Electron) = -1.0

"""
Concrete type for the *positrons* as a particle species. Mostly used for dispatch. 

!!! note "particle interface"
    Besides being a subtype of `AntiFermion`, `Positron` have

    ```julia
    mass(::Positron) = 1.0
    charge(::Positron) = 1.0
    ```
    
"""
struct Positron <: AntiFermion end
mass(::Positron) = 1.0
charge(::Positron) = 1.0

####
# particle information carts - Bosons
####
"""
Abstract base types for particle species that act like bosons in the sense of particle statistics. 
    
!!! note "particle interface"
    Every concrete subtype of `BosonLike` has `is_boson(::BosonLike) = true`.
"""
abstract type BosonLike <: AbstractParticleType end

is_boson(::BosonLike) = true

"""
Abstract base type for bosons as distinct from its anti-particle counterpart `AntiBoson`. 
    
!!! note "particle interface"
    All subtypes of `Boson` have
    ```julia 
    is_boson(::Boson) = true
    is_particle(::Boson) = true
    is_anti_particle(::Boson) = false
    ```

"""
abstract type Boson <: BosonLike end
Base.@pure is_particle(::Boson) = true
Base.@pure is_anti_particle(::Boson) = false

"""
Abstract base type for anti-bosons as distinct from its particle counterpart `Bosons`. 
    
!!! note "particle interface"
    All subtypes of `AntiBoson` have
    ```julia 
    is_boson(::AntiBoson) = true
    is_particle(::AntiBoson) = true
    is_anti_particle(::AntiBoson) = false
    ```

"""
abstract type AntiBoson <: BosonLike end
Base.@pure is_particle(::AntiBoson) = false
Base.@pure is_anti_particle(::AntiBoson) = true

"""
Abstract base type for majorana-bosons, i.e. bosons which are their own anti-particles.

!!! note "particle interface"
    All subtypes of `MajoranaBoson` have 
    ```julia 
    is_boson(::MajoranaBoson) = true
    is_particle(::MajoranaBoson) = true
    is_anti_particle(::MajoranaBoson) = true
    ```
    
"""
abstract type MajoranaBoson <: BosonLike end
Base.@pure is_particle(::MajoranaBoson) = true
Base.@pure is_anti_particle(::MajoranaBoson) = true

# TODO: is traits for massless and uncharged particles
# 		This makes the function calls more specialized
"""
Concrete type for the *photons* as a particle species. Mostly used for dispatch. 

!!! note "particle interface"
    Besides being a subtype of `MajoranaBoson`, `Photon` has

    ```julia
    mass(::Photon) = 0.0
    charge(::Photon) = 0.0
    ```
    
"""
struct Photon <: MajoranaBoson end
mass(::Photon) = 0.0
charge(::Photon) = 0.0
