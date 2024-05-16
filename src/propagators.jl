
###############
# Propagators
#
# This file contains implementations for the fermion and boson propagators. 
####

"""

    propagator(particle::AbstractParticleType, mom::QEDbase.AbstractFourMomentum, [mass::Real])

Return the propagator of a particle for a given four-momentum. If `mass` is passed, the respective propagator for massive particles is used, if not, it is assumed the particle passed in is massless.

!!! note "Convention"
    
    There are two types of implementations for propagators given in `QEDProcesses`: 
    For a `BosonLike` particle with four-momentum ``k`` and mass ``m``, the propagator is given as 

    ```math
    D(k) = \\frac{1}{k^2 - m^2}.
    ```

    For a `FermionLike` particle with four-momentum ``p`` and mass ``m``, the propagator is given as

    ```math
    S(p) = \\frac{\\gamma^\\mu p_\\mu + mass}{p^2 - m^2}.
    ```

!!! warning
    
    This function does not throw when the given particle is off-shell. If an off-shell particle is passed, the function `propagator` returns `Inf`.

"""
function propagator end

function _scalar_propagator(K::QEDbase.AbstractFourMomentum, mass::Real)
    return one(mass) / (K * K - mass^2)
end

function _scalar_propagator(K::QEDbase.AbstractFourMomentum)
    return one(getT(K)) / (K * K)
end

function _fermion_propagator(P::QEDbase.AbstractFourMomentum, mass::Real)
    return (slashed(P) + mass * one(DiracMatrix)) * _scalar_propagator(P, mass)
end

function _fermion_propagator(P::QEDbase.AbstractFourMomentum)
    return (slashed(P)) * _scalar_propagator(P)
end

function propagator(particle_type::BosonLike, K::QEDbase.AbstractFourMomentum)
    return _scalar_propagator(K, mass(particle_type))
end

function propagator(particle_type::Photon, K::QEDbase.AbstractFourMomentum)
    return _scalar_propagator(K)
end

function propagator(particle_type::FermionLike, P::QEDbase.AbstractFourMomentum)
    return _fermion_propagator(P, mass(particle_type))
end
