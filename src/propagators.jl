###############
# Propagators
#
# This file contains implementations for the fermion and boson propagators. 
####

import QEDbase: propagator

function _scalar_propagator(K::QEDbase.AbstractFourMomentum, mass::Real)
    return one(mass) / (K * K - mass^2)
end

function _scalar_propagator(K::QEDbase.AbstractFourMomentum)
    return one(QEDbase.getT(K)) / (K * K)
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
