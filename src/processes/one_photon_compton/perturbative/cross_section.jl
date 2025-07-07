#####
# Perturbative one-photon Compton scattering
# Implementation of the cross section interface
#####

function QEDbase._incident_flux(in_psp::InPhaseSpacePoint{<:Compton, PerturbativeQED})
    return momentum(in_psp, Incoming(), 1) * momentum(in_psp, Incoming(), 2)
end

function QEDbase._matrix_element(psp::PhaseSpacePoint{<:Compton, PerturbativeQED})
    in_ps = momenta(psp, Incoming())
    out_ps = momenta(psp, Outgoing())
    return _pert_compton_matrix_element(psp.proc, in_ps, out_ps)
end

"""
    _averaging_norm(::Type{<:Number}, proc::Compton)

!!! note "Convention"

    We average over the initial spins and pols, and sum over final.
"""
function QEDbase._averaging_norm(::Type{T}, proc::Compton) where {T <: Number}
    return one(T) / incoming_multiplicity(proc)
end

@inline function _all_onshell(psp::PhaseSpacePoint{<:Compton})
    return _is_onshell(incoming_particles(psp.proc)[1], momentum(psp, Incoming(), 1)) &&
        _is_onshell(incoming_particles(psp.proc)[2], momentum(psp, Incoming(), 2)) &&
        _is_onshell(outgoing_particles(psp.proc)[1], momentum(psp, Outgoing(), 1)) &&
        _is_onshell(outgoing_particles(psp.proc)[2], momentum(psp, Outgoing(), 2))
end

@inline function QEDbase._is_in_phasespace(psp::PhaseSpacePoint{<:Compton, PerturbativeQED})
    @inbounds if (
            !isapprox(
                momentum(psp, Incoming(), 1) + momentum(psp, Incoming(), 2),
                momentum(psp, Outgoing(), 1) + momentum(psp, Outgoing(), 2);
                rtol = sqrt(eps(momentum_eltype(psp)))
            )
        )
        return false
    end
    return _all_onshell(psp)
end

@inline function QEDbase._phase_space_factor(
        psp::PhaseSpacePoint{<:Compton, PerturbativeQED}
    )
    in_ps = momenta(psp, Incoming())
    out_ps = momenta(psp, Outgoing())
    return _pert_compton_ps_fac(psp.psl, in_ps[2], out_ps[2])
end

#######
# Matrix elements
#######

@inline function _pert_compton_matrix_element(
        proc::Compton, in_ps::NTuple{N, T}, out_ps::NTuple{M, T}
    ) where {N, M, T <: AbstractFourMomentum}
    in_electron_mom = in_ps[1]
    in_photon_mom = in_ps[2]
    out_electron_mom = out_ps[1]
    out_photon_mom = out_ps[2]

    in_electron_state = base_state(Electron(), Incoming(), in_electron_mom, proc.in_spin)
    in_photon_state = base_state(Photon(), Incoming(), in_photon_mom, proc.in_pol)

    out_electron_state = base_state(Electron(), Outgoing(), out_electron_mom, proc.out_spin)
    out_photon_state = base_state(Photon(), Outgoing(), out_photon_mom, proc.out_pol)

    return _pert_compton_matrix_element(
        in_electron_mom,
        in_electron_state,
        in_photon_mom,
        in_photon_state,
        out_electron_mom,
        out_electron_state,
        out_photon_mom,
        out_photon_state,
    )
end

function _pert_compton_matrix_element(
        in_electron_mom::T,
        in_electron_state,
        in_photon_mom::T,
        in_photon_state,
        out_electron_mom::T,
        out_electron_state,
        out_photon_mom::T,
        out_photon_state,
    ) where {T <: AbstractFourMomentum}
    base_states_comb = Iterators.product(
        QEDbase._as_svec(in_electron_state),
        QEDbase._as_svec(in_photon_state),
        QEDbase._as_svec(out_electron_state),
        QEDbase._as_svec(out_photon_state),
    )

    matrix_elements::NTuple{length(base_states_comb), Complex{eltype(T)}} = (
        (
            _pert_compton_matrix_element_single(
                    in_electron_mom,
                    in_el,
                    in_photon_mom,
                    in_ph,
                    out_electron_mom,
                    out_el,
                    out_photon_mom,
                    out_ph,
                ) for (in_el, in_ph, out_el, out_ph) in base_states_comb
        )...,
    )

    return matrix_elements
end

function _pert_compton_matrix_element_single(
        in_electron_mom::T,
        in_electron_state::BiSpinor,
        in_photon_mom::T,
        in_photon_state::SLorentzVector,
        out_electron_mom::T,
        out_electron_state::AdjointBiSpinor,
        out_photon_mom::T,
        out_photon_state::SLorentzVector,
    ) where {T <: AbstractFourMomentum}
    in_ph_slashed = slashed(in_photon_state)
    out_ph_slashed = slashed(out_photon_state)

    prop1 = QEDcore._fermion_propagator(
        in_photon_mom + in_electron_mom, mass(eltype(T), Electron())
    )
    prop2 = QEDcore._fermion_propagator(
        in_electron_mom - out_photon_mom, mass(eltype(T), Electron())
    )

    # TODO: fermion propagator is not yet in QEDbase
    inner_diagram_1 = (out_ph_slashed * (prop1 * in_ph_slashed))
    inner_diagram_2 = (in_ph_slashed * (prop2 * out_ph_slashed))

    result = out_electron_state * (inner_diagram_1 + inner_diagram_2) * in_electron_state

    # TODO: find (preferably unitful) global provider for physical constants
    # elementary charge
    return ELEMENTARY_CHARGE_SQUARE * result
end

#######
# Phase space factors
#######

function _pert_compton_ps_fac(
        in_psl::ComptonSphericalLayout{<:ComptonRestSystem}, in_photon_mom::T, out_photon_mom::T
    ) where {T <: AbstractFourMomentum}
    omega = getE(in_photon_mom)
    omega_prime = getE(out_photon_mom)
    return omega_prime^2 /
        (16 * convert(eltype(T), pi)^2 * omega * mass(eltype(T), Electron()))
end
