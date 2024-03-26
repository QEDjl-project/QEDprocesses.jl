#####
# Perturbative one-photon Compton scattering
# Implementation of the cross section interface
#####

function _incident_flux(
    proc::Compton, model::PerturbativeQED, in_ps::AbstractVector{T}
) where {T<:QEDbase.AbstractFourMomentum}
    return prod(in_ps)
end

function _matrix_element(
    proc::Compton,
    model::PerturbativeQED,
    in_ps::AbstractVector{T},
    out_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return _pert_compton_matrix_element(proc, in_ps, out_ps)
end

"""

!!! note "Convention"

    We average over the initial spins and pols, and sum over final.
"""
function _averaging_norm(proc::Compton)
    normalizations = number_of_spin_pol.(_in_spin_and_pol(proc))
    return inv(prod(normalizations))
end

function _all_onshell(
    proc::Compton, in_ps::AbstractVector{T}, out_ps::AbstractVector{T}
) where {T<:QEDbase.AbstractFourMomentum}
    sq_in_moms = getMass2.(in_ps)
    sq_out_moms = getMass2.(out_ps)
    sq_in_masses = mass.(incoming_particles(proc)) .^ 2
    sq_out_masses = mass.(outgoing_particles(proc)) .^ 2
    return isapprox(sq_in_moms, SVector(sq_in_masses)) &&
           isapprox(sq_out_moms, SVector(sq_out_masses))
end

function _is_in_phasespace(
    proc::Compton,
    model::PerturbativeQED,
    in_ps_def::AbstractPhasespaceDefinition,
    in_ps::AbstractVector{T},
    out_ps_def::AbstractPhasespaceDefinition,
    out_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return (!isapprox(sum(in_ps), sum(out_ps))) ? false : _all_onshell(proc, in_ps, out_ps)
end

@inline function _phase_space_factor(
    proc::Compton,
    model::PerturbativeQED,
    in_ps_def::AbstractPhasespaceDefinition,
    in_ps::AbstractVector{T},
    out_ps_def::AbstractPhasespaceDefinition,
    out_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    return _pert_compton_ps_fac(in_ps_def, in_ps[2], out_ps_def, out_ps[2])
end

#######
# Matrix elements
#######

@inline function _pert_compton_matrix_element(
    proc::Compton, in_ps::AbstractVector{T}, out_ps::AbstractVector{T}
) where {T<:QEDbase.AbstractFourMomentum}
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
) where {T<:QEDbase.AbstractFourMomentum}
    base_states_comb = Iterators.product(
        QEDbase._as_svec(in_electron_state),
        QEDbase._as_svec(in_photon_state),
        QEDbase._as_svec(out_electron_state),
        QEDbase._as_svec(out_photon_state),
    )

    matrix_elements = Vector{ComplexF64}()
    sizehint!(matrix_elements, length(base_states_comb))
    for (in_el, in_ph, out_el, out_ph) in base_states_comb
        push!(
            matrix_elements,
            _pert_compton_matrix_element_single(
                in_electron_mom,
                in_el,
                in_photon_mom,
                in_ph,
                out_electron_mom,
                out_el,
                out_photon_mom,
                out_ph,
            ),
        )
    end

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
) where {T<:QEDbase.AbstractFourMomentum}
    in_ph_slashed = slashed(in_photon_state)
    out_ph_slashed = slashed(out_photon_state)

    prop1 = _fermion_propagator(in_photon_mom + in_electron_mom, mass(Electron()))
    prop2 = _fermion_propagator(in_electron_mom - out_photon_mom, mass(Electron()))

    # TODO: fermion propagator is not yet in QEDbase
    diagram_1 =
        out_electron_state *
        (out_ph_slashed * (prop1 * (in_ph_slashed * in_electron_state)))
    diagram_2 =
        out_electron_state *
        (in_ph_slashed * (prop2 * (out_ph_slashed * in_electron_state)))

    result = diagram_1 + diagram_2

    # TODO: find (preferably unitful) global provider for physical constants
    # elementary charge
    return ELEMENTARY_CHARGE_SQUARE * result
end

#######
# Phase space factors
#######

function _pert_compton_ps_fac(
    in_ps_def::PhasespaceDefinition{inCS,ElectronRestFrame},
    in_photon_mom,
    out_ps_def::PhasespaceDefinition{SphericalCoordinateSystem},
    out_photon_mom,
) where {inCS}
    # TODO
    omega = getE(in_photon_mom)
    omega_prime = getE(out_photon_mom)
    return omega_prime^2 / (16 * pi^2 * omega * mass(Electron()))
end
