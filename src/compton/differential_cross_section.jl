
function _assert_valid_input(
    ::Compton,
    ::PerturbativeQED,
    in_phase_space::AbstractArray{NumericType},
    out_phase_space::AbstractArray{NumericType},
) where {NumericType<:QEDbase.AbstractFourMomentum}
    @assert is_onshell(out_phase_space[1])
    @assert is_onshell(out_phase_space[2], 1.0)
    return nothing
end

function _differential_cross_section(
    process::Compton,
    model::PerturbativeQED,
    in_phase_space::AbstractVector{NumericType},
    out_phase_space::AbstractVector{NumericType},
)::Float64 where {NumericType<:QEDbase.AbstractFourMomentum}
    if (!isapprox(sum(in_phase_space), sum(out_phase_space); rtol=sqrt(eps())))
        return zero(Float64)
    end

    photon_in = in_phase_space[1]
    electron_in = in_phase_space[2]
    photon_out = out_phase_space[1]
    electron_out = out_phase_space[2]

    # get base states of the particles
    photon_in_bstate = base_state(Photon(), Incoming(), photon_in, _spin_or_pol(process, Photon(), Incoming()))
    electron_in_bstate = base_state(Electron(), Incoming(), electron_in, _spin_or_pol(process, Electron(), Incoming()))
    photon_out_bstate = base_state(Photon(), Outgoing(), photon_out, _spin_or_pol(process, Photon(), Outgoing()))
    electron_out_bstate = base_state(Electron(), Outgoing(), electron_out, _spin_or_pol(process, Electron(), Outgoing()))

    # if the particles had AllSpin or AllPol, the base states can be vectors and we need to consider every combination of the base states with each other
    base_states_comb = Iterators.product(photon_in_bstate, electron_in_bstate, photon_out_bstate, electron_out_bstate)
    matrix_elements = Vector{ComplexF64}()
    sizehint!(matrix_elements, length(base_states_comb))
    for (phin, ein, phout, eout) in base_states_comb
        push!(matrix_elements, _perturbative_compton_matrix(phin, ein, phout, eout))
    end

    matrix_elements_sq = abs2.(matrix_elements)

    # average over incoming polarizations/spins, but sum over outgoing pols/spins
    normalization = 1.0 / (length(photon_in_bstate) * length(electron_in_bstate))
    I = photon_in * electron_in

    return I * normalization * sum(matrix_elements_sq) * _phase_space_factor(photon_in, electron_in, photon_out, electron_out)
end

function _perturbative_compton_matrix(ph_in_bstate::SLorentzVector{ComplexF64}, el_in_bstate::BiSpinor, ph_out_bstate::SLorentzVector{ComplexF64}, el_out_bstate::AdjointBiSpinor)
    # TODO
    return zero(ComplexF64)
end

function _phase_space_factor(ph_in::NumericType, el_in::NumericType, ph_out::NumericType, el_out::NumericType) where {NumericType<:QEDbase.AbstractFourMomentum}
    # TODO
    return zero(ComplexF64)
end

function _post_process(
    process::Compton,
    model::PerturbativeQED,
    in_phase_space::AbstractArray{NumericType},
    out_phase_space::AbstractArray{NumericType},
    result::Real,
) where {NumericType<:QEDbase.AbstractFourMomentum}
    # generally nothing to be done here
    return result
end
