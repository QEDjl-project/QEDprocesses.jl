
function _assert_valid_dcs_input(
    process::Compton,
    ::PerturbativeQED,
    in_phase_space::AbstractArray{NumericType},
    out_phase_space::AbstractArray{NumericType},
) where {NumericType<:QEDbase.AbstractFourMomentum}
    is_onshell(out_phase_space[1], mass(outcoming_particles(process)[1])) || throw(
        InvalidInputError(
            "First out-phasespace particle must be an on-shell Photon\nValue: $(out_phase_space[1])",
        ),
    )
    is_onshell(out_phase_space[2], mass(outgoing_particles(process)[2])) || throw(
        InvalidInputError(
            "Second out-phasespace particle must be an on-shell Electron\nValue: $(out_phase_space[2])",
        ),
    )
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

    matrix_elements_sq = _matrix_el_sq(
        process, model, photon_in, electron_in, photon_out, electron_out
    )

    # average over incoming polarizations/spins, but sum over outgoing pols/spins
    normalization = 1.0 / (length(photon_in_bstate) * length(electron_in_bstate))
    I = photon_in * electron_in

    return I *
           normalization *
           sum(matrix_elements_sq) *
           _phase_space_factor(photon_in, electron_in, photon_out, electron_out)
end

function _perturbative_compton_matrix(
    ph_in::NumericType,
    el_in::NumericType,
    ph_out::NumericType,
    el_out::NumericType,
    ph_in_bstate::SLorentzVector{ComplexF64},
    el_in_bstate::BiSpinor,
    ph_out_bstate::SLorentzVector{ComplexF64},
    el_out_bstate::AdjointBiSpinor,
) where {NumericType<:QEDbase.AbstractFourMomentum}
    ph_in_slashed = slashed(ph_in_bstate)
    ph_out_slashed = slashed(ph_out_bstate)

    # TODO: fermion propagator is not yet in QEDbase
    diagram_1 =
        ph_out_slashed *
        _fermion_propagator(ph_in + el_in, mass(Electron())) *
        ph_in_slashed
    diagram_2 =
        ph_in_slashed *
        _fermion_propagator(el_in - ph_out, mass(Electron())) *
        ph_out_slashed

    result = diagram_1 + diagram_2
    result = result * el_in_bstate
    result = el_out_bstate * result

    # TODO: find (preferably unitful) global provider for physical constants
    # elementary charge
    alpha = 1 / 137.035999084
    e = sqrt(4 * pi * alpha)

    return -e * e * result
end

function _phase_space_factor(
    ph_in::NumericType, el_in::NumericType, ph_out::NumericType, el_out::NumericType
) where {NumericType<:QEDbase.AbstractFourMomentum}
    # TODO
    return zero(ComplexF64)
end

function _post_process_dcs(
    process::Compton,
    model::PerturbativeQED,
    in_phase_space::AbstractArray{NumericType},
    out_phase_space::AbstractArray{NumericType},
    result::Real,
) where {NumericType<:QEDbase.AbstractFourMomentum}
    # generally nothing to be done here
    return result
end

function _matrix_el(
    process::Compton{InPol,InSpin,OutPol,OutSpin},
    model::PerturbativeQED,
    photon_in::NumericType,
    electron_in::NumericType,
    photon_out::NumericType,
    electron_out::NumericType,
) where {
    InPol<:AbstractPolarization,
    InSpin<:AbstractSpin,
    OutPol<:AbstractPolarization,
    OutSpin<:AbstractSpin,
    NumericType<:QEDbase.AbstractFourMomentum,
}
    # get base states of the particles
    photon_in_bstate = Vector{SLorentzVector{ComplexF64}}(
        base_state(
            Photon(), Incoming(), photon_in, _spin_or_pol(process, Photon(), Incoming())
        ),
    )
    electron_in_bstate = Vector{BiSpinor}(
        base_state(
            Electron(),
            Incoming(),
            electron_in,
            _spin_or_pol(process, Electron(), Incoming()),
        ),
    )
    photon_out_bstate = Vector{SLorentzVector{ComplexF64}}(
        base_state(
            Photon(), Outgoing(), photon_out, _spin_or_pol(process, Photon(), Outgoing())
        ),
    )
    electron_out_bstate = Vector{AdjointBiSpinor}(
        base_state(
            Electron(),
            Outgoing(),
            electron_out,
            _spin_or_pol(process, Electron(), Outgoing()),
        ),
    )

    # if the particles had AllSpin or AllPol, the base states can be vectors and we need to consider every combination of the base states with each other
    base_states_comb = Iterators.product(
        photon_in_bstate, electron_in_bstate, photon_out_bstate, electron_out_bstate
    )
    matrix_elements = Vector{ComplexF64}()
    sizehint!(matrix_elements, length(base_states_comb))
    for (ph_in, el_in, ph_out, el_out) in base_states_comb
        push!(
            matrix_elements,
            _perturbative_compton_matrix(
                photon_in,
                electron_in,
                photon_out,
                electron_out,
                ph_in,
                el_in,
                ph_out,
                el_out,
            ),
        )
    end

    return matrix_elements
end

function _matrix_el_sq(
    process::Compton{InPol,InSpin,OutPol,OutSpin},
    model::PerturbativeQED,
    photon_in::NumericType,
    electron_in::NumericType,
    photon_out::NumericType,
    electron_out::NumericType,
) where {
    NumericType<:QEDbase.AbstractFourMomentum,
    InPol<:AbstractPolarization,
    InSpin<:AbstractSpin,
    OutPol<:AbstractPolarization,
    OutSpin<:AbstractSpin,
}
    return abs2.(
        _matrix_el(process, model, photon_in, electron_in, photon_out, electron_out)
    )
end
