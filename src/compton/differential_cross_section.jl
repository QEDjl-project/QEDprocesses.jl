
function _assert_valid_input(
    ::Compton,
    ::PerturbativeQED,
    in_phase_space::AbstractArray{NumericType},
    out_phase_space::AbstractArray{NumericType},
) where {NumericType<:QEDbase.AbstractFourMomentum}
    @assert sum(in_phase_space) == sum(out_phase_space) "The sums of energies and momenta of in_phase_space and out_phase_space are not equal ($(sum(in_phase_space)) versus $(sum(out_phase_space)))"
    return nothing
end

function _differential_cross_section(
    ::Compton{NoPolarization},
    model::PerturbativeQED,
    in_phase_space::AbstractVector{NumericType},
    out_phase_space::AbstractVector{NumericType},
) where {NumericType<:QEDbase.AbstractFourMomentum}
    # TODO
    return 1.0
end

function _differential_cross_section(
    ::Compton{XPolarization},
    model::PerturbativeQED,
    in_phase_space::AbstractVector{NumericType},
    out_phase_space::AbstractVector{NumericType},
) where {NumericType<:QEDbase.AbstractFourMomentum}
    # TODO
    return 1.0
end

function _differential_cross_section(
    ::Compton{YPolarization},
    model::PerturbativeQED,
    in_phase_space::AbstractVector{NumericType},
    out_phase_space::AbstractVector{NumericType},
) where {NumericType<:QEDbase.AbstractFourMomentum}
    # TODO
    return 1.0
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
