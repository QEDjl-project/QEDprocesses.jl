
"""
    _generate_momenta(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVector{T},
        out_phase_space_def::AbstractPhasespaceDefinition,
        out_phase_space::AbstractVector{T},
    ) where {T<:Real}
"""
function _generate_momenta(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta = _generate_incoming_momenta(proc, model, in_phase_space_def, in_phase_space)
    out_momenta = _generate_outgoing_momenta(
        proc, model, out_phase_space_def, out_phase_space
    )

    return in_momenta, out_momenta
end
