
###########
# Total probability
###########

# total probability on a phase space point
# based on coordinates
function _total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta = _generate_incoming_momenta(proc, model, phase_space_def, in_phase_space)
    return _total_probability(proc, model, phase_space_def, in_momenta)
end

# total probability on several phase space points
function _total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
    res = Vector{eltype(T)}(undef, size(in_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i] = _total_probability(
            proc, model, phase_space_def, view(in_phase_space, :, i)
        )
    end
    return res
end

"""
    total_probability(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractMatrix{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

Return the total probability of a given model and process combination, evaluated at the particle momenta.
"""
function total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    _check_in_phase_space_dimension(proc, model, in_phase_space)

    return _total_probability(proc, model, phase_space_def, in_phase_space)
end

"""
    total_probability(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractMatrix{T},
    ) where {T<:Real}

Return the total probability of a given model and process combination, evaluated at the coordinates.
"""
function total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    _check_in_phase_space_dimension(proc, model, in_phase_space)

    return _total_probability(proc, model, phase_space_def, in_phase_space)
end
