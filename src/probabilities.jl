############
# scattering probabilities
#
# This file contains implementations of the scattering probability based on the
# process interface with and without input validation and/or phase space
# constraint.
############

# differential probability without energy momentum conservation check
# single in phase space points/ single out phase space point
# based on four-momenta
function _unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    matrix_elements_sq = _matrix_element_square(
        proc, model, in_phase_space, out_phase_space
    )

    normalization = _averaging_norm(proc)

    ps_fac = _phase_space_factor(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )

    return normalization * sum(matrix_elements_sq) * ps_fac
end

# differential probability without energy momentum conservation check
# based on coordinates
function _unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta, out_momenta = _generate_momenta(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
    return _unsafe_differential_probability(
        proc, model, in_phase_space_def, in_momenta, out_phase_space_def, out_momenta
    )
end

# differential probability without energy momentum conservation check
# single in phase space points/ several out phase space point
function _unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
    res = Vector{eltype(T)}(undef, size(out_phase_space, 2))
    for i in 1:size(out_phase_space, 2)
        res[i] = _unsafe_differential_probability(
            proc,
            model,
            in_phase_space_def,
            in_phase_space,
            out_phase_space_def,
            view(out_phase_space, :, i),
        )
    end
    return res
end

# differential probability without energy momentum conservation check
# several in phase space points/ one or several out phase space point
function _unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:AbstractPhasespaceElement}
    res = Matrix{eltype(T)}(undef, size(in_phase_space, 2), size(out_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i, :] .= _unsafe_differential_probability(
            proc,
            model,
            in_phase_space_def,
            view(in_phase_space, :, i),
            out_phase_space_def,
            out_phase_space,
        )
    end
    return res
end

"""
    unsafe_differential_probability(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
        out_phase_space_def::AbstractPhasespaceDefinition,
        out_phase_space::AbstractVecOrMat{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

Return differential probability evaluated at the four-momenta without checking if the given phase space(s) are physical.
"""
function unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    _check_in_phase_space_dimension(proc, model, in_phase_space)
    _check_out_phase_space_dimension(proc, model, out_phase_space)

    return _unsafe_differential_probability(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

"""
    unsafe_differential_probability(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
        out_phase_space_def::AbstractPhasespaceDefinition,
        out_phase_space::AbstractVecOrMat{T},
    ) where {T<:Real}

Return differential probability evaluated at the coordinates without checking if the given phase space(s) are physical.
"""
function unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    _check_in_phase_space_dimension(proc, model, in_phase_space)
    _check_out_phase_space_dimension(proc, model, out_phase_space)

    return _unsafe_differential_probability(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

# differential probability with energy momentum conservation check
# one in phase space point/ one out phase space point
# based on four-momenta
function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:AbstractPhasespaceElement}
    if !_is_in_phasespace(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
        return zero(eltype(T))
    end

    return _unsafe_differential_probability(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

# differential probability with energy momentum conservation check
# one in phase space point/ one out phase space point
# based on coordinates
function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta, out_momenta = _generate_momenta(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
    return _differential_probability(
        proc, model, in_phase_space_def, in_momenta, out_phase_space_def, out_momenta
    )
end

# differential probability with energy momentum conservation check
# one in phase space points/ several out phase space point
function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta, out_momenta = _generate_momenta(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
    return _differential_probability(
        proc, model, in_phase_space_def, in_momenta, out_phase_space_def, out_momenta
    )
end

function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
    res = Vector{eltype(T)}(undef, size(out_phase_space, 2))
    for i in 1:size(out_phase_space, 2)
        res[i] = _differential_probability(
            proc,
            model,
            in_phase_space_def,
            in_phase_space,
            out_phase_space_def,
            view(out_phase_space, :, i),
        )
    end
    return res
end

# differential probability with energy momentum conservation check
# several in phase space points/ one or several out phase space point
function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:AbstractPhasespaceElement}
    res = Matrix{eltype(T)}(undef, size(in_phase_space, 2), size(out_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i, :] .= _differential_probability(
            proc,
            model,
            in_phase_space_def,
            view(in_phase_space, :, i),
            out_phase_space_def,
            out_phase_space,
        )
    end
    return res
end

"""
    differential_probability(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
        out_phase_space_def::AbstractPhasespaceDefinition,
        out_phase_space::AbstractVecOrMat{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

If the given phase spaces are physical, return differential probability evaluated at the four-momenta. Zero otherwise.
"""
function differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    _check_in_phase_space_dimension(proc, model, in_phase_space)
    _check_out_phase_space_dimension(proc, model, out_phase_space)

    return _differential_probability(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

"""
    differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}

If the given phase spaces are physical, return differential probability evaluated at the coordinates. Zero otherwise.
"""
function differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    _check_in_phase_space_dimension(proc, model, in_phase_space)
    _check_out_phase_space_dimension(proc, model, out_phase_space)

    return _differential_probability(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

###########
# Total probability
###########

# total probability on a phase space point
# based on coordinates
function _total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta = _generate_incoming_momenta(proc, model, in_phase_space_def, in_phase_space)
    return _total_probability(proc, model, in_phase_space_def, in_momenta)
end

# total probability on several phase space points
function _total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta = _generate_incoming_momenta(proc, model, in_phase_space_def, in_phase_space)
    return _total_probability(proc, model, in_phase_space_def, in_momenta)
end

function _total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
    res = Vector{eltype(T)}(undef, size(in_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i] = _total_probability(
            proc, model, in_phase_space_def, view(in_phase_space, :, i)
        )
    end
    return res
end

"""
    total_probability(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractMatrix{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

Return the total probability of a given model and process combination, evaluated at the particle momenta.
"""
function total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    _check_in_phase_space_dimension(proc, model, in_phase_space)

    return _total_probability(proc, model, in_phase_space_def, in_phase_space)
end

"""
    total_probability(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractMatrix{T},
    ) where {T<:Real}

Return the total probability of a given model and process combination, evaluated at the coordinates.
"""
function total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    _check_in_phase_space_dimension(proc, model, in_phase_space)

    return _total_probability(proc, model, in_phase_space_def, in_phase_space)
end
