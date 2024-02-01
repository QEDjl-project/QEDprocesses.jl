############
# scattering probabilities
#
# This file contains implementations of the scattering probability based on the
# process interface with and without input validation and/or phase space
# constraint.
############

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

function _unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
) where {T<:QEDbase.AbstractFourMomentum}
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

function _unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
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

Return differential probability without checking if the given phase space(s) are physical.
"""
function unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    size(in_phase_space, 1) == number_incoming_particles(proc) || throw(
        DimensionMismatch(
            "The number of incoming particles <$(number_incoming_particles(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

    size(out_phase_space, 1) == number_outgoing_particles(proc) || throw(
        DimensionMismatch(
            "The number of outgoing particles <$(number_outgoing_particles(proc))> is inconsistent with input size <$(size(out_phase_space,1))>",
        ),
    )

    return _unsafe_differential_probability(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
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

function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
) where {T<:QEDbase.AbstractFourMomentum}
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

function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
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

Return the differential cross section if the given phase spaces are physical, and zero otherwise. 
"""
function differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    size(in_phase_space, 1) == number_incoming_particles(proc) || throw(
        DimensionMismatch(
            "The number of incoming particles <$(number_incoming_particles(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

    size(out_phase_space, 1) == number_outgoing_particles(proc) || throw(
        DimensionMismatch(
            "The number of outgoing particles <$(number_outgoing_particles(proc))> is inconsistent with input size <$(size(out_phase_space,1))>",
        ),
    )

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

function _total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
) where {T<:QEDbase.AbstractFourMomentum}
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

Return the total probability of a given model and process combination.
"""
function total_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    size(in_phase_space, 1) == number_incoming_particles(proc) || throw(
        DimensionMismatch(
            "The number of incoming particles <$(number_incoming_particles(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

    return _total_probability(proc, model, in_phase_space_def, in_phase_space)
end
