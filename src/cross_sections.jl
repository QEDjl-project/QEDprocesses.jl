########################
# differential cross sections and probabilities.
#
# This file contains default implementations for differential and total cross
# sections based on the scattering process interface
########################

############
#
# differential cross sections
#
############

# differential cross sections without energy momentum conservation check
# single in phase space point/ single out phase space point
function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    I = 1 / (4 * _incident_flux(proc, model, in_phase_space))

    return I * _unsafe_differential_probability(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

# differential cross sections without energy momentum conservation check
# single in phase space point/ several out phase space points
function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
    res = Vector{eltype(T)}(undef, size(out_phase_space, 2))
    for i in 1:size(out_phase_space, 2)
        res[i] = _unsafe_differential_cross_section(
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

# differential cross sections without energy momentum conservation check
# several in phase space points/ one or several out phase space points
function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:AbstractPhasespaceElement}
    res = Matrix{eltype(T)}(undef, size(in_phase_space, 2), size(out_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i, :] .= _unsafe_differential_cross_section(
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

    function unsafe_differential_cross_section(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
        out_phase_space_def::AbstractPhasespaceDefinition,
        out_phase_space::AbstractVecOrMat{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

Return the differential cross section without checking if the given phase space(s) are physical.
"""
function unsafe_differential_cross_section(
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

    return _unsafe_differential_cross_section(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

function unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    size(in_phase_space, 1) == in_phase_space_dimension(proc, model) || throw(
        DimensionMismatch(
            "The dimension of the in-phase-space <$(in_phase_space_dimension(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

    size(out_phase_space, 1) == out_phase_space_dimension(proc) || throw(
        DimensionMismatch(
            "The dimension of the out-phase-space <$(out_phase_space_dimension(proc))> is inconsistent with input size <$(size(out_phase_space,1))>",
        ),
    )

    return _unsafe_differential_cross_section(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

# differential cross sections with energy momentum conservation check
# single in phase space point/ single out phase space point
function _differential_cross_section(
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

    return _unsafe_differential_cross_section(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

# differential cross sections with energy momentum conservation check
# single in phase space point/ several out phase space points
function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:Real}
    in_momenta, out_momenta = _generate_momenta(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
    return _differential_cross_section(
        proc, model, in_phase_space_def, in_momenta, out_phase_space_def, out_momenta
    )
end

function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
    res = Vector{eltype(T)}(undef, size(out_phase_space, 2))
    for i in 1:size(out_phase_space, 2)
        res[i] = _differential_cross_section(
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

# differential cross sections with energy momentum conservation check
# several in phase space points/ one or several out phase space points
function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:AbstractPhasespaceElement}
    res = Matrix{eltype(T)}(undef, size(in_phase_space, 2), size(out_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i, :] .= _differential_cross_section(
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
    differential_cross_section(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
        out_phase_space_def::AbstractPhasespaceDefinition,
        out_phase_space::AbstractVecOrMat{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

If the given phase spaces are physical, return differential cross section. Zero otherwise

"""
function differential_cross_section(
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

    return _differential_cross_section(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

function differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    size(in_phase_space, 1) == in_phase_space_dimension(proc, model) || throw(
        DimensionMismatch(
            "The dimension of the in-phase-space <$(in_phase_space_dimension(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

    size(out_phase_space, 1) == out_phase_space_dimension(proc) || throw(
        DimensionMismatch(
            "The dimension of the out-phase-space <$(out_phase_space_dimension(proc))> is inconsistent with input size <$(size(out_phase_space,1))>",
        ),
    )

    return _differential_cross_section(
        proc,
        model,
        in_phase_space_def,
        in_phase_space,
        out_phase_space_def,
        out_phase_space,
    )
end

############
# Total cross sections
############

# total cross section on single phase space point
function _total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    I = 1 / (4 * _incident_flux(proc, model, in_phase_space))

    return I * _total_probability(proc, model, in_phase_space_def, in_phase_space)
end

# total cross section on several phase space points
function _total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
) where {T<:QEDbase.AbstractFourMomentum}
    res = Vector{eltype(T)}(undef, size(in_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i] = _total_cross_section(
            proc, model, in_phase_space_def, view(in_phase_space, :, i)
        )
    end
    return res
end

"""
    total_cross_section(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

Return the total cross section for a given combination of scattering process and compute model.
"""
function total_cross_section(
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

    return _total_cross_section(proc, model, in_phase_space_def, in_phase_space)
end

