########################
# differential and total cross sections.
#
# This file contains default implementations for differential and total cross
# sections based on the scattering process interface
########################

############
# differential cross sections
############

# differential cross sections without energy momentum conservation check
# single in phase space point/ single out phase space point
# based on four-momenta
function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    I = 1 / (4 * _incident_flux(proc, model, in_phase_space))

    return I * _unsafe_differential_probability(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
end

# differential cross sections without energy momentum conservation check
# single in phase space point/ single out phase space point
# based on coordinates
function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta, out_momenta = _generate_momenta(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
    return _unsafe_differential_cross_section(
        proc, model, phase_space_def, in_momenta, out_momenta
    )
end

# differential cross sections without energy momentum conservation check
# single in phase space point/ several out phase space points
function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
    res = Vector{eltype(T)}(undef, size(out_phase_space, 2))
    for i in 1:size(out_phase_space, 2)
        res[i] = _unsafe_differential_cross_section(
            proc,
            model,
            phase_space_def,
            in_phase_space,
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
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space::AbstractVecOrMat{T},
) where {T<:AbstractPhasespaceElement}
    res = Matrix{eltype(T)}(undef, size(in_phase_space, 2), size(out_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i, :] .= _unsafe_differential_cross_section(
            proc,
            model,
            phase_space_def,
            view(in_phase_space, :, i),
            out_phase_space,
        )
    end
    return res
end

"""

    function unsafe_differential_cross_section(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
        out_phase_space::AbstractVecOrMat{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

Return the differential cross section evaluated at the four-momenta without checking if the given phase space(s) are physical.
"""
function unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    _check_in_phase_space_dimension(proc, model, in_phase_space)
    _check_out_phase_space_dimension(proc, model, out_phase_space)

    return _unsafe_differential_cross_section(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
end

"""
    unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}

Return the differential cross section evaluated at the coordinates without checking if the given phase space(s) are physical.
"""
function unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    _check_in_phase_space_dimension(proc, model, in_phase_space)
    _check_out_phase_space_dimension(proc, model, out_phase_space)

    return _unsafe_differential_cross_section(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
end

function unsafe_differential_cross_section(phase_space_point::PhaseSpacePoint)
  return _unsafe_differential_cross_section(
    phase_space_point.proc,
    phase_space_point.model,
    phase_space_point.ps_def,
    momentum.(phase_space_point.in_particles),
    momentum.(phase_space_point.out_particles)
  )
end

# differential cross sections with energy momentum conservation check
# single in phase space point/ single out phase space point
function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    if !_is_in_phasespace(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
        return zero(eltype(T))
    end

    return _unsafe_differential_cross_section(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
end

# differential cross sections with energy momentum conservation check
# single in phase space point/ several out phase space points
function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:Real}
    in_momenta, out_momenta = _generate_momenta(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
    return _differential_cross_section(
        proc, model, phase_space_def, in_momenta, out_momenta
    )
end

function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
    res = Vector{eltype(T)}(undef, size(out_phase_space, 2))
    for i in 1:size(out_phase_space, 2)
        res[i] = _differential_cross_section(
            proc,
            model,
            phase_space_def,
            in_phase_space,
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
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space::AbstractVecOrMat{T},
) where {T<:AbstractPhasespaceElement}
    res = Matrix{eltype(T)}(undef, size(in_phase_space, 2), size(out_phase_space, 2))
    for i in 1:size(in_phase_space, 2)
        res[i, :] .= _differential_cross_section(
            proc,
            model,
            phase_space_def,
            view(in_phase_space, :, i),
            out_phase_space,
        )
    end
    return res
end

"""
    differential_cross_section(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
        out_phase_space::AbstractVecOrMat{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

If the given phase spaces are physical, return differential cross section evaluated at the four-momenta. Zero otherwise.

"""
function differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    _check_in_phase_space_dimension(proc, model, in_phase_space)
    _check_out_phase_space_dimension(proc, model, out_phase_space)

    return _differential_cross_section(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
end

"""
    differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}

If the given phase spaces are physical, return differential cross section evaluated at the coordinates. Zero otherwise.
"""
function differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    _check_in_phase_space_dimension(proc, model, in_phase_space)
    _check_out_phase_space_dimension(proc, model, out_phase_space)

    return _differential_cross_section(
        proc,
        model,
        phase_space_def,
        in_phase_space,
        out_phase_space,
    )
end

function differential_cross_section(phase_space_point::PhaseSpacePoint)
  return _differential_cross_section(
    phase_space_point.proc,
    phase_space_point.model,
    phase_space_point.ps_def,
    momentum.(phase_space_point.in_particles),
    momentum.(phase_space_point.out_particles)
  )
end

############
# Total cross sections
############

# total cross section on single phase space point
# based on four-momenta
function _total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    I = 1 / (4 * _incident_flux(proc, model, in_phase_space))

    return I * _total_probability(proc, model, phase_space_def, in_phase_space)
end

# total cross section on single phase space point
# based on coordinates
function _total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta = _generate_incoming_momenta(proc, model, phase_space_def, in_phase_space)
    return _total_cross_section(proc, model, phase_space_def, in_momenta)
end

# total cross section on several phase space points
function _total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
) where {T<:AbstractPhasespaceElement}
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

Return the total cross section for a given combination of scattering process and compute model, evaluated at the particle momenta.
"""
function total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    _check_in_phase_space_dimension(proc, model, in_phase_space)

    return _total_cross_section(proc, model, in_phase_space_def, in_phase_space)
end

"""
    total_cross_section(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition,
        in_phase_space_def::AbstractPhasespaceDefinition,
        in_phase_space::AbstractVecOrMat{T},
    ) where {T<:Real}

Return the total cross section for a given combination of scattering process and compute model, evaluated at the coordinates.
"""
function total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    _check_in_phase_space_dimension(proc, model, in_phase_space)

    return _total_cross_section(proc, model, in_phase_space_def, in_phase_space)
end
