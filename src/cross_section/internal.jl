
# convenience function for internal use only
# differential probability without energy momentum conservation check
# single in phase space points/ single out phase space point
# based on four-momenta
function _unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    psp = generate_phase_space(proc,model,phase_space_def,in_phase_space,out_phase_space)
    return unsafe_differential_probability(psp)
end

# convenience function for internal use only 
# differential probability without energy momentum conservation check
# based on coordinates
function _unsafe_differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta, out_momenta = _generate_momenta(
        proc, model, phase_space_def, in_phase_space, out_phase_space
    )
    return _unsafe_differential_probability(
        proc, model, phase_space_def, in_momenta, out_momenta
    )
end

# convenience function for internal use only 
# differential probability with energy momentum conservation check
# one in phase space point/ one out phase space point
# based on four-momenta
function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    psp = generate_phase_space(proc,model,phase_space_def,in_phase_space,out_phase_space)
    return differential_probability(psp)
end

# convenience function for internal use only 
# differential probability with energy momentum conservation check
# one in phase space point/ one out phase space point
# based on coordinates
function _differential_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta, out_momenta = _generate_momenta(
        proc, model, phase_space_def, in_phase_space, out_phase_space
    )
    return _differential_probability(proc, model, phase_space_def, in_momenta, out_momenta)
end

# convenience function for internal use only 
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
    psp = generate_phase_space(proc,model,phase_space_def,in_phase_space,out_phase_space)
    return unsafe_differential_cross_section(psp)
end

# convenience function for internal use only 
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
        proc, model, phase_space_def, in_phase_space, out_phase_space
    )
    return _unsafe_differential_cross_section(
        proc, model, phase_space_def, in_momenta, out_momenta
    )
end

# convenience function for internal use only 
# differential cross sections with energy momentum conservation check
# single in phase space point/ single out phase space point
function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    psp = generate_phase_space(proc,model,phase_space_def,in_phase_space,out_phase_space)
    return differential_cross_section(psp)
end

# convenience function for internal use only 
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
        proc, model, phase_space_def, in_phase_space, out_phase_space
    )
    return _differential_cross_section(
        proc, model, phase_space_def, in_momenta, out_momenta
    )
end
