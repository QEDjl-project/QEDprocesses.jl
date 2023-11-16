########################
# differential and total cross sections.
#
# This file contains default implementations for differential and total cross
# sections based on the scattering process interface
########################

############
<<<<<<< HEAD
=======
#
# scattering probability
#
############

function _unsafe_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    

    matrix_elements_sq =
        _matrix_element_square(proc, model, in_phase_space, out_phase_space)

    normalization = averaging_norm(proc)

    ps_fac = _phase_space_factor(proc,model,in_phasespace_def,in_phase_space,out_phasespace_def,out_phase_space)

    return normalization *
           sum(matrix_elements_sq) *
            ps_fac
end

function _unsafe_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    res = Vector{eltype(T)}(undef,size(out_phase_space,2))
    for i in 1:size(out_phase_space,2)
        res[i] = _unsafe_probability(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,view(out_phase_space,:,i))
    end
    return res
end

function _unsafe_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    res = Vector{eltype(T)}(undef,size(in_phase_space,2))
    for i in 1:size(in_phase_space,2)
        res[i] = _unsafe_probability(proc,model,in_phase_space_def,view(in_phase_space,:,i),out_phase_space_def,out_phase_space)
    end
    return res
end

function _unsafe_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    res = Matrix{eltype(T)}(undef,size(in_phase_space,2), size(out_phase_space,2))
    for i in 1:size(in_phase_space,2)
        for j in 1:size(out_phase_space,2)
            res[i,j] = _unsafe_probability(proc,model,in_phase_space_def,view(in_phase_space,:,i),out_phase_space_def,view(out_phase_space,:,j))
        end
    end
    return res
end

function unsafe_probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    size(in_ps,1)==number_incoming_paricles(proc) || throw(
        InvalidInputError("The number of incoming particles <$(number_incoming_paricles(proc))> is inconsistent with input size <$(size(in_ps,1))>"),
    )
        
    size(out_ps,1)==number_outgoing_paricles(proc) || throw(
        InvalidInputError("The number of outgoing particles <$(number_outgoing_paricles(proc))> is inconsistent with input size <$(size(out_ps,1))>"),
    )
        
    return _unsafe_probability(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,out_phase_space)
end
    
function _probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}

    # consider wrapping the unchecked diffCS in a function
    if (!isapprox(sum(in_phase_space), sum(out_phase_space); rtol = sqrt(eps())))
        return zero(eltype(T))
    end

    return _unsafe_probability(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,out_phase_space)
end

function _probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    res = Vector{eltype(T)}(undef,size(out_phase_space,2))
    for i in 1:size(out_phase_space,2)
        res[i] = _probability(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,view(out_phase_space,:,i))
    end
    return res
end

function _probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    res = Vector{eltype(T)}(undef,size(in_phase_space,2))
    for i in 1:size(in_phase_space,2)
        res[i] = _probability(proc,model,in_phase_space_def,view(in_phase_space,:,i),out_phase_space_def,out_phase_space)
    end
    return res
end

function _probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    res = Matrix{eltype(T)}(undef,size(in_phase_space,2), size(out_phase_space,2))
    for i in 1:size(in_phase_space,2)
        for j in 1:size(out_phase_space,2)
            res[i,j] = _probability(proc,model,in_phase_space_def,view(in_phase_space,:,i),out_phase_space_def,view(out_phase_space,:,j))
        end
    end
    return res
end


function probability(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    size(in_ps,1)==number_incoming_paricles(proc) || throw(
        InvalidInputError("The number of incoming particles <$(number_incoming_paricles(proc))> is inconsistent with input size <$(size(in_ps,1))>"),
    )
        
    size(out_ps,1)==number_outgoing_paricles(proc) || throw(
        InvalidInputError("The number of outgoing particles <$(number_outgoing_paricles(proc))> is inconsistent with input size <$(size(out_ps,1))>"),
    )
        
    return _probability(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,out_phase_space)
end

############
#
>>>>>>> 148f33c (fixed typos in function signatures and docstrings)
# differential cross sections
############

# differential cross sections without energy momentum conservation check
# single in phase space point/ single out phase space point
# based on four-momenta
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
# single in phase space point/ single out phase space point
# based on coordinates
function _unsafe_differential_cross_section(
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
    return _unsafe_differential_cross_section(
        proc, model, in_phase_space_def, in_momenta, out_phase_space_def, out_momenta
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

Return the differential cross section evaluated at the four-momenta without checking if the given phase space(s) are physical.
"""
function unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
<<<<<<< HEAD
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:QEDbase.AbstractFourMomentum}
    size(in_phase_space, 1) == number_incoming_particles(proc) || throw(
        DimensionMismatch(
            "The number of incoming particles <$(number_incoming_particles(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
=======
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    size(in_ps,1)==number_incoming_paricles(proc) || throw(
        InvalidInputError("The number of incoming particles <$(number_incoming_paricles(proc))> is inconsistent with input size <$(size(in_ps,1))>"),
>>>>>>> 148f33c (fixed typos in function signatures and docstrings)
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

"""
    unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}

Return the differential cross section evaluated at the coordinates without checking if the given phase space(s) are physical.
"""
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
            "The dimension of the in-phase-space <$(in_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

    size(out_phase_space, 1) == out_phase_space_dimension(proc, model) || throw(
        DimensionMismatch(
            "The dimension of the out-phase-space <$(out_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(out_phase_space,1))>",
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

If the given phase spaces are physical, return differential cross section evaluated at the four-momenta. Zero otherwise.

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

"""
    differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
    out_phase_space_def::AbstractPhasespaceDefinition,
    out_phase_space::AbstractVecOrMat{T},
) where {T<:Real}

If the given phase spaces are physical, return differential cross section evaluated at the coordinates. Zero otherwise.
"""
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
            "The dimension of the in-phase-space <$(in_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

    size(out_phase_space, 1) == out_phase_space_dimension(proc, model) || throw(
        DimensionMismatch(
            "The dimension of the out-phase-space <$(out_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(out_phase_space,1))>",
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
# based on four-momenta
function _total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    I = 1 / (4 * _incident_flux(proc, model, in_phase_space))

    return I * _total_probability(proc, model, in_phase_space_def, in_phase_space)
end

# total cross section on single phase space point
# based on coordinates
function _total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:Real}
    in_momenta = _generate_incoming_momenta(proc, model, in_phase_space_def, in_phase_space)
    return _total_cross_section(proc, model, in_phase_space_def, in_momenta)
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
    size(in_phase_space, 1) == number_incoming_particles(proc) || throw(
        DimensionMismatch(
            "The number of incoming particles <$(number_incoming_particles(proc))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

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
function QEDprocesses.total_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVecOrMat{T},
) where {T<:Real}
    size(in_phase_space, 1) == in_phase_space_dimension(proc, model) || throw(
        DimensionMismatch(
            "The dimension of the in-phase-space <$(in_phase_space_dimension(proc,model))> is inconsistent with input size <$(size(in_phase_space,1))>",
        ),
    )

    return _total_cross_section(proc, model, in_phase_space_def, in_phase_space)
end
