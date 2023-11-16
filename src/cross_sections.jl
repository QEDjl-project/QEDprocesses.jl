########################
# differential cross sections and probabilities.
#
# This file contains default implementations for differential and total cross
# sections based on the scattering process interface
########################


############
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
# differential cross sections
#
############

# differential cross sections without energy momentum conservation check
function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}

    I = 1/(4*incident_flux(proc,model,in_phase_space))

    return I * _unsafe_probability(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,out_phase_space)
end

function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    res = Vector{eltype(T)}(undef,size(out_phase_space,2))
    for i in 1:size(out_phase_space,2)
        res[i] = _unsafe_differential_cross_section(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,view(out_phase_space,:,i))
    end
    return res
end

function _unsafe_differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    res = Vector{eltype(T)}(undef,size(in_phase_space,2))
    for i in 1:size(in_phase_space,2)
        res[i] = _unsafe_differential_cross_section(proc,model,in_phase_space_def,view(in_phase_space,:,i),out_phase_space_def,out_phase_space)
    end
    return res
end

function _unsafe_differential_cross_section(
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
            res[i,j] = _unsafe_differential_cross_section(proc,model,in_phase_space_def,view(in_phase_space,:,i),out_phase_space_def,view(out_phase_space,:,j))
        end
    end
    return res
end

function unsafe_differential_cross_section(
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
        
    return _unsafe_differential_cross_section(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,out_phase_space)
end

# differential cross sections with energy momentum conservation check
function _differential_cross_section(
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

    return _unsafe_differential_cross_section(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,out_phase_space)
end

function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractMatrix{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    res = Vector{eltype(T)}(undef,size(out_phase_space,2))
    for i in 1:size(out_phase_space,2)
        res[i] = _differential_cross_section(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,view(out_phase_space,:,i))
    end
    return res
end

function _differential_cross_section(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_phase_space_def::AbstactPhasespaceDefinition,
    in_phase_space::AbstractMatrix{T},
    out_phase_space_def::AbstactPhasespaceDefinition,
    out_phase_space::AbstractVector{T},
)::Float64 where {T<:QEDbase.AbstractFourMomentum}
    
    res = Vector{eltype(T)}(undef,size(in_phase_space,2))
    for i in 1:size(in_phase_space,2)
        res[i] = _differential_cross_section(proc,model,in_phase_space_def,view(in_phase_space,:,i),out_phase_space_def,out_phase_space)
    end
    return res
end

function _differential_cross_section(
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
            res[i,j] = _differential_cross_section(proc,model,in_phase_space_def,view(in_phase_space,:,i),out_phase_space_def,view(out_phase_space,:,j))
        end
    end
    return res
end


function differential_cross_section(
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
        
    return _differential_cross_section(proc,model,in_phase_space_def,in_phase_space,out_phase_space_def,out_phase_space)
end

