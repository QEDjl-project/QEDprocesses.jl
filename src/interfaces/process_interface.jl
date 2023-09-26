###############
# The process interface
#
# In this file, we define the interface for working with scattering processes in
# general.
# 
# This file is part of `QEDprocesses.jl` which is by itself part of the `QED.jl`
# ecosystem.
#
###############

"""
Abstract base type for definitions of scattering processes. It is the root type for the 
process interface, which assumes that every subtype of `AbstractProcessDefinition`
implements at least 

```Julia
incoming_particles(proc_def::AbstractProcessDefinition)
outgoing_particles(proc_def::AbstractProcessDefinition)
```

which return a tuple of the incoming and outgoing particles, respectively.
"""
abstract type AbstractProcessDefinition end

"""

    incoming_particles(proc_def::AbstractProcessDefinition)

Interface function for scattering processes. Return a tuple of the incoming particles for the given process definition.
This function needs to be given to implement the scattering process interface.
"""
function incoming_particles end

"""

    outgoing_particles(proc_def::AbstractProcessDefinition)

Interface function for scattering processes. Return the tuple of outgoing particles for the given process definition.
This function needs to be given to implement the scattering process interface.
"""
function outgoing_particles end

"""

    $(TYPEDSIGNATURES)

Return the number of incoming particles of a given process. 
"""
@inline function number_incoming_particles(proc_def::AbstractProcessDefinition)
    return length(incoming_particles(proc_def))
end

"""

    $(TYPEDSIGNATURES)

Return the number of outgoing particles of a given process. 
"""
@inline function number_outgoing_particles(proc_def::AbstractProcessDefinition)
    return length(outgoing_particles(proc_def))
end

"""

    _differential_cross_section(
        proc_def::AbstractProcessDefinition,
        model_def::AbstractModelDefinition,
        init_phasespace::AbstractVector{T},
        final_phasespace::AbstractVector{T},
    ) where {T<:QEDbase.AbstractFourMomentum}

Interface function for the combination of scattering processes and models. Return the differential cross section of a 
given process and model for a passed initial and final phase space. The elements of the `AbstractVector` representing the phase spaces 
are the momenta of the respective particles. The implementation of this function for a concrete process and model must not 
check if the length of the passed phase spaces match the respective number of particles. 

!!! note "differential cross section interface"

    Given an implementation of this method, the following *unsafe* generic implementations are provided:

    ```julia

    _differential_cross_section(proc_def,model_def,init_phasespace::AbstractVector{T},finial_phasespace::AbstractMatrix{T})
    _differential_cross_section(proc_def,model_def,init_phasespace::AbstractMatrix{T},finial_phasespace::AbstractVector{T})
    _differential_cross_section(proc_def,model_def,init_phasespace::AbstractMatrix{T},finial_phasespace::AbstractMatrix{T})

    ```

    where `T<:QEDbase.AbstractFourMomentum`. Although, any combinations of initial and final phase space types given by *single set of points* (AbstractVector{T}) and *mutiple set of points* (AbstractMatrix{T}) 
    is implemented. Furthermore, a safe version of `_differential_cross_section` is also implemented: [`differential_cross_section`](@ref).

!!! note "unsafe implementation"
    
    Each instance of this function does not check the validity of the input. 
    Therefore, these functions are not exported and should be used with caution. To add a method in order to implement the cross section interface, 
    it is recommented to directly use `QEDprocesses._differential_cross_section` instead of globally `using QEDprocesses: _differential_cross_section`.


"""
function _differential_cross_section end

"""
    
    differential_cross_section(
        proc_def::AbstractProcessDefinition,
        model_def::AbstractModelDefinition,
        init_phasespace::Union{AbstractVector{T},AbstractMatrix{T}},
        final_phasespace::Union{AbstractVector{T},AbstractMatrix{T}},
    ) where {T<:QEDbase.AbstractFourMomentum}

Return the differential cross section for a given combination of a scattering process 
and model definition evaluated on the passed inital and final phase space points. 

This function will eventually call the respective interface function [`_differential_cross_section`](@ref).
"""
function differential_cross_section(
    proc_def::AbstractProcessDefinition,
    model_def::AbstractModelDefinition,
    init_phasespace::Union{AbstractVector{T},AbstractMatrix{T}},
    final_phasespace::Union{AbstractVector{T},AbstractMatrix{T}},
    ) where {T<:QEDbase.AbstractFourMomentum}
    size(init_phasespace, 1) == number_incoming_particles(proc_def) || throw(
        DimensionMismatch("The number of momenta in the initial phasespace <{length(init_phasespace)}> does not match the number of incoming particles of the process <{number_incoming_pariticles(proc_def)}>."),
    )
    size(final_phasespace, 1) == number_outgoing_particles(proc_def) || throw(
        DimensionMismatch("The number of momenta in the final phasespace <{length(final_phasespace)}> does not match the number of outgoing particles of the process <{number_outgoing_pariticles(proc_def)}>."),
    )
    return _differential_cross_section(proc_def, model_def, init_phasespace, final_phasespace)
end

# returns diffCS for single `initPS` and several `finalPS` points without input-check
function _differential_cross_section(
    proc_def::AbstractProcessDefinition,
    model_def::AbstractModelDefinition,
    init_phasespace::AbstractVector{T},
    final_phasespace::AbstractMatrix{T},
) where {T<:QEDbase.AbstractFourMomentum}
    res = Vector{_base_component_type(init_phasespace)}(undef, size(final_phasespace, 2))
    for i = 1:size(final_phasespace, 2)
        res[i] = _differential_cross_section(
            proc_def,
            model_def,
            init_phasespace,
            view(final_phasespace, :, i),
        )
    end
    return res
end

function _differential_cross_section(
    proc_def::AbstractProcessDefinition,
    model_def::AbstractModelDefinition,
    init_phasespace::AbstractMatrix{T},
    final_phasespace::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    res = Vector{_base_component_type(init_phasespace)}(undef, size(init_phasespace, 2))
    for i = 1:size(init_phasespace, 2)
        res[i] = _differential_cross_section(
            proc_def,
            model_def,
            view(init_phasespace, :, i),
            final_phasespace,
        )
    end
    return res
end

function _differential_cross_section(
    proc_def::AbstractProcessDefinition,
    model_def::AbstractModelDefinition,
    init_phasespace::AbstractMatrix{T},
    final_phasespace::AbstractMatrix{T},
) where {T<:QEDbase.AbstractFourMomentum}
    res = Matrix{_base_component_type(init_phasespace)}(
        undef,
        size(init_phasespace, 2),
        size(final_phasespace, 2),
    )
    for init_idx = 1:size(init_phasespace, 2)
        for final_idx = 1:size(final_phasespace, 2)
            res[init_idx, final_idx] = _differential_cross_section(
                proc_def,
                model_def,
                view(init_phasespace, :, init_idx),
                view(final_phasespace, :, final_idx),
            )
        end
    end
    return res
end

"""

    _total_cross_section(
        proc_def::AbstractProcessDefinition,
        model_def::AbstractModelDefinition,
        init_phasespace::AbstractVector{T},
    ) where {T<:QEDbase.AbstractFourMomentum} end

Interface function for the combination of scattering processes and models. Return the total cross section of a 
given process and model for a passed initial phase space. The elements of the `AbstractVector` representing the initial phase space
are the momenta of the respective particles. The implementation of this function for a concrete process and model must not 
check if the length of the passed initial phase spaces match number of incoming particles. 

!!! note "cross section interface"

    Given an implementation of this method, the following *unsafe* generic implementation is provided:

    ```julia

    _total_cross_section(proc_def,model_def,init_phasespace::AbstractMatrix{T})

    ```

    where `T<:QEDbase.AbstractFourMomentum`. Although, `_total_cross_section` is also implemented for a vector of initial phase space points.
    Furthermore, a safe version of `_total_cross_section` is also implemented: [`total_cross_section`](@ref).


!!! note 
    
    Each instance of this function does not check the validity of the input. 
    This function is not exported and should be used with caution. To add a method in order to implement the cross section interface, 
    it is recommented to directly use `QEDprocesses._total_cross_section` instead of globally `using QEDprocesses: _total_cross_section`.

"""
function _total_cross_section end

function _total_cross_section(
    proc_def::AbstractProcessDefinition,
    model_def::AbstractModelDefinition,
    init_phasespace::AbstractMatrix{T},
) where {T<:QEDbase.AbstractFourMomentum}
    res = Vector{_base_component_type(init_phasespace)}(undef, size(init_phasespace, 2))
    for i = 1:size(init_phasespace, 2)
        res[i] = _total_cross_section(proc_def, model_def, view(init_phasespace, :, i))
    end
    return res
end

"""

    total_cross_section(
        proc_def::AbstractProcessDefinition,
        model_def::AbstractModelDefinition,
        init_phasespace::Union{AbstractVector{T},AbstractMatrix{T}},
    ) where {T<:QEDbase.AbstractFourMomentum}

Return the total cross section for a combination of a scattering process and a compute model evaluated on a given initial phase space. 

This function will eventually call the respective interface function [`_total_cross_section`](@ref).

"""
function total_cross_section(
    proc_def::AbstractProcessDefinition,
    model_def::AbstractModelDefinition,
    init_phasespace::Union{AbstractVector{T},AbstractMatrix{T}},
) where {T<:QEDbase.AbstractFourMomentum}
    size(init_phasespace, 1) == number_incoming_particles(proc_def) || throw(
        DimensionMismatch("The number of momenta in the initial phasespace <{length(init_phasespace)}> does not match the number of incoming particles of the process <{number_incoming_pariticles(proc_def)}>."),
    )
    return _total_cross_section(proc_def, model_def, init_phasespace)
end
