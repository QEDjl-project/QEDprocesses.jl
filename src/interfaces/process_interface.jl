###############
# The process interface
#
# In this file, we define the interface for working with scattering processes in
# general.
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

Furthermore, to calculate scatteing probabilities and differential cross sections, the following 
interface functions need to be implemented

```Julia
    _incident_flux(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition, 
        in_ps::AbstractVector{T}
        ) where {T<:QEDbase.AbstractFourMomentum}

    _matrix_element(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition, 
        in_ps::AbstractVector{T}
        out_ps::AbstractVector{T}
        ) where {T<:QEDbase.AbstractFourMomentum}

    _averaging_norm(
        proc::AbstractProcessDefinition
        )

    _is_in_phasespace(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition, 
        in_ps_def::InPhasespaceDefinition,
        in_ps::AbstractVector{T}
        out_ps_def::OutPhasespaceDefinition,
        out_ps::AbstractVector{T}
    )

    _phase_space_factor(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition, 
        in_ps_def::InPhasespaceDefinition,
        in_ps::AbstractVector{T}
        out_ps_def::OutPhasespaceDefinition,
        out_ps::AbstractVector{T}
        ) where {T<:QEDbase.AbstractFourMomentum}
```

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
    _incident_flux(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition, 
        in_ps::AbstractVector{T}
        ) where {T<:QEDbase.AbstractFourMomentum}

Interface function, which returns the incident flux of the given scatteing process for a given incoming phase-space.

"""
function _incident_flux end

"""
    _matrix_element(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition, 
        in_ps::AbstractVector{T}
        out_ps::AbstractVector{T}
        ) where {T<:QEDbase.AbstractFourMomentum}

Interface function, which returns a tuple of scattering matrix elements for each spin and polarisation combination of `proc`. 
"""
function _matrix_element end

function _matrix_element_square(
    proc::AbstractProcessDefinition,
    model::AbstractModelDefinition,
    in_ps::AbstractVector{T},
    out_ps::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    mat_el = _matrix_element(proc, model, in_ps, out_ps)
    return abs2.(mat_el)
end

"""
    _averaging_norm(
        proc::AbstractProcessDefinition
        )

Interface function, which returns a normalization for the avering squared matrix elements over spins and polarizations. 
"""
function _averaging_norm end

"""

    _is_in_phasespace(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition, 
        in_ps_def::InPhasespaceDefinition,
        in_ps::AbstractVector{T}
        out_ps_def::OutPhasespaceDefinition,
        out_ps::AbstractVector{T}
        ) where {T<:QEDbase.AbstractFourMomentum}

Interface function, which returns `true`, if the combination of the given incoming and outgoing phase space
is physical, i.e. all momenta are on-shell and some sort of energy-momentum conservation holds.
"""
function _is_in_phasespace end

"""
    _phase_space_factor(
        proc::AbstractProcessDefinition,
        model::AbstractModelDefinition, 
        in_ps_def::InPhasespaceDefinition,
        in_ps::AbstractVector{T}
        out_ps_def::OutPhasespaceDefinition,
        out_ps::AbstractVector{T}
        ) where {T<:QEDbase.AbstractFourMomentum}

Interface function, which returns the pre-differential factor of the invariant phase space intergral measure. 

!!! note "Convention"

    It is assumed, that this function returns the value of 

    ```math
    \\mathrm{d}\\Pi_n:= \\prod_{i=1}^N \\frac{\\mathrm{p}^3p_i}{(2\\pi)^3 2 p_i^0} H(P_t, p_1, \\dots, p_N),
    ```
where ``H(\\dots)`` is the constraining of the phase space, e.g. ``\\delta^{(4)}(P_t - \\sum_i p_i)``.  
"""
function _phase_space_factor end

#######################
#
# utility functions
#
#######################

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

#
#     _total_cross_section(
#         proc_def::AbstractProcessDefinition,
#         model_def::AbstractModelDefinition,
#         in_phasespace::AbstractVector{T},
#     ) where {T<:QEDbase.AbstractFourMomentum} end
#
# Interface function for the combination of scattering processes and physical models. Return the total cross section of a 
# given process and model for a passed initial phase space. The elements of the `AbstractVector` representing the initial phase space
# are the momenta of the respective particles. The implementation of this function for a concrete process and model must not 
# check if the length of the passed initial phase spaces match number of incoming particles. 
#
# !!! note "cross section interface"
#
#     Given an implementation of this method, the following *unsafe* generic implementation is provided:
#
#     ```julia
#
#     _total_cross_section(proc_def,model_def,in_phasespace::AbstractMatrix{T})
#
#     ```
#
#     where `T<:QEDbase.AbstractFourMomentum`. Although, `_total_cross_section` is also implemented for a vector of initial phase space points.
#     Furthermore, a safe version of `_total_cross_section` is also implemented: [`total_cross_section`](@ref).
#
#
# !!! note 
#     
#     Each instance of this function does not check the validity of the input. 
#     This function is not exported and should be used with caution. To add a method in order to implement the cross section interface, 
#     it is recommented to directly use `QEDprocesses._total_cross_section` instead of globally `using QEDprocesses: _total_cross_section`.
#
# """
# function _total_cross_section end
#
# function _total_cross_section(
#     proc_def::AbstractProcessDefinition,
#     model_def::AbstractModelDefinition,
#     in_phasespace::AbstractMatrix{T},
# ) where {T<:QEDbase.AbstractFourMomentum}
#     res = Vector{_base_component_type(in_phasespace)}(undef, size(in_phasespace, 2))
#     for i in 1:size(in_phasespace, 2)
#         res[i] = _total_cross_section(proc_def, model_def, view(in_phasespace, :, i))
#     end
#     return res
# end
#
# """
#
#     total_cross_section(
#         proc_def::AbstractProcessDefinition,
#         model_def::AbstractModelDefinition,
#         in_phasespace::Union{AbstractVector{T},AbstractMatrix{T}},
#     ) where {T<:QEDbase.AbstractFourMomentum}
#
# Return the total cross section for a combination of a scattering process and a physical model evaluated on a given initial phase space. 
#
# This function will eventually call the respective interface function [`_total_cross_section`](@ref).
#
# """
# function total_cross_section(
#     proc_def::AbstractProcessDefinition,
#     model_def::AbstractModelDefinition,
#     in_phasespace::Union{AbstractVector{T},AbstractMatrix{T}},
# ) where {T<:QEDbase.AbstractFourMomentum}
#     size(in_phasespace, 1) == number_incoming_particles(proc_def) || throw(
#         DimensionMismatch(
#             "The number of momenta in the initial phasespace <{length(in_phasespace)}> does not match the number of incoming particles of the process <{number_incoming_particles(proc_def)}>.",
#         ),
#     )
#     return _total_cross_section(proc_def, model_def, in_phasespace)
# end
