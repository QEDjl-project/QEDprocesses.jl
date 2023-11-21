###########
# utility functions
#
# This file contains small helper and utility functions used throughout the package.
###############

"""

    _base_component_type(array_of_lv::AbstractArray{LV}) where {LV<:QEDbase.AbstractLorentzVector}

Return the type of the components of given Lorentz vectors, which are by themself elements of an 
`AbstractArray`.

# Examples
```julia
julia> using QEDbase
julia> using QEDprocesses
julia> v = Vector{SFourMomentum}(undef,10)
julia> QEDprocesses._base_component_type(v)
Float64
```
"""
function _base_component_type(::AbstractArray{LV}) where {LV<:QEDbase.AbstractLorentzVector}
    return eltype(LV)
end


