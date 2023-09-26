###############
# The process setup 
#
# In this file, we define the interface for process setups. A `setup` means
# here, a collection of setup-data needed to evaluate a dedicated quantity on a given
# running data. Despite that, the decomposition into setup and running data is
# arbitrary, this will be used for cases where a subset of the variables a
# quantity depends on is kept constant.  
# 
# This file is part of `QEDprocesses.jl` which is by itself part of the `QED.jl`
# ecosystem.
#
###############

# base type for any setup, will most properbly be defined somewhere else,
# e.g. QEDbase, therefore it is not exported
abstract type AbstractComputationSetup end

# convenience function to check if an object is a computation setup
_is_computation_setup(::AbstractComputationSetup) = true

"""

    _is_valid_input(stp::AbstractComputationSetup, input::Any)

Interface function, which returns true, if the constraints of the `input` associated with the quantity of `stp` are met.
This function will be called to validate the input of [`compute`](@ref) before calling [`_compute`](@ref).

!!! note "Default implementation"
    
    Since no input validation is equivalent to every input being valid, this functions returns `true` per default. 
    This behaviour needs to be overwritten, if this is not the case.

"""
@inline function _input_validation(stp::AbstractComputationSetup, input)
    return true
end

"""
    
    function _post_computation(stp::AbstractComputationSetup, input::Any, result::Any)
    
Interface functions which is called in [`compute`](@ref) after [`_compute`](@ref) was called. This function is dedicated to 
finalize the result of a computation. 

!!! note "default implementation"

    Since in the case of no post computation, the result of [`_compute`](@ref) is not changed, this function returns `result` per default.

"""
@inline function _post_computation(stp::AbstractComputationSetup, input, result)
    return result
end

"""
    
    _compute(stp::AbstractComputationSetup, input::Any)

Interface function which returns the value of the associated quantity evaluated on `input`, which can be anything the associated quantity is defined to be feasable for.

!!! note "unsafe implementation"

    This function should be implemented without any input validation or post computations (see [`_input_validation`](@ref) and [`_post_computation`](@ref)). The latter two are performed while calling 
    the safe version of this function [`compute`](@ref).

"""
function _compute end



function compute(stp::AbstractComputationSetup, input)
    _input_validation(stp,input) || error("InvalidInputError: there is something wrong with the input!\n setup:$stp \n input:$input")
    raw_result = _compute(stp,input)
    return _post_computation(stp, input,raw_result)
end

"""
Abstract base type for setups related to a combination scattering processes and compute models.  
Every subtype of `AbstractProcessSetup` is expected to implement at least the following 
interface functions

```Julia
scattering_process(::AbstractProcessSetup) 
compute_model(::AbstractProcessSetup) 
```
where `x` is the running input, i.e. a *single point* the quantity is evaluated on.

Derived from the these interface functions, the following delegations are implemented

```Julia
number_incoming_particles(::AbstractProcessSetup)
number_outgoing_particles(::AbstractProcessSetup)
```
"""
abstract type AbstractProcessSetup <: AbstractComputationSetup end

"""
    
    scattering_process(stp::AbstractProcessSetup)

Interface function, which returns the scattering process associated with `stp`,
i.e. an object which is a subtype of `AbstractProcessDefinition`. 
"""
function scattering_process end

"""

    compute_model(stp::AbstractProcessSetup)

Interface function, which returns the compute model associated with `stp`, i.e.
an object which is subtype of `AbstractModelDefinition`
"""
function compute_model end

@inline number_incoming_particles(stp::AbstractProcessSetup) = number_incoming_particles(scattering_process(stp))
@inline number_outgoing_particles(stp::AbstractProcessSetup) = number_outgoing_particles(scattering_process(stp))

