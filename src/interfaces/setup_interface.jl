###############
# The process setup 
#
# In this file, we define the interface for general computation and process setups.
# 
# This file is part of `QEDprocesses.jl` which is by itself part of the `QED.jl`
# ecosystem.
#
###############

"""
Abstract base type for computation setups.  A *setup* means
here, a collection of setup-data needed to evaluate a dedicated quantity on given
running data. Therefore, each setup is associated to a single quantity, which one may compute using the setup data and the funning data. 
Despite that, the decomposition into setup and running data is
arbitrary, and this can be used for cases where a subset of the variables a
quantity depends on is kept constant. 

!!! note "Computation setup interface"
    
    The computation performed using a computation setup is separated into three steps:
        
        1. input validation,
        2. actual computation,
        3. post computation,

    where every step has its own interface function (see [`compute`](@ref) for details). 
    
    ## input validation

    Every subtype of `AbstractComputationSetup` implements the interface function
    
    ```Julia
    _input_validation(stp::AbstractComputationSetup, input) # default: true
    ```
    
    which returns true, if the `input` is valid for the computation of the associated quantity (see [`_input_validation`](@ref) for more details). 
    The default implementation returns always `true`; provide a custom implementation if a different behavior is required.

    ## Actual computation
    
    Every subtype of `AbstractComputationSetup` should at least implement the following required interface function
    
    ```Julia
    _compute(stp::AbstractComputationSetup, input) 
    ```
    
    which computes the value of the associated quantity for a given `input` (see [`_compute`](@ref) for more details).


    ## Post computation

    Every subtype of `AbstractComputationSetup` implement the following interface function
    
    ```Julia
    _post_computation(stp::AbstractComputationSetup, input) 
    ```
    
    which performs *computations* after the actual computation, e.g. conversions or normalisations (see [`_post_computation`](@ref) for more details).


    
"""
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

"""

    compute(stp::AbstractComputationSetup, input::Any)

Return the value of the quantity associated with `stp` for a given `input`. 
In addition to the actual call of the associated unsafe version [`_compute`](@ref),
input validation (using [`_input_validation`](@ref)) and post computation
(using [`_post_computation`](@ref)) are wrapped around the calculation (see [`AbstractComputationSetup`](@ref) for details).
"""
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

