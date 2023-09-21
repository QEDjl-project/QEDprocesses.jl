###############
# The model interface
#
# In this file, we define the interface of working with compute models in
# general.
# 
# This file is part of `QEDprocesses.jl` which is by itself part of the `QED.jl`
# ecosystem.
#
###############
# root type for models
"""
Abstract base type for all compute model definitions in the context of scattering processes. Every subtype of `AbstractModelDefinition` is associated with a fundamental interaction. 
Therefore, one needs to implement the following soft interface function

```Julia
fundamental_interaction_type(::AbstractModelDefinition)
```
"""
abstract type AbstractModelDefinition end

"""

    fundamental_interaction_type(models_def::AbstractModelDefinition)

Return the fundamental interaction associated with the passed model definition.
"""
function fundamental_interaction_type end
