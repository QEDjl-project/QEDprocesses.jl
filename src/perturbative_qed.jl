###############
# Perturbative QED
#
# In this file, we define the Perturbative QED type and
# implement the AbstractModelDefinition interface to provide its functionality.
###############

struct PerturbativeQED <: AbstractModelDefinition end

fundamental_interaction_type(::PerturbativeQED) = :electromagnetic
