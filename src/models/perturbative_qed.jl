
struct PerturbativeQED <: AbstractModelDefinition end

fundamental_interaction_type(::PerturbativeQED) = :electromagnetic

"""
    in_phase_space_dimension(proc::AbstractProcessDefinition, ::PerturbativeQED)


Return the number of degrees of freedom to determine the incoming phase space for processes in PerturbativeQED. 

!!! note "Convention"

    The current implementation only supports the case, where two of the incoming particles collide heads-on. 

"""
function in_phase_space_dimension(proc::AbstractProcessDefinition, ::PerturbativeQED)
    return 3 * number_incoming_particles(proc) - 4 - 1
end

function out_phase_space_dimension(proc::AbstractProcessDefinition, ::PerturbativeQED)
    return 3 * number_outgoing_particles(proc) - 4
end
