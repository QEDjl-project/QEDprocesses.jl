
struct PerturbativeQED <: AbstractPerturbativeModel end

QEDbase.fundamental_interaction_type(::PerturbativeQED) = :electromagnetic

"""
    in_phase_space_dimension(proc::AbstractProcessDefinition, ::PerturbativeQED)

Return the number of degrees of freedom to determine the incoming phase space for processes in PerturbativeQED. 

!!! note "Convention"

    The current implementation only supports the case where two of the incoming particles collide head-on. 
"""
function QEDbase.in_phase_space_dimension(
    proc::AbstractProcessDefinition, ::PerturbativeQED
)
    return 3 * number_incoming_particles(proc) - 4 - 1
end

function QEDbase.out_phase_space_dimension(
    proc::AbstractProcessDefinition, ::PerturbativeQED
)
    return 3 * number_outgoing_particles(proc) - 4
end

function Base.show(io::IO, ::PerturbativeQED)
    print(io, "perturbative QED")
    return nothing
end

"""
    isphysical(proc::AbstractProcessDefinition, model::PerturbativeQED)

A utility function that returns whether a given `AbstractProcessDefinition` conserves the number and charge of fermions and has at least 2 participating particles.
"""
function isphysical(proc::AbstractProcessDefinition, ::PerturbativeQED)
    return (
        number_particles(proc, Incoming(), Electron()) +
        number_particles(proc, Outgoing(), Positron()) ==
        number_particles(proc, Incoming(), Positron()) +
        number_particles(proc, Outgoing(), Electron())
    ) && number_particles(proc, Incoming()) + number_particles(proc, Outgoing()) >= 2
end
