mutable struct GenericQEDProcess{INT,OUTT,INSP,OUTSP} <: AbstractProcessDefinition where {
    INT<:Tuple,OUTT<:Tuple,INSP<:Tuple,OUTSP<:Tuple
}
    incoming_particles::INT
    outgoing_particles::OUTT

    incoming_spins_pols::INSP
    outgoing_spins_pols::OUTSP

    matrix_element_squared::Function

    function GenericQEDProcess(
        in_particles::INT, out_particles::OUTT, in_sp::INSP, out_sp::OUTSP
    ) where {INT<:Tuple,OUTT<:Tuple,INSP<:Tuple,OUTSP<:Tuple}
        _assert_particle_type_tuple(in_particles)
        _assert_particle_type_tuple(out_particles)

        return new{INT,OUTT,INSP,OUTSP}(
            in_particles, out_particles, in_sp, out_sp, _ -> error("unimplemented")
        )
    end
end

function QEDbase.incoming_particles(proc::GenericQEDProcess)
    return proc.incoming_particles
end
function QEDbase.outgoing_particles(proc::GenericQEDProcess)
    return proc.outgoing_particles
end
function QEDbase.incoming_spin_pols(proc::GenericQEDProcess)
    return proc.incoming_spins_pols
end
function QEDbase.outgoing_spin_pols(proc::GenericQEDProcess)
    return proc.outgoing_spins_pols
end

"""
    isphysical(proc::GenericQEDProcess)

A utility function that returns whether a given GenericQEDProcess conserves the number and charge of fermions and has at least 2 participating particles.
"""
function isphysical(proc::GenericQEDProcess)
    return (
        number_particles(proc, Incoming(), Electron()) +
        number_particles(proc, Outgoing(), Positron()) ==
        number_particles(proc, Incoming(), Positron()) +
        number_particles(proc, Outgoing(), Electron())
    ) && number_particles(proc, Incoming()) + number_particles(proc, Outgoing()) >= 2
end
