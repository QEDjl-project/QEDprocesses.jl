"""
    GenericQEDProcess <: AbstractProcessDefinition

An implementation of the `AbstractProcessDefinition` interface for a generic particle process in QED
with any number of incoming and outgoing particles, and any combination of spins or polarizations for the particles set.

The [`isphysical`](@ref) function can be used to check whether the process is possible in QED.

!!! note
    The computation of cross sections and probabilities is currently unimplemented.
"""
struct GenericQEDProcess{INT,OUTT,INSP,OUTSP} <:
       AbstractProcessDefinition where {INT<:Tuple,OUTT<:Tuple,INSP<:Tuple,OUTSP<:Tuple}
    incoming_particles::INT
    outgoing_particles::OUTT

    incoming_spin_pols::INSP
    outgoing_spin_pols::OUTSP

    """
        GenericQEDProcess(
            in_particles::Tuple{AbstractParticleType},
            out_particles::Tuple{AbstractParticleType},
            in_sp::Tuple{AbstractSpinOrPolarization},
            out_sp::Tuple{AbstractSpinOrPolarization}
        )

    Constructor for a GenericQEDProcess with the given incoming and outgoing particles and their respective spins and pols.
    The constructor asserts that the particles are compatible with their respective spins and polarizations. If the assertion fails, an 
    `InvalidInputError` is thrown.
    """
    function GenericQEDProcess(
        in_particles::INT, out_particles::OUTT, in_spin_pols::INSP, out_spin_pols::OUTSP
    ) where {INT<:Tuple,OUTT<:Tuple,INSP<:Tuple,OUTSP<:Tuple}
        _assert_particle_type_tuple(in_particles)
        _assert_particle_type_tuple(out_particles)

        _assert_spin_pol_particle_compatability(in_particles, in_spin_pols)
        _assert_spin_pol_particle_compatability(out_particles, out_spin_pols)

        return new{INT,OUTT,INSP,OUTSP}(
            in_particles, out_particles, in_spin_pols, out_spin_pols
        )
    end
end

"""
    GenericQEDProcess(in_particles::Tuple{AbstractParticleType}, out_particles::Tuple{AbstractParticleType})

Constructor for a GenericQEDProcess, setting `AllPol` and `AllSpin` for every boson and fermion, respectively.
"""
function GenericQEDProcess(
    in_particles::INT, out_particles::OUTT
) where {INT<:Tuple,OUTT<:Tuple}
    # this will be called again by the default constructor, but it produces a nicer warning here
    # than the following spin/pol generation failing because is_fermion or is_boson isn't defined on not allowed types
    _assert_particle_type_tuple(in_particles)
    _assert_particle_type_tuple(out_particles)

    in_spin_pols = ntuple(
        x -> is_fermion(in_particles[x]) ? AllSpin() : AllPolarization(),
        length(in_particles),
    )
    out_spin_pols = ntuple(
        x -> is_fermion(out_particles[x]) ? AllSpin() : AllPolarization(),
        length(out_particles),
    )
    return GenericQEDProcess(in_particles, out_particles, in_spin_pols, out_spin_pols)
end

function QEDbase.incoming_particles(proc::GenericQEDProcess)
    return proc.incoming_particles
end
function QEDbase.outgoing_particles(proc::GenericQEDProcess)
    return proc.outgoing_particles
end
function QEDbase.incoming_spin_pols(proc::GenericQEDProcess)
    return proc.incoming_spin_pols
end
function QEDbase.outgoing_spin_pols(proc::GenericQEDProcess)
    return proc.outgoing_spin_pols
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

function Base.show(io::IO, proc::GenericQEDProcess)
    print(io, "generic QED process \"")
    for p in incoming_particles(proc)
        print(io, _particle_to_letter(p))
    end
    print(io, " -> ")
    for p in outgoing_particles(proc)
        print(io, _particle_to_letter(p))
    end
    print(io, "\"")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", proc::GenericQEDProcess)
    println(io, "generic QED process")
    for dir in (Incoming(), Outgoing())
        first = true
        for (p, sp) in zip(particles(proc, dir), spin_pols(proc, dir))
            if !first
                print(io, ", ")
            else
                print(io, "    $(dir): ")
                first = false
            end
            print(io, "$(p) ($(sp))")
        end
        println(io)
    end
    return nothing
end
