"""
    ScatteringProcess <: AbstractProcessDefinition

Generic implementation for scattering processes of arbitrary particles. Currently, only calculations in combination with `PerturbativeQED` are supported. 
However, this is supposed to describe scattering processes with any number of incoming and outgoing particles, and any combination of spins or polarizations for the particles.

The [`isphysical`](@ref) function can be used to check whether the process is possible in perturbative QED.

!!! warning  
    The computation of cross sections and probabilities is currently unimplemented.

## Constructors

    ScatteringProcess(
        in_particles::Tuple{AbstractParticleType},
        out_particles::Tuple{AbstractParticleType},
        [in_sp::Tuple{AbstractSpinOrPolarization},
        out_sp::Tuple{AbstractSpinOrPolarization}]
    )

Constructor for a ScatteringProcess with the given incoming and outgoing particles and their respective spins and pols.
The constructor asserts that the particles are compatible with their respective spins and polarizations. If the assertion fails, an 
`InvalidInputError` is thrown.

The `in_sp` and `out_sp` parameters can be omitted in which case all spins and polarizations will be set to `AllSpin` and `AllPol` for every fermion and boson, respectively.
"""
struct ScatteringProcess{INT,OUTT,INSP,OUTSP} <:
       AbstractProcessDefinition where {INT<:Tuple,OUTT<:Tuple,INSP<:Tuple,OUTSP<:Tuple}
    incoming_particles::INT
    outgoing_particles::OUTT

    incoming_spin_pols::INSP
    outgoing_spin_pols::OUTSP

    function ScatteringProcess(
        in_particles::NTuple{I,AbstractParticleType},
        out_particles::NTuple{O,AbstractParticleType},
        in_spin_pols::NTuple{I,AbstractSpinOrPolarization},
        out_spin_pols::NTuple{O,AbstractSpinOrPolarization},
    ) where {I,O}
        _assert_spin_pol_particle_compatability(in_particles, in_spin_pols)
        _assert_spin_pol_particle_compatability(out_particles, out_spin_pols)

        return new{
            typeof(in_particles),
            typeof(out_particles),
            typeof(in_spin_pols),
            typeof(out_spin_pols),
        }(
            in_particles, out_particles, in_spin_pols, out_spin_pols
        )
    end
end

function ScatteringProcess(
    in_particles::NTuple{I,AbstractParticleType},
    out_particles::NTuple{O,AbstractParticleType},
) where {I,O}
    in_spin_pols = ntuple(
        x -> is_fermion(in_particles[x]) ? AllSpin() : AllPolarization(),
        length(in_particles),
    )
    out_spin_pols = ntuple(
        x -> is_fermion(out_particles[x]) ? AllSpin() : AllPolarization(),
        length(out_particles),
    )
    return ScatteringProcess(in_particles, out_particles, in_spin_pols, out_spin_pols)
end

function QEDbase.incoming_particles(proc::ScatteringProcess)
    return proc.incoming_particles
end
function QEDbase.outgoing_particles(proc::ScatteringProcess)
    return proc.outgoing_particles
end
function QEDbase.incoming_spin_pols(proc::ScatteringProcess)
    return proc.incoming_spin_pols
end
function QEDbase.outgoing_spin_pols(proc::ScatteringProcess)
    return proc.outgoing_spin_pols
end

function Base.show(io::IO, proc::ScatteringProcess)
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

function Base.show(io::IO, ::MIME"text/plain", proc::ScatteringProcess)
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
