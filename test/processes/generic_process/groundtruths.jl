POLS = [PolX(), PolY(), AllPol()]
SPINS = [SpinUp(), SpinDown(), AllSpin()]

function _groundtruth_is_physical(proc::ScatteringProcess, ::PerturbativeQED)
    incoming_electrons = number_particles(proc, Incoming(), Electron())
    incoming_positrons = number_particles(proc, Incoming(), Positron())
    outgoing_electrons = number_particles(proc, Outgoing(), Electron())
    outgoing_positrons = number_particles(proc, Outgoing(), Positron())

    return incoming_electrons + outgoing_positrons ==
           outgoing_electrons + incoming_positrons
end

function _groundtruth_spin_pols(particles)
    return ntuple(
        x -> is_fermion(particles[x]) ? AllSpin() : AllPolarization(), length(particles)
    )
end

function _random_spin_pols(RNG, particles)
    return ntuple(
        x -> is_fermion(particles[x]) ? rand(RNG, SPINS) : rand(RNG, POLS),
        length(particles),
    )
end
