
struct Compton{Pol} <: AbstractScatteringProcess where {Pol<:AbstractPolarization}
    polarization::Pol
end

const ComptonDCS = DifferentialCrossSection{
    Compton{Pol},
    PerturbativeQED,
    PhaseSpace,
} where {
    Pol<:AbstractPolarization,
    PhaseSpace<:AbstractArray{SFourMomentum},
}

function Compton()
    return Compton(NoPolarization())
end

function incoming_particles(::AbstractScatteringProcess)
    return (Photon(), Electron())
end

function outgoing_particles(::AbstractScatteringProcess)
    return (Photon(), Electron())
end

function polarization(process::Compton)
    return process.polarization
end
