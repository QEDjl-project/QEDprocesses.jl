
"""
    Compton{InPhotonPol,InElectronSpin,OutPhotonPol,OutElectronSpin} <: AbstractScatteringProcess where {
        InPhotonPol<:AbstractPolarization,
        InElectronSpin<:AbstractSpin,
        OutPhotonPol<:AbstractPolarization,
        OutElectronSpin<:AbstractSpin,
    }

Type for a Compton scattering process. The Compton process is parametrized with the types of [`AbstractPolarization`](@ref) and [`AbstractSpin`](@ref) used for the incoming and outgoing particles.
This type implements the [`AbstractScatteringProcess`](@ref) interface. It can be used with [`DifferentialCrossSection`](@ref) to calculate differential cross sections of Compton scattering events.
"""
struct Compton{InPhotonPol,InElectronSpin,OutPhotonPol,OutElectronSpin} <: AbstractScatteringProcess where {
    InPhotonPol<:AbstractPolarization,
    InElectronSpin<:AbstractSpin,
    OutPhotonPol<:AbstractPolarization,
    OutElectronSpin<:AbstractSpin,
}
    polarizations::Tuple{InPhotonPol,OutPhotonPol}
    spins::Tuple{InElectronSpin,OutElectronSpin}
end

"""
    ComptonDCS

A convenience type specialization of [`DifferentialCrossSection`](@ref) with [`Compton`](@ref) set as the process type and [`PerturbativeQED`](@ref) as the model.
"""
const ComptonDCS = DifferentialCrossSection{
    Compton{P1,S1,P2,S2},
    PerturbativeQED,
    PhaseSpace,
} where {P1<:AbstractPolarization,S1<:AbstractSpin,P2<:AbstractPolarization,S2<:AbstractSpin,PhaseSpace<:AbstractArray{SFourMomentum}}

"""
    $(TYPEDSIGNATURES)

Default constructor of a [`Compton`](@ref) scattering process, using [`AllPolarization`](@ref) and [`AllSpin`](@ref) for all particles.
"""
function Compton()
    return Compton((AllPol(), AllPol()), (AllSpin(), AllSpin()))
end

"""
    $(TYPEDSIGNATURES)

Return the incoming particles of a [`Compton`](@ref) scattering process which is a [`Photon`](@ref) and an [`Electron`](@ref).

!!! note "Input order"
    The order of the types in the tuple this function returns is also the order that the [`SFourMomentum`](@ref) of the particles must be given in when calling [`compute`](@ref).
"""
function incoming_particles(::Compton)
    return (Photon(), Electron())
end

"""
    $(TYPEDSIGNATURES)

Return the outgoing particles of a [`Compton`](@ref) scattering process which is a [`Photon`](@ref) and an [`Electron`](@ref).

!!! note "Input order"
    The order of the types in the tuple this function returns is also the order that the [`SFourMomentum`](@ref) of the particles must be given in when calling [`compute`](@ref).
"""
function outgoing_particles(::Compton)
    return (Photon(), Electron())
end

"""
    _spin_or_pol(process::Compton, ::Particle, ::ParticleDirection)

Return the implementation of [`AbstractPolarization`](@ref) or [`AbstractSpin`] set for this [`Compton`](@ref) scattering process, depending on the particle type ([`Electron`](@ref) or [`Photon`](@ref) and its [`ParticleDirection`](@ref)).
"""
function _spin_or_pol end

function _spin_or_pol(process::Compton, ::Electron, ::Incoming)
    return process.spin[1]
end

function _spin_or_pol(process::Compton, ::Electron, ::Outgoing)
    return process.spin[2]
end

function _spin_or_pol(process::Compton, ::Photon, ::Incoming)
    return process.polarization[1]
end

function _spin_or_pol(process::Compton, ::Photon, ::Outgoing)
    return process.polarization[2]
end
