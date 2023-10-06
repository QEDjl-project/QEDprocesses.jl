
"""
    Compton{Pol} <: AbstractScatteringProcess where {Pol<:AbstractPolarization}

Type for a Compton scattering process. The Compton process is parametrized with the type of [`AbstractPolarization`](@ref) used.
This type implements the [`AbstractScatteringProcess`](@ref) interface. It can be used with [`DifferentialCrossSection`](@ref) to calculate differential cross sections of Compton scattering events.
"""
struct Compton{Pol} <: AbstractScatteringProcess where {Pol<:AbstractPolarization}
    polarizations::Tuple{Pol,Pol}
    spins::Tuple{AbstractSpin,AbstractSpin}
end

"""
    ComptonDCS

A convenience type specialization of [`DifferentialCrossSection`](@ref) with [`Compton`](@ref) set as the process type and [`PerturbativeQED`](@ref) as the model.
"""
const ComptonDCS = DifferentialCrossSection{
    Compton{Pol},
    PerturbativeQED,
    PhaseSpace,
} where {Pol<:AbstractPolarization,PhaseSpace<:AbstractArray{SFourMomentum}}

"""
    $(TYPEDSIGNATURES)

Default constructor of a [`Compton`](@ref) scattering process, using [`NoPolarization`](@ref).
"""
function Compton()
    return Compton(NoPolarization())
end

"""
    $(TYPEDSIGNATURES)

Return the incoming particles of a [`Compton`](@ref) scattering process which is a [`Photon`](@ref) and an [`Electron`](@ref).
"""
function incoming_particles(::Compton)
    return (Photon(), Electron())
end

"""
    $(TYPEDSIGNATURES)

Return the outgoing particles of a [`Compton`](@ref) scattering process which is a [`Photon`](@ref) and an [`Electron`](@ref).
"""
function outgoing_particles(::Compton)
    return (Photon(), Electron())
end

"""
    $(TYPEDSIGNATURES)

Return the implementation of [`AbstractPolarization`](@ref) set for this [`Compton`](@ref) scattering process.
"""
function polarization(process::Compton)
    return process.polarization
end
