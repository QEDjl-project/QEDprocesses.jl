"""
    Compton(
        in_spin [= AllSpin()]
        in_pol [= AllPol()]
        out_spin [= AllSpin()]
        out_pol [= AllPol()]
    )
"""
struct Compton{InElectronSpin,InPhotonPol,OutElectronSpin,OutPhotonPol} <:
       AbstractProcessDefinition where {
    InElectronSpin<:AbstractSpin,
    InPhotonPol<:AbstractPolarization,
    OutElectronSpin<:AbstractSpin,
    OutPhotonPol<:AbstractPolarization,
}
    in_spin::InElectronSpin
    in_pol::InPhotonPol

    out_spin::OutElectronSpin
    out_pol::OutPhotonPol
end

function Compton()
    return Compton(AllSpin(), AllPol(), AllSpin(), AllPol())
end

function Compton(in_pol::AbstractPolarization)
    return Compton(AllSpin(), in_pol, AllSpin(), AllPol())
end
function Compton(in_pol::AbstractPolarization, out_pol::AbstractPolarization)
    return Compton(AllSpin(), in_pol, AllSpin(), out_pol)
end

_polarizations(proc::Compton) = (proc.in_pol, proc.out_pol)
_spins(proc::Compton) = (proc.in_spin, proc.out_spin)
_in_spin_and_pol(proc::Compton) = (proc.in_spin, proc.in_pol)
_out_spin_and_pol(proc::Compton) = (proc.out_spin, proc.out_pol)

function QEDbase.incoming_particles(::Compton)
    return (Electron(), Photon())
end

function QEDbase.outgoing_particles(::Compton)
    return (Electron(), Photon())
end

function QEDbase.incoming_spin_pols(proc::Compton)
    return (proc.in_spin, proc.in_pol)
end

function QEDbase.outgoing_spin_pols(proc::Compton)
    return (proc.out_spin, proc.out_pol)
end

function _spin_or_pol(process::Compton, ::Electron, ::Incoming)
    return process.in_spin
end

function _spin_or_pol(process::Compton, ::Electron, ::Outgoing)
    return process.out_spin
end

function _spin_or_pol(process::Compton, ::Photon, ::Incoming)
    return process.in_pol
end

function _spin_or_pol(process::Compton, ::Photon, ::Outgoing)
    return process.out_pol
end

function Base.show(io::IO, ::Compton)
    print(io, "one-photon Compton scattering")
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", proc::Compton)
    println(io, "one-photon Compton scattering")
    println(io, "    incoming: electron ($(proc.in_spin)), photon ($(proc.in_pol))")
    println(io, "    outgoing: electron ($(proc.out_spin)), photon ($(proc.out_pol))")
    return nothing
end
