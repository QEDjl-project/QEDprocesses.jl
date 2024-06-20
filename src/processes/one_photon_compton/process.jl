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
    InElectronSpin<:QEDbase.AbstractSpin,
    InPhotonPol<:QEDbase.AbstractPolarization,
    OutElectronSpin<:QEDbase.AbstractSpin,
    OutPhotonPol<:QEDbase.AbstractPolarization,
}
    in_spin::InElectronSpin
    in_pol::InPhotonPol

    out_spin::OutElectronSpin
    out_pol::OutPhotonPol
end

function Compton()
    return Compton(QEDbase.AllSpin(), QEDbase.AllPol(), QEDbase.AllSpin(), QEDbase.AllPol())
end

function Compton(in_pol::QEDbase.AbstractPolarization)
    return Compton(QEDbase.AllSpin(), in_pol, QEDbase.AllSpin(), QEDbase.AllPol())
end
function Compton(
    in_pol::QEDbase.AbstractPolarization, out_pol::QEDbase.AbstractPolarization
)
    return Compton(QEDbase.AllSpin(), in_pol, QEDbase.AllSpin(), out_pol)
end

_polarizations(proc::Compton) = (proc.in_pol, proc.out_pol)
_spins(proc::Compton) = (proc.in_spin, proc.out_spin)
_in_spin_and_pol(proc::Compton) = (proc.in_spin, proc.in_pol)
_out_spin_and_pol(proc::Compton) = (proc.out_spin, proc.out_pol)

function QEDprocesses.incoming_particles(::Compton)
    return (Electron(), Photon())
end

function QEDprocesses.outgoing_particles(::Compton)
    return (Electron(), Photon())
end

function _spin_or_pol(process::Compton, ::Electron, ::QEDbase.Incoming)
    return process.in_spin
end

function _spin_or_pol(process::Compton, ::Electron, ::QEDbase.Outgoing)
    return process.out_spin
end

function _spin_or_pol(process::Compton, ::Photon, ::QEDbase.Incoming)
    return process.in_pol
end

function _spin_or_pol(process::Compton, ::Photon, ::QEDbase.Outgoing)
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
