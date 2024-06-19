#############
# Patches for `QEDbase.jl`
# remove if this went into `QEDbase.jl`
#
#############

# fix: https://github.com/QEDjl-project/QEDbase.jl/pull/61
Base.show(io::IO, ::Electron) = print(io, "electron")
Base.show(io::IO, ::Positron) = print(io, "positron")
Base.show(io::IO, ::Photon) = print(io, "photon")
Base.show(io::IO, ::QEDbase.Incoming) = print(io, "incoming")
Base.show(io::IO, ::QEDbase.Outgoing) = print(io, "outgoing")
Base.show(io::IO, ::QEDbase.PolX) = print(io, "x-polarized")
Base.show(io::IO, ::QEDbase.PolY) = print(io, "y-polarized")
Base.show(io::IO, ::QEDbase.AllPol) = print(io, "all polarizations")
Base.show(io::IO, ::QEDbase.SpinUp) = print(io, "spin up")
Base.show(io::IO, ::QEDbase.SpinDown) = print(io, "spin down")
Base.show(io::IO, ::QEDbase.AllSpin) = print(io, "all spins")
 
#=
# fix: https://github.com/QEDjl-project/QEDbase.jl/pull/62
Broadcast.broadcastable(dir::QEDbase.Incoming) = Ref(dir)
Broadcast.broadcastable(dir::QEDbase.Outgoing) = Ref(dir)
Broadcast.broadcastable(part::QEDbase.AbstractParticleType) = Ref(part)
Broadcast.broadcastable(spin_or_pol::QEDbase.AbstractSpinOrPolarization) = Ref(spin_or_pol)

=#
# fix: https://github.com/QEDjl-project/QEDbase.jl/pull/63
number_of_spin_pol(::QEDbase.AbstractDefinitePolarization) = 1
number_of_spin_pol(::QEDbase.AbstractDefiniteSpin) = 1
number_of_spin_pol(::QEDbase.AbstractIndefinitePolarization) = 2
number_of_spin_pol(::QEDbase.AbstractIndefiniteSpin) = 2
