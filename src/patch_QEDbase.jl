#############
# Patches for `QEDbase.jl`
# remove if this went into `QEDbase.jl`
#
#############

# fix: https://github.com/QEDjl-project/QEDbase.jl/pull/61
Base.show(io::IO, ::Electron) = print(io, "electron")
Base.show(io::IO, ::Positron) = print(io, "positron")
Base.show(io::IO, ::Photon) = print(io, "photon")
Base.show(io::IO, ::Incoming) = print(io, "incoming")
Base.show(io::IO, ::Outgoing) = print(io, "outgoing")
Base.show(io::IO, ::PolX) = print(io, "x-polarized")
Base.show(io::IO, ::PolY) = print(io, "y-polarized")
Base.show(io::IO, ::AllPol) = print(io, "all polarizations")
Base.show(io::IO, ::SpinUp) = print(io, "spin up")
Base.show(io::IO, ::SpinDown) = print(io, "spin down")
Base.show(io::IO, ::AllSpin) = print(io, "all spins")

# fix: https://github.com/QEDjl-project/QEDbase.jl/pull/62
Broadcast.broadcastable(dir::Incoming) = Ref(dir)
Broadcast.broadcastable(dir::Outgoing) = Ref(dir)
Broadcast.broadcastable(part::AbstractParticleType) = Ref(part)
Broadcast.broadcastable(spin_or_pol::AbstractSpinOrPolarization) = Ref(spin_or_pol)

# fix: https://github.com/QEDjl-project/QEDbase.jl/pull/63
number_of_spin_pol(::AbstractDefinitePolarization) = 1
number_of_spin_pol(::AbstractDefiniteSpin) = 1
number_of_spin_pol(::AbstractIndefinitePolarization) = 2
number_of_spin_pol(::AbstractIndefiniteSpin) = 2
