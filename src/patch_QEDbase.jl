#############
# Patches for `QEDbase.jl`
# remove if this went into `QEDbase.jl`
#
#############

# fix: https://github.com/QEDjl-project/QEDbase.jl/pull/61
Base.show(io::IO, ::MIME"text/plain", ::Electron) = print(io, "electron")
Base.show(io::IO, ::MIME"text/plain", ::Positron) = print(io, "positron")
Base.show(io::IO, ::MIME"text/plain", ::Photon) = print(io, "photon")
Base.show(io::IO, ::MIME"text/plain", ::Incoming) = print(io, "incoming")
Base.show(io::IO, ::MIME"text/plain", ::Outgoing) = print(io, "outgoing")
Base.show(io::IO, ::MIME"text/plain", ::PolX) = print(io, "x-polarized")
Base.show(io::IO, ::MIME"text/plain", ::PolY) = print(io, "y-polarized")
Base.show(io::IO, ::MIME"text/plain", ::AllPol) = print(io, "all polarizations")
Base.show(io::IO, ::MIME"text/plain", ::SpinUp) = print(io, "spin up")
Base.show(io::IO, ::MIME"text/plain", ::SpinDown) = print(io, "spin down")
Base.show(io::IO, ::MIME"text/plain", ::AllSpin) = print(io, "all spins")

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
