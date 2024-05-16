#############
# Patches for `QEDbase.jl`
# remove if this went into `QEDbase.jl`
#
#############

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
