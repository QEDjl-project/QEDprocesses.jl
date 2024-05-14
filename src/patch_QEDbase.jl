#############
# Patches for `QEDbase.jl`
# remove if this went into `QEDbase.jl`
#
# fix will be provided here: https://github.com/QEDjl-project/QEDbase.jl/pull/62
#############

Broadcast.broadcastable(dir::Incoming) = Ref(dir)
Broadcast.broadcastable(dir::Outgoing) = Ref(dir)
Broadcast.broadcastable(part::AbstractParticleType) = Ref(part)
Broadcast.broadcastable(spin_or_pol::AbstractSpinOrPolarization) = Ref(spin_or_pol)
