
# function barrier for base_state
#
# WARN: This file can be removed if the implementations of `base_state` are removed from
# `QEDbase`.
#
# TL;DR: 
# The final state of the restructuing is, that QEDbase exports the symbol `base_state` and 
# QEDcore implements the concrete functions. However, currently, QEDbase provides both,
# the symbol and the implementation. Hiding this implementation in the QEDbase namespace 
# causes the usage of QEDbase.base_state implementation and not the QEDcore.base_state. 
# This is problematic, because the return type of `base_state` must be from `QEDcore`,
# which is not the case for `QEDbase.base_state`. Therefore, we need to recast the returns
# of `QEDbase.base_state` into types from `QEDcore`.
#
# This should not break if we update QEDbase by removing the concrete implementation of
# `base_state` in `QEDbase`, because `QEDcore` adds the implementations to the symbol 
# `QEDbase.base_state` anyway. 

base_state(p::QEDbase.Fermion,d::QEDbase.Incoming,mom,spin) = BiSpinor(QEDbase.base_state(p,d,mom,spin))
base_state(p::QEDbase.AntiFermion,d::QEDbase.Incoming,mom,spin) = AdjointBiSpinor(QEDbase.base_state(p,d,mom,spin))
base_state(p::QEDbase.Fermion,d::QEDbase.Outgoing,mom,spin) = AdjointBiSpinor(QEDbase.base_state(p,d,mom,spin))
base_state(p::QEDbase.AntiFermion,d::QEDbase.Outgoing,mom,spin) = BiSpinor(QEDbase.base_state(p,d,mom,spin))
base_state(p::QEDbase.BosonLike,d::QEDbase.ParticleDirection,mom,spin) = SLorentzVector(QEDbase.base_state(p,d,mom,spin))

base_state(p::QEDbase.Fermion,d::QEDbase.Incoming,mom,spin::QEDbase.AllSpin) = BiSpinor.(QEDbase.base_state(p,d,mom,spin))
base_state(p::QEDbase.AntiFermion,d::QEDbase.Incoming,mom,spin::QEDbase.AllSpin) = AdjointBiSpinor.(QEDbase.base_state(p,d,mom,spin))
base_state(p::QEDbase.Fermion,d::QEDbase.Outgoing,mom,spin::QEDbase.AllSpin) = AdjointBiSpinor.(QEDbase.base_state(p,d,mom,spin))
base_state(p::QEDbase.AntiFermion,d::QEDbase.Outgoing,mom,spin::QEDbase.AllSpin) = BiSpinor.(QEDbase.base_state(p,d,mom,spin))
base_state(p::QEDbase.BosonLike,d::QEDbase.ParticleDirection,mom,spin::QEDbase.AllPol) = SLorentzVector.(QEDbase.base_state(p,d,mom,spin))
