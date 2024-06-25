@inline function _pert_omega_prime(omega, cth; mass=1.0)
    return omega / (1 + omega / mass * (1 - cth))
end

function generate_momenta(
    proc::Compton,
    model::PerturbativeQED,
    in_ps_def::PhasespaceDefinition{SphericalCoordinateSystem,ElectronRestFrame},
    in_ps::NTuple{N,T},
    out_ps::NTuple{M,T},
) where {N,M,T<:Real}
    return QEDbase._generate_momenta(proc, model, in_ps_def, in_ps, out_ps)
end

function QEDbase._generate_incoming_momenta(
    proc::Compton,
    model::PerturbativeQED,
    in_ps_def::PhasespaceDefinition{SphericalCoordinateSystem,ElectronRestFrame},
    in_ps::NTuple{N,T},
) where {N,T<:Real}
    om = in_ps[1]

    P = SFourMomentum(one(om), zero(om), zero(om), zero(om))
    K = SFourMomentum(om, zero(om), zero(om), om)

    return P, K
end

function QEDbase._generate_momenta(
    proc::Compton,
    model::PerturbativeQED,
    in_ps_def::PhasespaceDefinition{SphericalCoordinateSystem,ElectronRestFrame},
    in_ps::NTuple{N,T},
    out_ps::NTuple{M,T},
) where {N,M,T<:Real}
    omega = in_ps[1]
    cth = out_ps[1]
    phi = out_ps[2]
    P, K, Pp, Kp = _generate_momenta_elab_sph(omega, cth, phi) # TODO: do this coord and frame dependent
    in_moms = (P, K)
    out_moms = (Pp, Kp)
    return in_moms, out_moms
end

function _generate_momenta_elab_sph(om, cth, phi, m=1.0)
    P = SFourMomentum(m, zero(m), zero(m), zero(m))
    K = SFourMomentum(om, zero(om), zero(om), om)
    omp = _pert_omega_prime(om, cth)
    sth = sqrt(1 - cth^2)
    sphi, cphi = sincos(phi)
    Kp = SFourMomentum(omp, omp * sth * cphi, omp * sth * sphi, omp * cth)
    Pp = P + K - Kp
    return P, K, Pp, Kp
end
