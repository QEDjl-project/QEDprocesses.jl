
function _pert_omega_prime_elab_sph(om,cth)
    om/(1+om*(1-cth))
end

"""
Klein-Nishina: differential cross section

source: Peskin, Schroeder. "Quantum field theory." (1995). eq: 5.91
note: we compute d sigma/ d Omega, and *not* d sigma/ d cos theta (the difference is a factor 2 pi)

"""
function _groundtruth_pert_compton_diffCS_spinsum_polsum_elab_sph(om, cth, mass)
    prefac = ALPHA_SQUARE / 2
    omp = _pert_omega_prime_elab_sph(om, cth)
    sth_sq = 1 - cth^2
    return prefac * (omp / om)^2 * (omp / om + om / omp - sth_sq)
end

function _groundtruth_pert_compton_diffCS_spinsum_xpol_elab_sph(omega, ctheta, phi, mass)
    om_prime = _pert_omega_prime_elab_sph(omega, ctheta)
    om_prime_over_om = om_prime / omega
    return 0.5 * ALPHA_SQUARE / mass^2 *
           om_prime_over_om^2 *
           (om_prime_over_om + 1.0 / om_prime_over_om - 2 * (1 - ctheta^2) * cos(phi)^2)
end

function _groundtruth_pert_compton_diffCS_spinsum_ypol_elab_sph(omega, ctheta, phi, mass)
    return _groundtruth_pert_compton_diffCS_spinsum_xpol_elab_sph(omega, ctheta, phi + pi / 2, mass)
end
