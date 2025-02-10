
# incoming phase space layout
# alias for electron rest system (based on QEDcore.TwoBodyRestSystem)

const ComptonRestSystem{COORD} = TwoBodyRestSystem{1,COORD} where {COORD}
ComptonRestSystem(coord) = ComptonRestSystem(Val(particle_index(coord)), coord)
function ComptonRestSystem(
    run_idx::Val{1}, coord::COORD
) where {COORD<:QEDcore.AbstractSingleParticleCoordinate}
    throw(
        ArgumentError(
            "the first incoming particle of Compton is an electron, which has no degrees of freedom in the Compton rest frame.",
        ),
    )
end
function ComptonRestSystem(
    run_idx::Val{2}, coord::COORD
) where {COORD<:QEDcore.AbstractSingleParticleCoordinate{2}}
    return TwoBodyRestSystem{1}(coord)
end
ComptonRestSystem(coord::CMSEnergy) = TwoBodyRestSystem{1}(coord)
function ComptonRestSystem(::Rapidity)
    throw(
        ArgumentError(
            "there is no finite rapidity for a photon in the electron rest system"
        ),
    )
end
ComptonRestSystem() = ComptonRestSystem(Energy(2))

# outgoing phase space layout
# spherical coordinates in electron rest frame

@inline function _pert_Compton_omega_prime(pt, cth)
    Et = getE(pt)
    rho2_t = getMag2(pt)
    s = Et^2 - rho2_t

    return (s - 1) / (2 * (Et - sqrt(rho2_t) * cth))
end

struct ComptonSphericalLayout{INPSL<:AbstractTwoBodyInPhaseSpaceLayout} <:
       AbstractOutPhaseSpaceLayout{INPSL}
    in_psl::INPSL
end

function ComptonSphericalLayout(::TwoBodyRestSystem{2})
    throw(
        ArgumentError(
            "the first incoming particle of Compton is an electron, which has no degrees of freedom in the Compton rest frame.",
        ),
    )
end

function QEDbase.phase_space_dimension(
    proc::Compton, model::PerturbativeQED, psl::ComptonSphericalLayout
)
    # cth, phi
    return 2
end

QEDbase.in_phase_space_layout(psl::ComptonSphericalLayout) = psl.in_psl

function QEDbase._build_momenta(
    proc::Compton,
    model::PerturbativeQED,
    psl::ComptonSphericalLayout,
    in_coords::NTuple{1,T},
    out_coords::NTuple{2,T},
) where {T<:Real}
    P, K = QEDbase._build_momenta(proc, model, in_phase_space_layout(psl), in_coords)
    Pt = P + K
    cth, phi = @inbounds out_coords
    omega_prime = _pert_Compton_omega_prime(Pt, cth)
    sth = sqrt(1 - cth^2)
    sphi, cphi = sincos(phi)

    Kp = SFourMomentum(
        omega_prime, omega_prime * sth * cphi, omega_prime * sth * sphi, omega_prime * cth
    )
    Pp = Pt - Kp

    return (P, K), (Pp, Kp)
end
