using Random
import QEDbase
using QEDcore
using QEDprocesses

RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

function _rand_momentum(rng::AbstractRNG)
    return SFourMomentum(rand(rng, 4))
end

groundtruth_propagator(::QEDbase.Photon, mom) = one(eltype(mom)) / (mom * mom)
function groundtruth_propagator(particle::QEDbase.FermionLike, mom)
    return (slashed(mom) + QEDbase.mass(particle) * one(DiracMatrix)) /
           (mom * mom - QEDbase.mass(particle)^2)
end

@testset "propagators" begin
    @testset "$P" for P in (QEDbase.Electron(), QEDbase.Positron(), QEDbase.Photon())
        mom = _rand_momentum(RNG)
        groundtruth = groundtruth_propagator(P, mom)
        test_prop = propagator(P, mom)
        @test isapprox(test_prop, groundtruth, atol=ATOL, rtol=RTOL)
    end
end
