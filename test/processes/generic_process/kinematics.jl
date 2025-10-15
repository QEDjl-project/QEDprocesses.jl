using QEDprocesses
using Random
using QEDcore

include("groundtruths.jl")
include("../../utils.jl")
include("../../test_implementation/random_coordinates.jl")

const RNG = MersenneTwister(1597463007)

MODEL = PerturbativeQED()
PARTICLES = [Electron(), Positron(), Photon()]
PSL = FlatPhaseSpaceLayout(TwoBodyRestSystem())

@testset "cross section interface ($(FLOAT_T))" for FLOAT_T in (Float64, Float32)
    n_out_particles = LARGE_TESTS() ? (2, 3, 4) : (2, 3)
    @testset "$n -> $m processes" for (n, m) in Base.product((2,), n_out_particles)
        IN_PARTICLES = Tuple(rand(RNG, PARTICLES, n))
        OUT_PARTICLES = Tuple(rand(RNG, PARTICLES, m))

        while (!isphysical(ScatteringProcess(IN_PARTICLES, OUT_PARTICLES), MODEL))
            IN_PARTICLES = Tuple(rand(RNG, PARTICLES, n))
            OUT_PARTICLES = Tuple(rand(RNG, PARTICLES, m))
        end

        @testset "spin/pols: $in_sp -> $out_sp" for (in_sp, out_sp) in Iterators.zip(_random_spin_pols(RNG, IN_PARTICLES, 5), _random_spin_pols(RNG, OUT_PARTICLES, 5))
            proc = ScatteringProcess(IN_PARTICLES, OUT_PARTICLES, in_sp, out_sp)
            psps = [_rand_psp(RNG, proc, MODEL, PSL) for _ in 1:100]

            differential_cross_section.(psps)
            differential_probability.(psps)

            # integration is not implemented yet
            #total_cross_section.(psps)
            #total_probability.(psps)
        end
    end
end
