using QEDcore, QEDprocesses
using Random

"""
Return a tuple of tuples of incoming and outgoing coordinates for a given process, model and ps_def that make up a physical phase space point.
"""
function _rand_coordinates(
        rng::AbstractRNG, ::PROCESS, ::MODEL, ::PSL, FLOAT_T = Float64
    ) where {PROCESS <: Compton, MODEL <: PerturbativeQED, PSL <: AbstractPhaseSpaceLayout}
    # TODO: this adds a small bit to the incoming energy to prevent precision problems that appear when it's close to 0
    return ((FLOAT_T(0.95) * rand(rng, FLOAT_T) + FLOAT_T(0.05),), (rand(rng, FLOAT_T), rand(rng, FLOAT_T)))
end

tuple_isapprox(::Tuple{}, ::Tuple{Vararg}; atol = 0.0, rtol = eps()) = false
tuple_isapprox(::Tuple{Vararg}, ::Tuple{}; atol = 0.0, rtol = eps()) = false
tuple_isapprox(::Tuple{}, ::Tuple{}; atol = 0.0, rtol = eps()) = true
function tuple_isapprox(
        a::Tuple{<:Number, Vararg}, b::Tuple{<:Number, Vararg}; atol = 0.0, rtol = eps()
    )
    return isapprox(a[1], b[1]; atol = atol, rtol = rtol) &&
        tuple_isapprox(a[2:end], b[2:end]; atol = atol, rtol = rtol)
end

# generate random on-shell photon
function _rand_mom(rng::AbstractRNG, ::Photon)
    w = rand(rng)
    return SFourMomentum(w, 0.0, 0.0, w)
end

# generate random on-shell particle with mass
function _rand_mom(rng::AbstractRNG, pt::AbstractParticleType)
    (x, y, z) = (rand(rng), rand(rng), rand(rng))
    E2 = x^2 + y^2 + z^2 + mass(pt)^2
    return SFourMomentum(sqrt(E2), x, y, z)
end

# for a process definition, generate a random phase space point for the given process
# once QEDevents.jl supports random phase space point generation, that should be used here instead
# since it generally generates one off-shell particle, it should not be used for actual calculations
function _rand_psp(rng::AbstractRNG, proc::AbstractProcessDefinition, model::AbstractModelDefinition, psl::AbstractPhaseSpaceLayout)
    in_momenta = ((_rand_mom(rng, p) for p in incoming_particles(proc))...,)
    out_momenta = ((_rand_mom(rng, p) for p in outgoing_particles(proc))...,)

    leftover = sum(in_momenta) - sum(out_momenta)

    local fermion_like_index = -1
    for index in eachindex(in_momenta)
        if (incoming_particles(proc)[index] isa FermionLike)
            fermion_like_index = index
            break
        end
    end

    if (fermion_like_index != -1)
        # this leads to the phase space point conserving momentum
        # it also makes the FermionLike off-shell
        in_momenta = ntuple(
            i -> i == fermion_like_index ? in_momenta[i] + leftover : in_momenta[i],
            length(in_momenta),
        )
    else
        # try out_momenta (there has to be a fermionlike somewhere, otherwise the process would not be valid)
        for index in eachindex(out_momenta)
            if (outgoing_particles(proc)[index] isa FermionLike)
                fermion_like_index = index
                break
            end
        end
        out_momenta = ntuple(
            i -> i == fermion_like_index ? out_momenta[i] + leftover : out_momenta[i],
            length(out_momenta),
        )
    end

    return PhaseSpacePoint(proc, model, psl, in_momenta, out_momenta)
end
