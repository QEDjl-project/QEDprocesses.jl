
function _assert_valid_input(
    ::Compton,
    ::PerturbativeQED,
    in_phase_space::AbstractArray{NumericType},
    out_phase_space::AbstractArray{NumericType},
) where {NumericType<:QEDbase.AbstractFourMomentum}
    @assert is_onshell(out_phase_space[1])
    @assert is_onshell(out_phase_space[1], 1.0)
    return nothing
end


function _differential_cross_section(
    process::Compton,
    model::PerturbativeQED,
    in_phase_space::AbstractVector{NumericType},
    out_phase_space::AbstractVector{NumericType},
)::Float64 where {NumericType<:QEDbase.AbstractFourMomentum}
    if (!isapprox(sum(in_phase_space), sum(out_phase_space); rtol = sqrt(eps())))
        return zero(Float64)
    end

    k = in_phase_space[1]
    kp = out_phase_space[1]
    p = in_phase_space[2]
    pp = out_phase_space[2]

    #=
    in_states = Vector{Tuple{AbstractParticle,NumericType}}()
    for (particle_type, in_state) in zip(incoming_particles(process), in_phase_space)
        push!(in_states, incoming_particle_state(particle_type, spin_or_pol(process, particle_type)))
    end
    =#

    p_state = incoming_particle_state(incoming_particles(process)[1], spins(process, 1), p)
    pp_state = outgoing_electron_state(spins(process, 2), pp)
    k_state = incoming_photon_state(pol(process, 1), k)
    kp_state = outgoing_photon_state(pol(process, 2), kp)

    matrix_elements = Vector{ComplexF64}()
    iterator = Iterators.product(p_state, pp_state, k_state, kp_state)
    sizehint!(matrix_elements, length(iterator))
    for (p1, p2, p3, p4) in iterator
        push!(matrix_elements, _perturbative_compton_matrix(p1, p2, p3, p4))
    end

    matrix_elements_sq = abs2.(matrix_elements)

    normalization = 1.0 / (length(p_state) * length(k_state))
    I = k * p

    return normalization * sum(matrix_elements_sq) * _phase_space_factor(p, pp, k, kp)
end

function _post_process(
    process::Compton,
    model::PerturbativeQED,
    in_phase_space::AbstractArray{NumericType},
    out_phase_space::AbstractArray{NumericType},
    result::Real,
) where {NumericType<:QEDbase.AbstractFourMomentum}
    # generally nothing to be done here
    return result
end
