
function _total_probability(
    proc::Compton,
    model::PerturbativeQED,
    in_phase_space_def::AbstractPhasespaceDefinition,
    in_phase_space::AbstractVector{T},
) where {T<:QEDbase.AbstractFourMomentum}
    omega = getE(in_phase_space[2])

    function func(x)
        return unsafe_differential_probability(
            proc, model, in_phase_space_def, [omega], [x, 0.0]
        )
    end

    tot_prob, _ = quadgk(func, -1, 1; rtol=sqrt(eps(omega)))

    tot_prob *= 2 * pi # phi integration is trivial
    return tot_prob
end
