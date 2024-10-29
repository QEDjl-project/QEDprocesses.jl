function QEDbase._total_probability(in_psp::InPhaseSpacePoint{<:Compton,PerturbativeQED})
    omega = getE(momentum(in_psp[Incoming(), 2]))

    function func(x)
        return unsafe_differential_probability(
            PhaseSpacePoint(in_psp.proc, in_psp.model, in_psp.ps_def, (omega,), (x, 0.0))
        )
    end

    tot_prob, _ = quadgk(func, -1, 1; rtol=sqrt(eps(omega)))

    tot_prob *= 2 * pi # phi integration is trivial
    return tot_prob
end
