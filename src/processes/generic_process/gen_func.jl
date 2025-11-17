using Memoization

@memoize function _generic_proc_graph(proc::PROC, target::Symbol) where {PROC <: ScatteringProcess}
    g = graph(proc; target = target)
    optimize_to_fixpoint!(ReductionOptimizer(), g)
    return g
end

@memoize function _mat_el_func(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc, :mat_el_sqsum)
    return compute_function(g, proc, cpu_st(), @__MODULE__)
end

@memoize function _mat_el_kernel(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc, :mat_el_sqsum)
    return kernel(g, proc, @__MODULE__)
end

@memoize function _diff_cs_kernel(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc, :diff_cs)
    return kernel(g, proc, @__MODULE__)
end

@memoize function _diff_prob_kernel(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc, :diff_prob)
    return kernel(g, proc, @__MODULE__)
end
