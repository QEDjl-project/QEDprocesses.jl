using Memoization

@memoize function _generic_proc_graph(proc::PROC) where {PROC <: ScatteringProcess}
    g = graph(proc)
    optimize_to_fixpoint!(ReductionOptimizer(), g)
    return g
end

@memoize function _mat_el_func(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc)
    return compute_function(g, proc, cpu_st(), @__MODULE__)
end

@memoize function _mat_el_kernel(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc)
    return kernel(g, proc, @__MODULE__)
end
