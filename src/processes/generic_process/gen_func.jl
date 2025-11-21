"""
    _generic_proc_graph(proc::ScatteringProcess, target::Symbol)

Returns the generated DAG for the given [`ScatteringProcess`](@ref) and the specific target, one of `:mat_el_sqsum`, `:diff_prob`, `:diff_cs`.

!!! note
    This function is memoized so it will cache the result for a unique set of arguments and not reevaluate.
"""
@memoize function _generic_proc_graph(proc::PROC, target::Symbol) where {PROC <: ScatteringProcess}
    g = graph(proc; target = target)
    optimize_to_fixpoint!(ReductionOptimizer(), g)
    return g
end

"""
    _mat_el_sq_sum_func(proc::ScatteringProcess, ::Type{PSP})

Returns the generated compute function ready to be called on a `PhaseSpacePoint::PSP`, returning the square sum of matrix elements.

!!! note
    This function is memoized so it will cache the result for a unique set of arguments and not reevaluate.
"""
@memoize function _mat_el_sq_sum_func(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc, :mat_el_sqsum)
    return compute_function(g, proc, cpu_st(), @__MODULE__)
end

"""
    _diff_cs_kernel(proc::ScatteringProcess, ::Type{PSP})

Returns a generated KernelAbstractions.jl kernel ready to be called to compute differential cross sections on any of KernelAbstractions.jl's backends.
The function signature is `diff_cs(out::AbstractVector{FLOAT_T}, in::AbstractVector{PhaseSpacePoint{...}})`, where the phase space point must be of type `PSP` and its underlying float type must be convertible to `FLOAT_T`.

!!! note
    This function is memoized so it will cache the result for a unique set of arguments and not reevaluate.

See also: [`_diff_prob_kernel`](@ref)
"""
@memoize function _diff_cs_kernel(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc, :diff_cs)
    return kernel(g, proc, @__MODULE__)
end

"""
    _diff_prob_kernel(proc::ScatteringProcess, ::Type{PSP})

Returns a generated KernelAbstractions.jl kernel ready to be called to compute differential probability on any of KernelAbstractions.jl's backends.
The function signature is `diff_cs(out::AbstractVector{FLOAT_T}, in::AbstractVector{PhaseSpacePoint{...}})`, where the phase space point must be of type `PSP` and its underlying float type must be convertible to `FLOAT_T`.

!!! note
    This function is memoized so it will cache the result for a unique set of arguments and not reevaluate.

See also: [`_diff_cs_kernel`](@ref)
"""
@memoize function _diff_prob_kernel(proc::PROC, ::Type{PSP}) where {
        PROC <: ScatteringProcess, PSP <: PhaseSpacePoint{PROC, PerturbativeQED},
    }
    g = _generic_proc_graph(proc, :diff_prob)
    return kernel(g, proc, @__MODULE__)
end
