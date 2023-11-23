
function _groundtruth_incident_flux(in_ps)
    s = sum(in_ps)
    return s * s
end

function _groundtruth_matrix_element(in_ps, out_ps)
    s_in = sum(in_ps)
    s_out = sum(out_ps)
    res = s_in * s_in + 1im * (s_out * s_out)
    return (res, 2 * res, 3 * res)
end

function _groundtruth_averaging_norm(proc)
    return 1.0 / (number_incoming_particles(proc) + number_outgoing_particles(proc))
end
function _groundtruth_is_in_phasespace(in_ps, out_ps)
    if in_ps[1] == SFourMomentum(zeros(4))
        return false
    end
    return true
end

function _groundtruth_phase_space_factor(in_ps, out_ps)
    en_in = getE.(in_ps)
    en_out = getE.(out_ps)
    return 1 / (prod(en_in) * prod(en_out))
end

function _groundtruth_unsafe_probability(proc, in_ps, out_ps)
    mat_el = _groundtruth_matrix_element(in_ps, out_ps)
    mat_el_sq = abs2.(mat_el)
    normalization = _groundtruth_averaging_norm(proc)
    ps_fac = _groundtruth_phase_space_factor(in_ps, out_ps)
    return sum(mat_el_sq) * ps_fac * normalization
end

function _groundtruth_unsafe_probability(
    proc, in_ps::AbstractVector, out_ps::AbstractMatrix
)
    res = Vector{Float64}(undef, size(out_ps, 2))
    for i in 1:size(out_ps, 2)
        res[i] = _groundtruth_unsafe_probability(proc, in_ps, view(out_ps, :, i))
    end
    return res
end

function _groundtruth_unsafe_probability(
    proc, in_ps::AbstractMatrix, out_ps::AbstractVector
)
    res = Vector{Float64}(undef, size(in_ps, 2))
    for i in 1:size(in_ps, 2)
        res[i] = _groundtruth_unsafe_probability(proc, view(in_ps, :, i), out_ps)
    end
    return res
end

function _groundtruth_unsafe_probability(
    proc, in_ps::AbstractMatrix, out_ps::AbstractMatrix
)
    res = Matrix{Float64}(undef, size(in_ps, 2), size(out_ps, 2))
    for i in 1:size(in_ps, 2)
        for j in 1:size(out_ps, 2)
            res[i, j] = _groundtruth_unsafe_probability(
                proc, view(in_ps, :, i), view(out_ps, :, j)
            )
        end
    end
    return res
end

function _groundtruth_safe_probability(proc, in_ps, out_ps)
    if !_groundtruth_is_in_phasespace(in_ps, out_ps)
        return zero(Float64)
    end
    return _groundtruth_unsafe_probability(proc, in_ps, out_ps)
end

function _groundtruth_safe_probability(proc, in_ps::AbstractVector, out_ps::AbstractMatrix)
    res = Vector{Float64}(undef, size(out_ps, 2))
    for i in 1:size(out_ps, 2)
        res[i] = _groundtruth_safe_probability(proc, in_ps, view(out_ps, :, i))
    end
    return res
end

function _groundtruth_safe_probability(proc, in_ps::AbstractMatrix, out_ps::AbstractVector)
    res = Vector{Float64}(undef, size(in_ps, 2))
    for i in 1:size(in_ps, 2)
        res[i] = _groundtruth_safe_probability(proc, view(in_ps, :, i), out_ps)
    end
    return res
end

function _groundtruth_safe_probability(proc, in_ps::AbstractMatrix, out_ps::AbstractMatrix)
    res = Matrix{Float64}(undef, size(in_ps, 2), size(out_ps, 2))
    for i in 1:size(in_ps, 2)
        for j in 1:size(out_ps, 2)
            res[i, j] = _groundtruth_safe_probability(
                proc, view(in_ps, :, i), view(out_ps, :, j)
            )
        end
    end
    return res
end

function _groundtruth_unsafe_diffCS(proc, in_ps, out_ps)
    init_flux = _groundtruth_incident_flux(in_ps)
    return _groundtruth_unsafe_probability(proc, in_ps, out_ps) / (4 * init_flux)
end

function _groundtruth_unsafe_diffCS(proc, in_ps::AbstractVector, out_ps::AbstractMatrix)
    res = Vector{Float64}(undef, size(out_ps, 2))
    for i in 1:size(out_ps, 2)
        res[i] = _groundtruth_unsafe_diffCS(proc, in_ps, view(out_ps, :, i))
    end
    return res
end

function _groundtruth_unsafe_diffCS(proc, in_ps::AbstractMatrix, out_ps::AbstractVector)
    res = Vector{Float64}(undef, size(in_ps, 2))
    for i in 1:size(in_ps, 2)
        res[i] = _groundtruth_unsafe_diffCS(proc, view(in_ps, :, i), out_ps)
    end
    return res
end

function _groundtruth_unsafe_diffCS(proc, in_ps::AbstractMatrix, out_ps::AbstractMatrix)
    res = Matrix{Float64}(undef, size(in_ps, 2), size(out_ps, 2))
    for i in 1:size(in_ps, 2)
        for j in 1:size(out_ps, 2)
            res[i, j] = _groundtruth_unsafe_diffCS(
                proc, view(in_ps, :, i), view(out_ps, :, j)
            )
        end
    end
    return res
end

function _groundtruth_safe_diffCS(proc, in_ps, out_ps)
    if !_groundtruth_is_in_phasespace(in_ps, out_ps)
        return zero(Float64)
    end
    return _groundtruth_unsafe_diffCS(proc, in_ps, out_ps)
end

function _groundtruth_safe_diffCS(proc, in_ps::AbstractVector, out_ps::AbstractMatrix)
    res = Vector{Float64}(undef, size(out_ps, 2))
    for i in 1:size(out_ps, 2)
        res[i] = _groundtruth_safe_diffCS(proc, in_ps, view(out_ps, :, i))
    end
    return res
end

function _groundtruth_safe_diffCS(proc, in_ps::AbstractMatrix, out_ps::AbstractVector)
    res = Vector{Float64}(undef, size(in_ps, 2))
    for i in 1:size(in_ps, 2)
        res[i] = _groundtruth_safe_diffCS(proc, view(in_ps, :, i), out_ps)
    end
    return res
end

function _groundtruth_safe_diffCS(proc, in_ps::AbstractMatrix, out_ps::AbstractMatrix)
    res = Matrix{Float64}(undef, size(in_ps, 2), size(out_ps, 2))
    for i in 1:size(in_ps, 2)
        for j in 1:size(out_ps, 2)
            res[i, j] = _groundtruth_safe_diffCS(
                proc, view(in_ps, :, i), view(out_ps, :, j)
            )
        end
    end
    return res
end
