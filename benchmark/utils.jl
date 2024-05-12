"""
    _find_y_lims(data)

Return the limits which encloses the data, i.e. the supremum and infimum of the data in powers of 10.
This is useful for fining limit of logarithmic data.

"""
function _find_y_lims(data::BenchmarkTools.BenchmarkGroup)
    t_min = time(minimum(data))
    data_min = minimum(values(t_min))
    ymin = 10^floor(log10(data_min))

    t_max = time(maximum(data))
    data_max = maximum(values(t_max))
    ymax = 10^(ceil(log10(data_max)))

    return (ymin, ymax)
end
