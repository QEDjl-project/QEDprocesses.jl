
"""
Return a vector of random four momenta, i.e. a random phase space point
"""
function _rand_momenta(rng::AbstractRNG, N)
    moms = Vector{SFourMomentum}(undef, N)
    for i in 1:N
        moms[i] = SFourMomentum(rand(rng, 4))
    end
    return moms
end
"""
Return a matrix of random four momenta, i.e. a collection of phase space points
"""
function _rand_momenta(rng::AbstractRNG, N1, N2)
    moms = Matrix{SFourMomentum}(undef, N1, N2)
    for i in 1:N1
        for j in 1:N2
            moms[i, j] = SFourMomentum(rand(rng, 4))
        end
    end
    return moms
end

"""
Return a random phase space point, which is failing the incoming phase space constraint, 
i.e. the first entry of the vector is the null momentum.
"""
function _rand_in_momenta_failed(rng::AbstractRNG, N)
    moms = _rand_momenta(rng, N)
    moms[1] = zero(SFourMomentum)
    return moms
end

"""
Return a random phase space point, which is failing the outgoing phase space constraint, 
i.e. the last entry of the vector is the unit momentum.
"""
function _rand_out_momenta_failed(rng::AbstractRNG, N)
    moms = _rand_momenta(rng, N)
    moms[end] = ones(SFourMomentum)
    return moms
end

"""
Return a collection of incoming phase space points, where the first point is failing the phase space constraint, 
i.e. the first entry of the matrix is the null momentum, but the others pass. 
"""
function _rand_in_momenta_failed_mix(rng::AbstractRNG, N1, N2)
    moms = _rand_momenta(rng, N1, N2)
    moms[1, 1] = zero(SFourMomentum)
    return moms
end

"""
Return a collection of incoming phase space points, where all points are failing the phase space constraint, 
i.e. their first entries are null momenta.
"""
function _rand_in_momenta_failed_all(rng::AbstractRNG, N1, N2)
    moms = _rand_momenta(rng, N1, N2)
    for n in 1:N2
        moms[1, n] = zero(SFourMomentum)
    end
    return moms
end

"""
Return a collection of outgoing phase space points, where the first point is failing the phase space constraint, 
i.e. the last entry of the matrix is the unit momentum, but the others pass. 
"""
function _rand_out_momenta_failed_mix(rng::AbstractRNG, N1, N2)
    moms = _rand_momenta(rng, N1, N2)
    moms[end, 1] = ones(SFourMomentum)
    return moms
end

"""
Return a collection of outgoing phase space points, where all points are failing the phase space constraint, 
i.e. their last entries are unit momenta.
"""
function _rand_out_momenta_failed_all(rng::AbstractRNG, N1, N2)
    moms = _rand_momenta(rng, N1, N2)
    for n in 1:N2
        moms[end, n] = ones(SFourMomentum)
    end
    return moms
end
