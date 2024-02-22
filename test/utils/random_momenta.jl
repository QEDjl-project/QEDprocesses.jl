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
Return a random phase space point, which is failing the phase space constraint, 
i.e. the first entry of the vector is the null momentum.
"""
function _rand_momenta_failed(rng::AbstractRNG, N)
    moms = _rand_momenta(rng, N)
    moms[1] = SFourMomentum(zeros(4))
    return moms
end

"""
Return a collection of phase space points, where the first point is failing the phase space constraint, 
i.e. the first entry of the matrix is the null momentum, but the others pass. 
"""
function _rand_momenta_failed_mix(rng::AbstractRNG, N1, N2)
    moms = _rand_momenta(rng, N1, N2)
    moms[1, 1] = SFourMomentum(zeros(4))
    return moms
end

"""
Return a collection of phase space points, where all points are failing the phase space constraint, 
i.e. their first entries are null momenta.
"""
function _rand_momenta_failed_all(rng::AbstractRNG, N1, N2)
    moms = _rand_momenta(rng, N1, N2)
    moms[1, 1] = SFourMomentum(zeros(4))
    moms[1, 2] = SFourMomentum(zeros(4))
    return moms
end
