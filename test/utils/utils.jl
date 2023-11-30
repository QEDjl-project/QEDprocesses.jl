include("random_momenta.jl")
include("groundtruths.jl")
include("test_implementations.jl")

# Check if any failed type is containted
_any_fail(x...) = true
_any_fail(::TestProcess, ::TestModel) = false
_any_fail(::TestProcess, ::TestModel, ::TestPhasespaceDef, ::TestPhasespaceDef) = false

# unrolls all elements of a list of four-momenta into vector of coordinates
function _unroll_moms(ps_moms::AbstractVector{T}) where {T<:QEDbase.AbstractFourMomentum}
     collect(Iterators.flatten(ps_moms))
end

function _unroll_moms(ps_moms::AbstractMatrix{T}) where {T<:QEDbase.AbstractFourMomentum}
end


# collect components of four-momenta from a vector of coordinates
function __furl_moms(ps_coords::AbstractVector{T}) where {T<:Real}
    SFourMomentum.(eachcol(reshape(ps_coords,4,:))) 
end

function _furl_moms(ps_coords::AbstractVector{T}) where {T<:Real}
    @assert length(ps_coords)%4==0
    __furl_moms(ps_coords)
end

function _furl_moms(ps_coords::AbstractMatrix{T}) where {T<:Real}
    @assert size(ps_coords,1)%4==0
    res = Matrix{SFourMomentum}(undef,Int(size(ps_coords,1)//4),size(ps_coords,2))
    for i in 1:size(ps_coords,2)
        res[:,i] .= __furl_moms(view(ps_coords,:,i))
    end
    return res
end
