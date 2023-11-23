include("random_momenta.jl")
include("groundtruths.jl")
include("test_implementations.jl")



# Check if any failed type is containted
_any_fail(x...) = true
_any_fail(::TestProcess, ::TestModel) = false
_any_fail(::TestProcess, ::TestModel, ::TestPhasespaceDef, ::TestPhasespaceDef) = false

