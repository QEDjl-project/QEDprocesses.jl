
# Check if any failed type is in the input 
_any_fail(x...) = true
_any_fail(::TestProcess, ::TestModel) = false
_any_fail(::TestProcess, ::TestModel, ::TestPhasespaceDef, ::TestPhasespaceDef) = false
