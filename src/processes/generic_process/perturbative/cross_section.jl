# None of these functions are currently implemented.
# They will be implemented using the QEDFeynmanDiagrams project when that is released.

#=
function QEDbase._incident_flux(psp::InPhaseSpacePoint{<:ScatteringProcess,PerturbativeQED}) end

function QEDbase._matrix_element_squared(psp::PhaseSpacePoint{<:ScatteringProcess,PerturbativeQED})
    proc = process(psp)
    # by using FeynmanDiagramGenerator.jl
    return proc.matrix_element_squared(psp)
end

function QEDbase._averaging_norm(proc::<:ScatteringProcess) end

function QEDbase._is_in_phasespace(psp::PhaseSpacePoint{<:ScatteringProcess,PerturbativeQED}) end

function QEDbase._phase_space_factor(
    psp::PhaseSpacePoint{<:ScatteringProcess,PerturbativeQED,CustomPhasespaceDefinition}
) end
=#
