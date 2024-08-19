# None of these functions are currently implemented.
# They will be implemented using the QEDFeynmanDiagrams project when that is released.

#=
function QEDbase._incident_flux(psp::InPhaseSpacePoint{<:GenericQEDProcess,PerturbativeQED}) end

function QEDbase._matrix_element(psp::PhaseSpacePoint{<:GenericQEDProcess,PerturbativeQED})
    proc = process(psp)
    # by using FeynmanDiagramGenerator.jl
    return sqrt(proc.matrix_element_squared(psp))
end

function QEDbase._averaging_norm(proc::<:GenericQEDProcess) end

function QEDbase._is_in_phasespace(psp::PhaseSpacePoint{<:GenericQEDProcess,PerturbativeQED}) end

function QEDbase._phase_space_factor(
    psp::PhaseSpacePoint{<:GenericQEDProcess,PerturbativeQED,CustomPhasespaceDefinition}
) end

function QEDbase._total_probability(in_psp::InPhaseSpacePoint{<:GenericQEDProcess,PerturbativeQED}) end
=#
