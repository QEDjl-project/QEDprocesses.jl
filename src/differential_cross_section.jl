###############
# DifferentialCrossSection
#
# In this file, we define the DifferentialCrossSection type and 
# implement the AbstractProcessSetup interface to provide its functionality.
###############

struct DifferentialCrossSection{Model,Process,PhaseSpace} <: AbstractProcessSetup
    model::Model
    process::Process
    in_phase_space::PhaseSpace

    function DifferentialCrossSection(
        proc::Process,
        model::Model,
        in_phase_space::PhaseSpace,
    ) where {
        Model<:AbstractModelDefinition,
        Process<:AbstractScatteringProcess,
        NumericType,
        PhaseSpace<:AbstractArray{NumericType},
    }
        is_onshell(in_phase_space[1], mass(incoming_particles(proc)[1])) || throw(
            InvalidInputError(
                "First in-phasespace particle must be an on-shell Photon\nValue: $(in_phase_space[1])",
            ),
        )
        is_onshell(in_phase_space[2], mass(incoming_particles(proc)[2])) || throw(
            InvalidInputError(
                "Second in-phasespace particle must be an on-shell Electron\nValue: $(in_phase_space[1])",
            ),
        )

        length(in_phase_space) == number_incoming_particles(proc) || throw(
            InvalidInputError(
                "In-phasespace dimension is inconsistent with input size ($(length(in_phase_space)) versus $(number_incoming_particles(proc)))",
            ),
        )
        return new{Model,Process,PhaseSpace}(model, proc, in_phase_space)
    end
end

"""
    _assert_valid_input(model::AbstractModelDefinition, process::AbstractScatteringProcess, phaseSpace::AbstractArray{NumericType}, input::AbstractArray{NumericType}) where {NumericType}

Interface function that must be implemented for valid model-process combinations. This function will be dispatched to when [`_assert_valid_input`](@ref) is called on a [`DifferentialCrossSection`](@ref), which in turn is called in [`compute`](@ref).

For details on what it should do, see [`_assert_valid_input`](@ref).
"""
function _assert_valid_dcs_input end

"""
    _post_process_dcs(model::AbstractModelDefinition, process::AbstractScatteringProcess, phaseSpace::AbstractArray{NumericType}, input::AbstractArray{NumericType}, result::Any) where {NumericType}

Interface function that must be implemented for valid model-process combinations. This function will be dispatched to when [`_post_processing`](@ref) is called on a [`DifferentialCrossSection`](@ref), which in turn is called in [`compute`](@ref).

For details on what it should do, see [`_post_processing`](@ref).
"""
function _post_process_dcs end

function scattering_process(dcs::DifferentialCrossSection)
    return dcs.process
end

function physical_model(dcs::DifferentialCrossSection)
    return dcs.model
end

function _assert_valid_input(
    dcs::DifferentialCrossSection,
    out_phase_space::AbstractArray{NumericType},
) where {NumericType}
    length(out_phase_space) == number_outgoing_particles(dcs.process) || throw(
        InvalidInputError(
            "Out-phasespace dimension is inconsistent with input size ($(length(out_phase_space)) versus $(number_outgoing_particles(dcs.process)))",
        ),
    )
    _assert_valid_dcs_input(dcs.process, dcs.model, dcs.in_phase_space, out_phase_space)
end

function _compute(dcs::DifferentialCrossSection, input)
    return _differential_cross_section(dcs.process, dcs.model, dcs.in_phase_space, input)
end

function _post_processing(dcs::DifferentialCrossSection, input, result)
    return _post_process_dcs(dcs.process, dcs.model, dcs.in_phase_space, input, result)
end
