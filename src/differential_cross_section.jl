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
        is_onshell(in_phase_space[1], mass(incoming_particles(proc)[1])) ||
            throw(InvalidInputError("Should be Photon"))
        is_onshell(in_phase_space[2], mass(incoming_particles(proc)[2]))
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
    _differential_cross_section(model::AbstractModelDefinition, process::AbstractScatteringProcess, phaseSpace::AbstractArray{NumericType}, input::AbstractArray{NumericType}) where {NumericType}

Interface function that must be implemented for valid model-process combinations. This function will be dispatched to when [`_compute`](@ref) is called on a [`DifferentialCrossSection`](@ref), which in turn is called in [`compute`](@ref).

For details on what it should do, see [`_compute`](@ref).
"""
function _differential_cross_section end

"""
    _post_process(model::AbstractModelDefinition, process::AbstractScatteringProcess, phaseSpace::AbstractArray{NumericType}, input::AbstractArray{NumericType}, result::Any) where {NumericType}

Interface function that must be implemented for valid model-process combinations. This function will be dispatched to when [`_post_process`](@ref) is called on a [`DifferentialCrossSection`](@ref), which in turn is called in [`compute`](@ref).

For details on what it should do, see [`_post_process`](@ref).
"""
function _post_process end

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
    isa(dcs.in_phase_space, AbstractArray{NumericType}) || throw(
        InvalidInputError(
            "Input type is '$(typeof(out_phase_space))' but '$(typeof(dcs.in_phase_space))' was expected",
        ),
    )
    length(dcs.in_phase_space) == number_incoming_particles(dcs.process) || throw(
        InvalidInputError(
            "In phase space dimension is inconsistent with input size ($(size(dcs.in_phase_space)) versus $(number_incoming_particles(dcs.process)))",
        ),
    )
    length(out_phase_space) == number_outgoing_particles(dcs.process) || throw(
        InvalidInputError(
            "Out phase space dimension is inconsistent with input size ($(size(out_phase_space, 1)) versus $(number_outgoing_particles(dcs.process)))",
        ),
    )
    _assert_valid_input(dcs.process, dcs.model, dcs.in_phase_space, out_phase_space)
end

function _compute(dcs::DifferentialCrossSection, input)
    return _differential_cross_section(dcs.process, dcs.model, dcs.in_phase_space, input)
end

function _post_processing(dcs::DifferentialCrossSection, input, result)
    return _post_process(dcs.process, dcs.model, dcs.in_phase_space, input, result)
end
