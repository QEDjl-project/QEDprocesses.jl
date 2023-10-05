###############
# DifferentialCrossSection
#
# In this file, we define the DifferentialCrossSection type and 
# implement the AbstractProcessSetup interface to provide its functionality.
###############

struct DifferentialCrossSection{Model,Process,PhaseSpace,NumericType} <:
       AbstractProcessSetup where {Model<:AbstractModelDefinition,Process<:AbstractScatteringProcess,PhaseSpace<:AbstractArray{NumericType}}
    model::Model
    process::Process
    inPhaseSpace::PhaseSpace

    function DifferentialCrossSection(
        model::Model,
        proc::Process,
        inPS::PhaseSpace,
    ) where {Model,Process,PhaseSpace<:AbstractArray{NumericType},NumericType}
        return new{Model,Process,PhaseSpace,NumericType}(model, proc, inPS)
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

function _assert_valid_input(dcs::DCS, input::AbstractArray{NumericType}) where {NumericType}
    isa(dcs.inPhaseSpace, AbstractArray{NumericType}) || throw(
        InvalidInputError("Input type is '$(typeof(input))' but '$(typeof(dcs.inPS))' was expected.")
    )
    size(inPS, 1) == in_phasespace_dimension(proc, model) || throw(
        InvalidInputError("In-phase-space dimension is inconsistent with input size ($(size(inPS, 1)) versus $(in_phasespace_dimension(proc, model)))."),
    )
    _assert_valid_input(dcs.model, dcs.process, dcs.inPhaseSpace, input)
end

function _compute(dcs::DifferentialCrossSection, input)
    _differential_cross_section(dcs.model, dcs.process, dcs.inPhaseSpace, input)
end

function _post_processing(dcs::DifferentialCrossSection, input, result)
    return _post_process(dcs.model, dcs.process, dcs.inPhaseSpace, input, result)
end
