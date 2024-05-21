var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = QEDprocesses","category":"page"},{"location":"#QEDprocesses","page":"Home","title":"QEDprocesses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for QEDprocesses.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [QEDprocesses]","category":"page"},{"location":"#QEDprocesses.AbstractComputationSetup","page":"Home","title":"QEDprocesses.AbstractComputationSetup","text":"Abstract base type for computation setups.  A setup means a collection of setup data needed to evaluate a dedicated quantity of given running data. Therefore, each setup is associated with a single quantity, which one may compute using the setup data and the running data.  Despite that, the decomposition into setup and running data is arbitrary, and this can be used for cases where a subset of the variables a quantity depends on is kept constant. \n\nnote: Computation setup interface\nThe computation performed using a computation setup is separated into three steps:1. input validation\n2. actual computation\n3. post processingwhere every step has its own interface function (see compute for details). Input validationEvery subtype of AbstractComputationSetup should implement the interface function_assert_valid_input(stp::AbstractComputationSetup, input)which should throw and an exception subtyped from AbstractInvalidInputException if the input is not valid for the computation of the associated quantity (see _assert_valid_input for more details).  The default implementation does nothing, i.e. every input is valid by default. Provide a custom implementation if a different behavior is required.Actual computationEvery subtype of AbstractComputationSetup must at least implement the required interface function_compute(stp::AbstractComputationSetup, input) which computes the value of the associated quantity for a given input (see _compute for more details).Post processingEvery subtype of AbstractComputationSetup should implement the interface function_post_processing(stp::AbstractComputationSetup, input, result) which performs task after the actual computation, e.g. conversions or normalizations (see _post_processing for more details).\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.AbstractInvalidInputException","page":"Home","title":"QEDprocesses.AbstractInvalidInputException","text":"Abstract base type for exceptions indicating invalid input. See InvalidInputError for a simple concrete implementation.  Concrete implementations should at least implement \n\n\nBase.showerror(io::IO, err::CustomInvalidError) where {CustomInvalidError<:AbstractInvalidInputException}\n\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.AbstractModelDefinition","page":"Home","title":"QEDprocesses.AbstractModelDefinition","text":"Abstract base type for all compute model definitions in the context of scattering processes. Every subtype of AbstractModelDefinition is associated with a fundamental interaction.  Therefore, one needs to implement the following soft interface function\n\nfundamental_interaction_type(::AbstractModelDefinition)\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.AbstractProcessDefinition","page":"Home","title":"QEDprocesses.AbstractProcessDefinition","text":"Abstract base type for definitions of scattering processes. It is the root type for the  process interface, which assumes that every subtype of AbstractProcessDefinition implements at least \n\nincoming_particles(proc_def::AbstractProcessDefinition)\noutgoing_particles(proc_def::AbstractProcessDefinition)\n\nwhich return a tuple of the incoming and outgoing particles, respectively.\n\nFurthermore, to calculate scattering probabilities and differential cross sections, the following  interface functions need to be implemented for every combination of CustomProcess<:AbstractProcessDefinition,  CustomModel<:AbstractModelDefinition, and CustomPhasespaceDefinition<:AbstractPhasespaceDefinition.\n\n    _incident_flux(\n        proc::AbstractProcessDefinition,\n        model::AbstractModelDefinition, \n        in_phase_space::AbstractVector{T}\n        ) where {T<:QEDbase.AbstractFourMomentum}\n\n    _matrix_element(psp::PhaseSpacePoint{CustomProcess,CustomModel})\n\n    _averaging_norm(proc::CustomProcess)\n\n    _is_in_phasespace(psp::PhaseSpacePoint{CustomProcess,CustomModel})\n\n    _phase_space_factor(psp::PhaseSpacePoint{CustomProcess,CustomModel,CustomPhasespaceDefinition})\n\nOptional is the implementation of \n\n\n    _total_probability(psp::PhaseSpacePoint{CustomProcess,CustomModel,CustomPhasespaceDefinition})\n\n\nto enable the calculation of total probabilities and cross sections.\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.AbstractProcessSetup","page":"Home","title":"QEDprocesses.AbstractProcessSetup","text":"Abstract base type for setups related to combining scattering processes and physical models.   Every subtype of AbstractProcessSetup must implement at least the following  interface functions:\n\nscattering_process(::AbstractProcessSetup) \nphysical_model(::AbstractProcessSetup) \n\nDerived from these interface functions, the following delegations are provided:\n\nnumber_incoming_particles(::AbstractProcessSetup)\nnumber_outgoing_particles(::AbstractProcessSetup)\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.Compton","page":"Home","title":"QEDprocesses.Compton","text":"Compton(\n    in_spin [= AllSpin()]\n    in_pol [= AllPol()]\n    out_spin [= AllSpin()]\n    out_pol [= AllPol()]\n)\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.InvalidInputError","page":"Home","title":"QEDprocesses.InvalidInputError","text":"InvalidInputError(msg::String)\n\nException which is thrown if a given input is invalid, e.g. passed to _assert_valid_input.\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.ParticleStateful","page":"Home","title":"QEDprocesses.ParticleStateful","text":"ParticleStateful <: AbstractParticle\n\nRepresentation of a particle with a state. It has four fields:\n\ndir::ParticleDirection: The direction of the particle, QEDbase.Incoming() or QEDbase.Outgoing().\nspecies::AbstractParticleType: The species of the particle, QEDbase.Electron(), QEDbase.Positron() etc.\nmom::AbstractFourMomentum: The momentum of the particle.\nspin_or_pol::AbstractSpinOrPolarization: The spin or polarization of the particle, QEDbase.SpinUp(), QEDbase.PolX() etc. Can only use spins with fermions and polarizations with bosons.AllSpinorAllPol` by default.\n\nOverloads for QEDbase.is_fermion, QEDbase.is_boson, QEDbase.is_particle, QEDbase.is_anti_particle, QEDbase.is_incoming, QEDbase.is_outgoing, QEDbase.mass, and QEDbase.charge are provided, delegating the call to the correct field and thus implementing the QEDbase.AbstractParticle interface.\n\nThe legality of the combination of species and spin_or_pol is checked on construction. If, for example, the construction of an Electron() with a polarization is attempted, an InvalidInputError is thrown.\n\njulia> using QEDbase; using QEDprocesses\n\njulia> ParticleStateful(Incoming(), Electron(), SFourMomentum(1, 0, 0, 0))\nParticleStateful: incoming electron\n    spin: all spins\n    momentum: [1.0, 0.0, 0.0, 0.0]\n\njulia> ParticleStateful(Outgoing(), Photon(), SFourMomentum(1, 0, 0, 0), PolX())\nParticleStateful: outgoing photon\n    polarization: x-polarized\n    momentum: [1.0, 0.0, 0.0, 0.0]\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.PhaseSpacePoint","page":"Home","title":"QEDprocesses.PhaseSpacePoint","text":"PhaseSpacePoint\n\nRepresentation of a point in the phase space of a process. Contains the process (AbstractProcessDefinition), the model (AbstractModelDefinition), the phase space definition ([AbstractPhasespaceDefinition]), and stateful incoming and outgoing particles (ParticleStateful).\n\nThe legality of the combination of the given process and the incoming and outgoing particles is checked on construction. If the numbers of particles mismatch, the types of particles mismatch (note that order is important), or incoming particles have an Outgoing direction, an error is thrown.\n\njulia> using QEDprocesses; using QEDbase\n\njulia> PhaseSpacePoint(\n            Compton(), \n            PerturbativeQED(), \n            PhasespaceDefinition(SphericalCoordinateSystem(), ElectronRestFrame()), \n            [\n                ParticleStateful(Incoming(), Electron(), SFourMomentum(1, 0, 0, 0)), \n                ParticleStateful(Incoming(), Photon(), SFourMomentum(1, 0, 0, 0))\n            ], \n            [\n                ParticleStateful(Outgoing(), Electron(), SFourMomentum(1, 0, 0, 0)), \n                ParticleStateful(Outgoing(), Photon(), SFourMomentum(1, 0, 0, 0))\n            ]\n        )\nPhaseSpacePoint:\n    process: one-photon Compton scattering\n    model: perturbative QED\n    phasespace definition: spherical coordinates in electron rest frame\n    incoming particles:\n     -> incoming electron (all spins): [1.0, 0.0, 0.0, 0.0]\n     -> incoming photon (all polarizations): [1.0, 0.0, 0.0, 0.0]\n    outgoing particles:\n     -> outgoing electron (all spins): [1.0, 0.0, 0.0, 0.0]\n     -> outgoing photon (all polarizations): [1.0, 0.0, 0.0, 0.0]\n\n\n\n\n\n","category":"type"},{"location":"#QEDprocesses.PhasespaceDefinition","page":"Home","title":"QEDprocesses.PhasespaceDefinition","text":"PhasespaceDefinition(coord_sys::AbstractCoordinateSystem, frame::AbstractFrameOfReference)\n\nConvenient type to dispatch on coordiante systems and frames of reference.\n\n\n\n\n\n","category":"type"},{"location":"#Base.getindex-Tuple{PhaseSpacePoint, QEDbase.Incoming, Int64}","page":"Home","title":"Base.getindex","text":"Base.getindex(psp::PhaseSpacePoint, dir::Incoming, n::Int)\n\nOverload for the array indexing operator []. Returns the nth incoming particle in this phase space point.\n\n\n\n\n\n","category":"method"},{"location":"#Base.getindex-Tuple{PhaseSpacePoint, QEDbase.Outgoing, Int64}","page":"Home","title":"Base.getindex","text":"Base.getindex(psp::PhaseSpacePoint, dir::Outgoing, n::Int)\n\nOverload for the array indexing operator []. Returns the nth outgoing particle in this phase space point.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses._assert_valid_input-Tuple{AbstractComputationSetup, Any}","page":"Home","title":"QEDprocesses._assert_valid_input","text":"_assert_valid_input(stp::AbstractComputationSetup, input::Any)\n\nInterface function, which asserts that the given input is valid, and throws an InvalidInputError if not.\n\nnote: default implementation\nBy default, every input is assumed to be valid. Therefore, this function does nothing.  To customize this behavior, add your own implementation of_assert_valid_input(stp::YourCustomSetup,input)which should throw an exception, which is a subtype of AbstractInvalidInputException. One may also use the concrete implementation InvalidInputError if the input is invalid instead of writing a custom exception type.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses._averaging_norm","page":"Home","title":"QEDprocesses._averaging_norm","text":"_averaging_norm(\n    proc::AbstractProcessDefinition\n    )\n\nInterface function, which returns a normalization for the averaging of the squared matrix elements over spins and polarizations. \n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses._averaging_norm-Tuple{Compton}","page":"Home","title":"QEDprocesses._averaging_norm","text":"note: Convention\nWe average over the initial spins and pols, and sum over final.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses._base_component_type-Union{Tuple{AbstractArray{LV}}, Tuple{LV}} where LV<:QEDbase.AbstractLorentzVector","page":"Home","title":"QEDprocesses._base_component_type","text":"_base_component_type(array_of_lv::AbstractArray{LV}) where {LV<:QEDbase.AbstractLorentzVector}\n\nReturn the type of the components of given Lorentz vectors, which are by themself elements of an  AbstractArray.\n\nExamples\n\njulia> using QEDbase\njulia> using QEDprocesses\njulia> v = Vector{SFourMomentum}(undef,10)\njulia> QEDprocesses._base_component_type(v)\nFloat64\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses._compute","page":"Home","title":"QEDprocesses._compute","text":"_compute(stp::AbstractComputationSetup, input::Any)\n\nInterface function that returns the value of the associated quantity evaluated on input, which can be anything the associated quantity is defined to be feasible for.\n\nnote: unsafe implementation\nThis function must be implemented for any subtype of AbstractComputationSetup. It should not do any input validation or post processing (see _assert_valid_input and _post_processing), as those two are performed while calling  the safe version of this function compute.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses._generate_incoming_momenta","page":"Home","title":"QEDprocesses._generate_incoming_momenta","text":"_generate_incoming_momenta\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition,\n    phase_space_def::AbstractPhasespaceDefinition,\n    in_phase_space::AbstractVector{T},\n) where {T<:Real}\n\nInterface function to generate the four-momenta of the incoming particles from coordinates for a given phase-space definition.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses._generate_momenta-Union{Tuple{T}, Tuple{AbstractProcessDefinition, AbstractModelDefinition, AbstractPhasespaceDefinition, AbstractVector{T}, AbstractVector{T}}} where T<:Real","page":"Home","title":"QEDprocesses._generate_momenta","text":"_generate_momenta(\nproc::AbstractProcessDefinition,\nmodel::AbstractModelDefinition,\nphase_space_def::AbstractPhasespaceDefinition,\nin_phase_space::AbstractVector{T},\nout_phase_space::AbstractVector{T},\n\n) where {T<:Real}\n\nReturn four-momenta for incoming and outgoing particles for given coordinate based phase-space points. \n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses._generate_outgoing_momenta","page":"Home","title":"QEDprocesses._generate_outgoing_momenta","text":"_generate_outgoing_momenta\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition,\n    phase_space_def::AbstractPhasespaceDefinition,\n    out_phase_space::AbstractVector{T},\n) where {T<:Real}\n\nInterface function to generate the four-momenta of the outgoing particles from coordinates for a given phase-space definition.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses._incident_flux","page":"Home","title":"QEDprocesses._incident_flux","text":"_incident_flux(\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition, \n    in_phase_space::AbstractVector{T}\n    ) where {T<:QEDbase.AbstractFourMomentum}\n\nInterface function, which returns the incident flux of the given scattering process for a given incoming phase-space.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses._is_in_phasespace","page":"Home","title":"QEDprocesses._is_in_phasespace","text":"_is_in_phasespace(PhaseSpacePoint{PROC,MODEL}) where {\n    PROC <: AbstractProcessDefinition,\n    MODEL <: AbstractModelDefinition,\n}\n\nInterface function which returns true if the combination of the given incoming and outgoing phase space is physical, i.e. all momenta are on-shell and some sort of energy-momentum conservation holds.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses._matrix_element","page":"Home","title":"QEDprocesses._matrix_element","text":"_matrix_element(PhaseSpacePoint{PROC,MODEL}) where {\n    PROC <: AbstractProcessDefinition,\n    MODEL <: AbstractModelDefinition,\n}\n\nInterface function which returns a tuple of scattering matrix elements for each spin and polarisation combination of proc. \n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses._phase_space_factor","page":"Home","title":"QEDprocesses._phase_space_factor","text":"_phase_space_factor(PhaseSpacePoint{PROC,MODEL,PSDEF}) where {\n    PROC <: AbstractProcessDefinition,\n    MODEL <: AbstractModelDefinition\n    PSDEF <: AbstractPhasespaceDefinition,\n}\n\nInterface function, which returns the pre-differential factor of the invariant phase space intergral measure. \n\nnote: Convention\nIt is assumed, that this function returns the value of mathrmdPi_n= prod_i=1^N fracmathrmd^3p_i(2pi)^3 2 p_i^0 H(P_t p_1 dots p_N)\n\nwhere H(dots) is a characteristic function (or distribution) which constrains the phase space, e.g. delta^(4)(P_t - sum_i p_i).  \n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses._post_processing-Tuple{AbstractComputationSetup, Any, Any}","page":"Home","title":"QEDprocesses._post_processing","text":"function _post_processing(stp::AbstractComputationSetup, input::Any, result::Any)\n\nInterface function, which is called in compute after _compute has been called. This function is dedicated to  finalize the result of a computation. \n\nnote: default implementation\nSince in the case of no post processing the result of _compute is unchanged, this function returns result by default.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses._total_probability","page":"Home","title":"QEDprocesses._total_probability","text":"_total_probability(\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition, \n    ps_def::AbstractPhasespaceDefinition,\n    in_phase_space::AbstractVector{T}\n) where {T<: QEDbase.AbstractFourMomentum}\n\nInterface function for the combination of a scattering process and a physical model. Return the total of a  given process and model for a passed initial phase space definition and point. The elements of in_phase_space, which represent the initial phase space, are the momenta of the respective particles. The implementation of this function for a concrete process and model must not check if the length of  the passed initial phase spaces match the number of incoming particles. \n\nnote: probability interface\nGiven an implementation of this method, the following generic implementation without input check is provided:\n_total_probability(proc_def,model_def,in_phasespace::AbstractMatrix{T})\nwhere T<:QEDbase.AbstractFourMomentum, i.e. _total_probability is also implemented for a vector of initial phase space points. Furthermore, a safe version of _total_probability is also implemented: total_probability.\n\nnote: total cross section\nGiven an implementaion of this method and _incident_flux, the respective functions for the total cross section are available, i.e. _total_cross_section (unsafe and not exported), and total_cross_section, respectively.\n\nnote: Note\nEach instance of this function does not check the validity of the input.  This function is not exported and should be used with caution. To add a method in order to implement the cross section interface,  it is recommented to directly use QEDprocesses._total_cross_section instead of globally using QEDprocesses: _total_cross_section.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.compute-Tuple{AbstractComputationSetup, Any}","page":"Home","title":"QEDprocesses.compute","text":"compute(stp::AbstractComputationSetup, input::Any)\n\nReturn the value of the quantity associated with stp for a given input.  In addition to the actual call of the associated unsafe version _compute, input validation ([_assert_valid_input]) and post processing  (using _post_processing) are wrapped around the calculation (see AbstractComputationSetup for details).\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.differential_cross_section-Tuple{PhaseSpacePoint}","page":"Home","title":"QEDprocesses.differential_cross_section","text":"differential_cross_section(phase_space_point::PhaseSpacePoint)\n\nIf the given phase spaces are physical, return differential cross section evaluated on a phase space point. Zero otherwise.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.differential_probability-Tuple{PhaseSpacePoint}","page":"Home","title":"QEDprocesses.differential_probability","text":"differential_probability(phase_space_point::PhaseSpacePoint)\n\nIf the given phase spaces are physical, return differential probability evaluated on a phase space point. Zero otherwise.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.fundamental_interaction_type","page":"Home","title":"QEDprocesses.fundamental_interaction_type","text":"fundamental_interaction_type(models_def::AbstractModelDefinition)\n\nReturn the fundamental interaction associated with the passed model definition.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.generate_phase_space-Union{Tuple{ElType}, Tuple{AbstractProcessDefinition, AbstractModelDefinition, AbstractPhasespaceDefinition, AbstractVector{ElType}, AbstractVector{ElType}}} where ElType<:QEDbase.AbstractFourMomentum","page":"Home","title":"QEDprocesses.generate_phase_space","text":"generate_phase_space(\n    proc::AbstractProcessDefinition, \n    model::AbstractModelDefinition, \n    ps_def::AbstractPhasespaceDefinition, \n    in_ps::AbstractVector{ElType}, \n    out_ps::AbstractVector{ElType},\n) where {ElType<:QEDbase.AbstractFourMomentum}\n\nReturn the respective phase space point for given momenta of incoming and outgoing particles regarding a given process.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.in_phase_space_dimension","page":"Home","title":"QEDprocesses.in_phase_space_dimension","text":"in_phase_space_dimension(\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition,\n    )\n\nTBW\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.in_phase_space_dimension-Tuple{AbstractProcessDefinition, PerturbativeQED}","page":"Home","title":"QEDprocesses.in_phase_space_dimension","text":"in_phase_space_dimension(proc::AbstractProcessDefinition, ::PerturbativeQED)\n\nReturn the number of degrees of freedom to determine the incoming phase space for processes in PerturbativeQED. \n\nnote: Convention\nThe current implementation only supports the case where two of the incoming particles collide head-on. \n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.incoming_particles","page":"Home","title":"QEDprocesses.incoming_particles","text":"incoming_particles(proc_def::AbstractProcessDefinition)\n\nInterface function for scattering processes. Return a tuple of the incoming particles for the given process definition. This function needs to be given to implement the scattering process interface.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.momentum-Tuple{PhaseSpacePoint, QEDbase.ParticleDirection, Int64}","page":"Home","title":"QEDprocesses.momentum","text":"momentum(psp::PhaseSpacePoint, dir::ParticleDirection, n::Int)\n\nReturns the momentum of the nth particle in the given PhaseSpacePoint which has direction dir. If n is outside the valid range for this phase space point, a BoundsError is thrown.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.number_incoming_particles-Tuple{AbstractProcessDefinition}","page":"Home","title":"QEDprocesses.number_incoming_particles","text":"number_incoming_particles(proc_def::AbstractProcessDefinition)\n\nReturn the number of incoming particles of a given process. \n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.number_outgoing_particles-Tuple{AbstractProcessDefinition}","page":"Home","title":"QEDprocesses.number_outgoing_particles","text":"number_outgoing_particles(proc_def::AbstractProcessDefinition)\n\nReturn the number of outgoing particles of a given process. \n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.out_phase_space_dimension","page":"Home","title":"QEDprocesses.out_phase_space_dimension","text":"out_phase_space_dimension(\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition,\n    )\n\nTBW\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.outgoing_particles","page":"Home","title":"QEDprocesses.outgoing_particles","text":"outgoing_particles(proc_def::AbstractProcessDefinition)\n\nInterface function for scattering processes. Return the tuple of outgoing particles for the given process definition. This function needs to be given to implement the scattering process interface.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.physical_model","page":"Home","title":"QEDprocesses.physical_model","text":"physical_model(stp::AbstractProcessSetup)\n\nInterface function that returns the physical model associated with stp, i.e. an object which is a subtype of AbstractModelDefinition.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.propagator","page":"Home","title":"QEDprocesses.propagator","text":"propagator(particle::AbstractParticleType, mom::QEDbase.AbstractFourMomentum, [mass::Real])\n\nReturn the propagator of a particle for a given four-momentum. If mass is passed, the respective propagator for massive particles is used, if not, it is assumed the particle passed in is massless.\n\nnote: Convention\nThere are two types of implementations for propagators given in QEDProcesses:  For a BosonLike particle with four-momentum k and mass m, the propagator is given as D(k) = frac1k^2 - m^2For a FermionLike particle with four-momentum p and mass m, the propagator is given asS(p) = fracgamma^mu p_mu + massp^2 - m^2\n\nwarning: Warning\nThis function does not throw when the given particle is off-shell. If an off-shell particle is passed, the function propagator returns Inf.\n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.scattering_process","page":"Home","title":"QEDprocesses.scattering_process","text":"scattering_process(stp::AbstractProcessSetup)\n\nInterface function that returns the scattering process associated with stp, i.e. an object which is a subtype of AbstractProcessDefinition. \n\n\n\n\n\n","category":"function"},{"location":"#QEDprocesses.total_cross_section-Union{Tuple{T}, Tuple{AbstractProcessDefinition, AbstractModelDefinition, AbstractPhasespaceDefinition, AbstractVecOrMat{T}}} where T<:QEDbase.AbstractFourMomentum","page":"Home","title":"QEDprocesses.total_cross_section","text":"total_cross_section(\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition,\n    in_phase_space_def::AbstractPhasespaceDefinition,\n    in_phase_space::AbstractVecOrMat{T},\n) where {T<:QEDbase.AbstractFourMomentum}\n\nReturn the total cross section for a given combination of scattering process and compute model, evaluated at the particle momenta.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.total_cross_section-Union{Tuple{T}, Tuple{AbstractProcessDefinition, AbstractModelDefinition, AbstractPhasespaceDefinition, AbstractVecOrMat{T}}} where T<:Real","page":"Home","title":"QEDprocesses.total_cross_section","text":"total_cross_section(\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition,\n    in_phase_space_def::AbstractPhasespaceDefinition,\n    in_phase_space::AbstractVecOrMat{T},\n) where {T<:Real}\n\nReturn the total cross section for a given combination of scattering process and compute model, evaluated at the coordinates.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.total_probability-Union{Tuple{T}, Tuple{AbstractProcessDefinition, AbstractModelDefinition, AbstractPhasespaceDefinition, AbstractVecOrMat{T}}} where T<:QEDbase.AbstractFourMomentum","page":"Home","title":"QEDprocesses.total_probability","text":"total_probability(\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition,\n    phase_space_def::AbstractPhasespaceDefinition,\n    in_phase_space::AbstractMatrix{T},\n) where {T<:QEDbase.AbstractFourMomentum}\n\nReturn the total probability of a given model and process combination, evaluated at the particle momenta.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.total_probability-Union{Tuple{T}, Tuple{AbstractProcessDefinition, AbstractModelDefinition, AbstractPhasespaceDefinition, AbstractVecOrMat{T}}} where T<:Real","page":"Home","title":"QEDprocesses.total_probability","text":"total_probability(\n    proc::AbstractProcessDefinition,\n    model::AbstractModelDefinition,\n    phase_space_def::AbstractPhasespaceDefinition,\n    in_phase_space::AbstractMatrix{T},\n) where {T<:Real}\n\nReturn the total probability of a given model and process combination, evaluated at the coordinates.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.unsafe_differential_cross_section-Tuple{PhaseSpacePoint}","page":"Home","title":"QEDprocesses.unsafe_differential_cross_section","text":"unsafe_differential_cross_section(phase_space_point::PhaseSpacePoint)\n\nReturn the differential cross section evaluated on a phase space point without checking if the given phase space is physical.\n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses.unsafe_differential_probability-Tuple{PhaseSpacePoint}","page":"Home","title":"QEDprocesses.unsafe_differential_probability","text":"unsafe_differential_probability(phase_space_point::PhaseSpacePoint)\n\nReturn differential probability evaluated on a phase space point without checking if the given phase space(s) are physical.\n\n\n\n\n\n","category":"method"}]
}
