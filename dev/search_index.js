var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = QEDprocesses","category":"page"},{"location":"#QEDprocesses","page":"Home","title":"QEDprocesses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for QEDprocesses.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [QEDprocesses]","category":"page"},{"location":"#QEDprocesses.Compton","page":"Home","title":"QEDprocesses.Compton","text":"Compton(\n    in_spin [= AllSpin()]\n    in_pol [= AllPol()]\n    out_spin [= AllSpin()]\n    out_pol [= AllPol()]\n)\n\n\n\n\n\n","category":"type"},{"location":"#QEDbase._averaging_norm-Tuple{Compton}","page":"Home","title":"QEDbase._averaging_norm","text":"_averaging_norm(proc::Compton)\n\nnote: Convention\nWe average over the initial spins and pols, and sum over final.\n\n\n\n\n\n","category":"method"},{"location":"#QEDbase.in_phase_space_dimension-Tuple{QEDbase.AbstractProcessDefinition, PerturbativeQED}","page":"Home","title":"QEDbase.in_phase_space_dimension","text":"in_phase_space_dimension(proc::AbstractProcessDefinition, ::PerturbativeQED)\n\nReturn the number of degrees of freedom to determine the incoming phase space for processes in PerturbativeQED. \n\nnote: Convention\nThe current implementation only supports the case where two of the incoming particles collide head-on. \n\n\n\n\n\n","category":"method"},{"location":"#QEDprocesses._base_component_type-Union{Tuple{AbstractArray{LV}}, Tuple{LV}} where LV<:QEDbase.AbstractLorentzVector","page":"Home","title":"QEDprocesses._base_component_type","text":"_base_component_type(array_of_lv::AbstractArray{LV}) where {LV<:AbstractLorentzVector}\n\nReturn the type of the components of given Lorentz vectors, which are by themself elements of an  AbstractArray.\n\nExamples\n\njulia> using QEDbase\njulia> using QEDprocesses\njulia> v = Vector{SFourMomentum}(undef,10)\njulia> QEDprocesses._base_component_type(v)\nFloat64\n\n\n\n\n\n","category":"method"}]
}
