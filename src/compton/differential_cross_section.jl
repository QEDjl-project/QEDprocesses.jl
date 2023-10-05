
const ComptonDCS = DifferentialCrossSection{
    Compton{Pol},
    PerturbativeQED,
    PhaseSpace,
    NumericType,
} where {
    Pol<:AbstractPolarization,
    NumericType<:SFourMomentum,
    PhaseSpace<:AbstractArray{SFourMomentum},
}
