
function Base.show(io::IO, ::QEDbase.AllPol)
    return println(io, "AllPol()")
end

function Base.show(io::IO, ::QEDbase.PolX)
    return println(io, "PolX()")
end

function Base.show(io::IO, ::QEDbase.PolY)
    return println(io, "PolY()")
end

number_of_spin_pol(::AbstractDefinitePolarization) = 1
number_of_spin_pol(::AbstractDefiniteSpin) = 1
number_of_spin_pol(::AbstractIndefinitePolarization) = 2
number_of_spin_pol(::AbstractIndefiniteSpin) = 2
