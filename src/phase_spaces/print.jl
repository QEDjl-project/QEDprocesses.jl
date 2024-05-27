function Base.show(io::IO, ::SphericalCoordinateSystem)
    print(io, "spherical coordinates")
    return nothing
end

function Base.show(io::IO, ::CenterOfMomentumFrame)
    print(io, "center-of-momentum frame")
    return nothing
end

function Base.show(io::IO, ::ElectronRestFrame)
    print(io, "electron rest frame")
    return nothing
end

function Base.show(io::IO, m::MIME"text/plain", ps_def::PhasespaceDefinition)
    println(io, "PhasespaceDefinition")
    println(io, "    coordinate system: $(ps_def.coord_sys)")
    println(io, "    frame: $(ps_def.frame)")
    return nothing
end

function Base.show(io::IO, ps_def::PhasespaceDefinition)
    print(io, "$(ps_def.coord_sys) in $(ps_def.frame)")
    return nothing
end

function Base.show(io::IO, particle::ParticleStateful)
    print(io, "$(particle.dir) $(particle.species): $(particle.mom)")
    return nothing
end

function Base.show(io::IO, m::MIME"text/plain", particle::ParticleStateful)
    println(io, "ParticleStateful: $(particle.dir) $(particle.species)")
    println(io, "    momentum: $(particle.mom)")
    return nothing
end

function Base.show(io::IO, psp::PhaseSpacePoint)
    print(io, "PhaseSpacePoint of $(psp.proc)")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", psp::PhaseSpacePoint)
    println(io, "PhaseSpacePoint:")
    println(io, "    process: $(psp.proc)")
    println(io, "    model: $(psp.model)")
    println(io, "    phasespace definition: $(psp.ps_def)")
    println(io, "    incoming particles:")
    for p in psp.in_particles
        println(io, "     -> $(p)")
    end
    println(io, "    outgoing particles:")
    for p in psp.out_particles
        println(io, "     -> $(p)")
    end
    return nothing
end
