
# dummy particles
struct TestFermion <: FermionLike end
QEDbase.mass(::TestFermion) = 1.0
QEDbase.charge(::TestFermion) = 2.0
struct TestBoson <: BosonLike end
QEDbase.mass(::TestBoson) = 0.0
QEDbase.charge(::TestBoson) = -2.0

const PARTICLE_SET = [TestFermion(), TestBoson()]

"""

    TestProcess(rng,incoming_particles,outgoing_particles)

"""
struct TestProcess{IP<:AbstractVector,OP<:AbstractVector} <: AbstractProcessDefinition
    incoming_particles::IP
    outgoing_particles::OP
end

function TestProcess(rng::AbstractRNG, N_in::Int, N_out::Int)
    in_particles = rand(rng, PARTICLE_SET, N_in)
    out_particles = rand(rng, PARTICLE_SET, N_out)
    return TestProcess(in_particles, out_particles)
end

QEDprocesses.incoming_particles(proc::TestProcess) = proc.incoming_particles
QEDprocesses.outgoing_particles(proc::TestProcess) = proc.outgoing_particles

struct TestPhasespaceDef <: AbstractPhasespaceDefinition end

struct TestModel <: AbstractModelDefinition end
QEDprocesses.fundamental_interaction_type(::TestModel) = :test_interaction
