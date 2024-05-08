using QEDbase
using QEDprocesses
using BenchmarkTools
using StaticArrays
using Random

include("test_implementation.jl")

const SUITE = BenchmarkGroup()


SUITE["stateful"]=BenchmarkGroup()

MOM_IN_1 = rand(SFourMomentum)
MOM_IN_1K = rand(SFourMomentum,1000)
MOM_IN_1M= rand(SFourMomentum,1000000)
SUITE["stateful"]["construction"] = BenchmarkGroup()
SUITE["stateful"]["construction"]["single"] = @benchmarkable ParticleStateful(Incoming(),TestFermion(),mom,AllSpin()) setup=(mom=$MOM_IN_1)
SUITE["stateful"]["construction"]["broadcast 1e3"] = @benchmarkable ParticleStateful.(Incoming(),TestFermion(),mom,AllSpin()) setup=(mom=$MOM_IN_1K)
SUITE["stateful"]["construction"]["broadcast 1e6"] = @benchmarkable ParticleStateful.(Incoming(),TestFermion(),mom,AllSpin()) setup=(mom=$MOM_IN_1M)

PARTICLE_STATEFUL = ParticleStateful(Incoming(), TestFermion(), MOM_IN_1,AllSpin())
SUITE["stateful"]["accessor"] = BenchmarkGroup()
SUITE["stateful"]["accessor"]["is_incoming"] = @benchmarkable QEDprocesses.is_incoming(particle) setup=(particle=$PARTICLE_STATEFUL)
SUITE["stateful"]["accessor"]["is_outgoing"] = @benchmarkable QEDprocesses.is_outgoing(particle) setup=(particle=$PARTICLE_STATEFUL)
SUITE["stateful"]["accessor"]["is_fermion"] = @benchmarkable is_fermion(particle) setup=(particle=$PARTICLE_STATEFUL)
SUITE["stateful"]["accessor"]["is_boson"] = @benchmarkable is_boson(particle) setup=(particle=$PARTICLE_STATEFUL)
SUITE["stateful"]["accessor"]["is_particle"] = @benchmarkable is_particle(particle) setup=(particle=$PARTICLE_STATEFUL)
SUITE["stateful"]["accessor"]["is_anti_particle"] = @benchmarkable is_anti_particle(particle) setup=(particle=$PARTICLE_STATEFUL)
SUITE["stateful"]["accessor"]["mass"]  = @benchmarkable mass(particle) setup=(particle=$PARTICLE_STATEFUL)
SUITE["stateful"]["accessor"]["charge"] = @benchmarkable charge(particle) setup=(particle=$PARTICLE_STATEFUL)


MODEL = TestModel()
PSDEF = TestPhasespaceDef()

PROCESS_STRINGS = Dict(
    "BF->BF" => TestProcess(SVector(TestBoson(), TestFermion()), SVector(TestBoson(), TestFermion())),
    "BB->BB" => TestProcess(SVector(TestBoson(), TestBoson()), SVector(TestBoson(), TestBoson())),
    "FF->FF" => TestProcess(SVector(TestFermion(), TestFermion()), SVector(TestFermion(), TestFermion())),
    "FF->BF" => TestProcess(SVector(TestFermion(), TestFermion()), SVector(TestBoson(), TestFermion())),
    "FB->FF" => TestProcess(SVector(TestBoson(), TestFermion()), SVector(TestFermion(), TestFermion())),
)

SUITE["phase space point"] = BenchmarkGroup()
SUITE["phase space point"]["generate"] = BenchmarkGroup()
for (proc_string,proc) in PROCESS_STRINGS
    in_moms = rand(SFourMomentum,number_incoming_particles(proc))
    out_moms = rand(SFourMomentum,number_outgoing_particles(proc))

    in_particles= incoming_particles(proc)
    out_particles= outgoing_particles(proc)

    in_particle_states = SVector(collect(ParticleStateful(Incoming(),in_particles[i],in_moms[i]) for i in 1:number_incoming_particles(proc))...)
    out_particle_states = SVector(collect(ParticleStateful(Outgoing(),out_particles[i],out_moms[i]) for i in 1:number_outgoing_particles(proc))...)

    SUITE["phase space point"]["construction"][proc_string]= @benchmarkable PhaseSpacePoint($proc,$MODEL,$PSDEF,$in_particle_states,$out_particle_states)
    SUITE["phase space point"]["generate"][proc_string]= @benchmarkable generate_phase_space($proc,$MODEL,$PSDEF,$in_moms,$out_moms)

    psp = generate_phase_space(proc,MODEL,PSDEF,in_moms,out_moms)
    SUITE["phase space point"]["momentum"][proc_string]= @benchmarkable momentum($psp,Incoming(),1)
    SUITE["phase space point"]["indexing"][proc_string]= @benchmarkable $psp[Incoming(),1]

end
 


