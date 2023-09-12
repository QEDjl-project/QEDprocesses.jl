using QEDprocesses

@testset "fermion likes" begin
    @testset "fermion" begin
        struct TestFermion <: Fermion end
        @test is_fermion(TestFermion())
        @test is_particle(TestFermion())
        @test !is_anti_particle(TestFermion())
    end

    @testset "antifermion" begin
        struct TestAntiFermion <: AntiFermion end
        @test is_fermion(TestAntiFermion())
        @test !is_particle(TestAntiFermion())
        @test is_anti_particle(TestAntiFermion())
    end

    @testset "majorana fermion" begin
        struct TestMajoranaFermion <: MajoranaFermion end
        @test is_fermion(TestMajoranaFermion())
        @test is_particle(TestMajoranaFermion())
        @test is_anti_particle(TestMajoranaFermion())
    end

    @testset "electron" begin
        @test is_fermion(Electron())
        @test is_particle(Electron())
        @test !is_anti_particle(Electron())
        @test mass(Electron()) == 1.0
        @test charge(Electron()) == -1.0
    end

    @testset "positron" begin
        @test is_fermion(Positron())
        @test !is_particle(Positron())
        @test is_anti_particle(Positron())
        @test mass(Positron()) == 1.0
        @test charge(Positron()) == 1.0
    end
end

@testset "boson likes" begin
    @testset "boson" begin
        struct TestBoson <: Boson end
        @test !is_fermion(TestBoson())
        @test is_boson(TestBoson())
        @test is_particle(TestBoson())
        @test !is_anti_particle(TestBoson())
    end

    @testset "antiboson" begin
        struct TestAntiBoson <: AntiBoson end
        @test !is_fermion(TestAntiBoson())
        @test is_boson(TestAntiBoson())
        @test !is_particle(TestAntiBoson())
        @test is_anti_particle(TestAntiBoson())
    end

    @testset "majorana boson" begin
        struct TestMajoranaBoson <: MajoranaBoson end
        @test !is_fermion(TestMajoranaBoson())
        @test is_boson(TestMajoranaBoson())
        @test is_particle(TestMajoranaBoson())
        @test is_anti_particle(TestMajoranaBoson())
    end

    @testset "photon" begin
        @test !is_fermion(Photon())
        @test is_boson(Photon())
        @test is_particle(Photon())
        @test is_anti_particle(Photon())
        @test charge(Photon()) == 0.0
        @test mass(Photon()) == 0.0
    end
end
