using CMInject
include("../src/StarkEffect.jl")
using Test

@testset "StarkEffect" begin
    # Validate types
    @test isa(getHamiltonian(1, 5, 1, 1, 1, 1, 1), AbstractMatrix)
    @test isa(getHamiltonian(1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1), AbstractMatrix)
    @test isa(calculateEnergies(1, 5, 1, 1, 1, 1, 1), Vector)
    @test isa(calculateEnergies(1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1), Vector)
    @test isa(calculateStarkCurves(1, 1, 10, 1, 5, 1, 1, 1, 1), Vector)
    @test isa(calculateStarkCurves(1, 1, 10, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1), Vector)
    @test all(curve -> isa(curve, StarkCurve), calculateStarkCurves(1, 1, 10, 1, 5, 1, 1, 1, 1))
    # Validate that Hamiltonians are hermitian
    @test all(x -> abs(x) < 0.000001, getHamiltonian(1, 5, 1, 1, 1, 1, 1)' -
              getHamiltonian(1, 5, 1, 1, 1, 1, 1))
    @test all(x -> abs(x) < 0.000001, getHamiltonian(1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1)' -
              getHamiltonian(1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1))
end
