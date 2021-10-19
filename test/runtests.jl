using CMInject
include("../src/StarkEffect.jl")
using Test

@testset "StarkEffect" begin
    # Validate that Hamiltonians are hermitian
    @test all(x -> abs(x) < 0.000001, getHamiltonian(1, 5, 1, 1, 1, 1, 1)' -
              getHamiltonian(1, 5, 1, 1, 1, 1, 1))
    @test all(x -> abs(x) < 0.000001, getHamiltonian(1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1)' -
              getHamiltonian(1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1))
end
