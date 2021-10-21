using CMInject
using Test

@testset "StarkEffect" begin
    # Validate generic types
    @test all(curve -> isa(curve, CMInject.AbstractInterpolation),
              CMInject.calculateStarkCurves(1, 1, 10, 1, 5, 1, 1, 1, 1))
    @test all(curve -> isa(curve, CMInject.AbstractInterpolation),
              CMInject.calculateStarkCurves(1e-3, 0.0, 1e-2, 0, 5, 0, 9e-26, 3e-26, 1e-29))

    # Validate that Hamiltonians are hermitian
    @test all(x -> abs(x) < 0.000001, 
              CMInject.StarkEffect.getHamiltonian(1, 5, 1, 1, 1, 1, 1)' -
              CMInject.StarkEffect.getHamiltonian(1, 5, 1, 1, 1, 1, 1))
    @test all(x -> abs(x) < 0.000001, 
              CMInject.StarkEffect.getHamiltonian(1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1)' -
              CMInject.StarkEffect.getHamiltonian(1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1))
end
