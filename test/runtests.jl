using CMInject
using Test

@testset "StarkEffect" begin
    # Validate types
    @test isa(getHamiltonian(1, 5, 1, 1, 1, 1, 1), AbstractMatrix)
    @test isa(calculateEnergies(1, 5, 1, 1, 1, 1, 1), Vector)
    @test isa(calculateStarkCurves(1, 1, 10, 1, 5, 1, 1, 1, 1), Vector)
    @test all(curve -> isa(curve, StarkCurve), calculateStarkCurves(1, 1, 10, 1, 5, 1, 1, 1, 1))
end
