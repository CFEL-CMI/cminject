using CMInject
using Test

@testset "Stark Simulation" begin
    fieldLines=readlines("test/example_field")
    initial = [parse(Float64, token) for token ∈ split(fieldLines[1], " ") if length(token) > 0]
    stepSizes = [parse(Float64, token) for token ∈ split(fieldLines[2], " ") if length(token) > 0]
    counts = [parse(Int, token) for token ∈ split(fieldLines[3], " ") if length(token) > 0]
    grid = [parse(Float64, fieldLines[4 + y + x * counts[2]]) for y = 0:(counts[2]-1), x = 0:(counts[1]-1)]
    itp = CMInject.interpolate(grid, CMInject.BSpline(CMInject.Quadratic(CMInject.Natural(CMInject.OnGrid()))))
    itpScaled = CMInject.itpscale(itp, initial[1]:stepSizes[1]:(initial[1]+(counts[1]-1)*stepSizes[1]), initial[2]:stepSizes[2]:(initial[2]+(counts[2]-1)*stepSizes[2]))
    dists = Dict(:x => CMInject.Normal(0, 1e-3), :y => CMInject.Dirac(-0.128), :vx => CMInject.Normal(0, 0.1), :vy => CMInject.Dirac(10.0), :m => CMInject.Dirac(2))
    stateDists = Dict(:J => CMInject.DiscreteUniform(2,5), :M => CMInject.DiscreteUniform(0,2))
    # TODO: Add molecule to git
    source = CMInject.StarkSamplingSource{CMInject.StarkParticle2D{Float64}, Float64}(dists, stateDists, "../cmistark/OCS.molecule")

    @test CMInject.generate(source, 1) != nothing

    field = CMInject.ElectricField2d(itpScaled)
    experiment = CMInject.Experiment(;source=source, n_particles=10, fields=(field,), detectors=Tuple(CMInject.SectionDetector{Float64,:y}.(-0.003:0.001:0.003, true)), solver=CMInject.EulerHeun(), time_span=(0.0, 0.03), time_step=1e-5, ensemble_alg=CMInject.EnsembleThreads())
    simRes = CMInject.simulate(experiment)
    data = simRes[2]

    @test simRes != nothing
    @test data != nothing
end
