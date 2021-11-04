using CMInject
using Test

@testset "Stark Simulation" begin
    itp = CMInject.interpolate([2^2+2^2 1^2+2^2 0^2+2^2 1^2+2^2 2^2+2^2; 2^2+1^2 1^2+1^2 0^2+1^2 1^2+1^2 2^2+1^2; 2^2+0^2 1^2+0^2 0^2+0^2 1^2+0^2 2^2+0^2; 2^2+1^2 1^2+1^2 0^2+1^2 1^2+1^2 2^2+1^2; 2^2+2^2 1^2+2^2 0^2+2^2 1^2+2^2 2^2+2^2], CMInject.BSpline(CMInject.Quadratic(CMInject.Natural(CMInject.OnGrid()))))
    itpScaled = CMInject.itpscale(itp, -2:1:2, -2:1:2)
    dists = Dict(:x => CMInject.Normal(0, 1e-3), :y => CMInject.Dirac(-0.128), :vx => CMInject.Normal(0, 0.1), :vy => CMInject.Dirac(10.0), :m => CMInject.Dirac(2))
    stateDists = Dict(:J => CMInject.DiscreteUniform(2,5), :M => CMInject.DiscreteUniform(0,2))
    # TODO: These are not appropriate properties for the current particle
    properties = Dict(:B => 1, :D => 2, :μ => 3)
    source = CMInject.StarkSamplingSource{CMInject.StarkParticle2D{Float64}, Float64}(dists, stateDists, properties)

    @test CMInject.generate(source, 1) != nothing

    field = CMInject.ElectricField2d(itpScaled)
    experiment = CMInject.Experiment(;source=source, n_particles=10, fields=(field,), detectors=Tuple(CMInject.SectionDetector{Float64,:y}.(-0.003:0.001:0.003, true)), solver=CMInject.EulerHeun(), time_span=(0.0, 0.03), time_step=1e-5, ensemble_alg=CMInject.EnsembleThreads())
    simRes = CMInject.simulate(experiment)
    data = simRes[2]

    @test simRes != nothing
    @test data != nothing
end
