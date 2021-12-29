@testset "Fly Stark Simulation" begin
    # TODO: Read and shift field
    grid = [1e7 2e7; 1e7 2e7]
    itp = CMInject.interpolate(grid, CMInject.BSpline(CMInject.Linear()))
    ext = CMInject.extrapolate(itp, 0)
    itpScaled = CMInject.itpscale(ext, 0:1:1, -0.5:1:0.5)

    dists = Dict(:x => CMInject.Normal(0, 0.0001*0.0001),
                 :y => CMInject.Dirac(0),
                 :vx => CMInject.Normal(0, 1.5*1.5),
                 :vy => CMInject.Normal(670, 7*7),
                 :m => CMInject.Dirac((4*12+7*1.0078+14.003+15.995)/9.223e18))
    stateDists = Dict(:J => CMInject.DiscreteUniform(0,0), :M => CMInject.DiscreteUniform(0,0))
    sourcePyrroleWater = CMInject.StarkSamplingSource{CMInject.StarkParticle2D{Float64}, Float64}(
                                     dists, stateDists, "test/pyrrole-water.molecule")

    @test CMInject.generate(sourcePyrroleWater, 1) != nothing

    field = CMInject.ElectricField2D(itpScaled)
    particles = 10
    @inline function detectHits(args...)
        x = args[2].x
        y = args[2].y
        # TODO: There are contradictoray information for these values
        skimmer1Y = 0.065
        skimmer2Y = 0.302
        knifeY    = 0.512
        knifeX    = 0
        skimmer3Y = 0.537
        skimmer1R = 0.003
        skimmer2R = 0.0015
        skimmer3R = 0.0015
        # TODO: Don't ignore deflector hits
        if (y ≥ skimmer1Y-0.01 && y ≤ skimmer1Y &&
            abs(x) ≤ skimmer1R)
            print("HIT\n")
            return true;
        end
        if (y ≥ skimmer2Y-0.01 && y ≤ skimmer2Y &&
            abs(x) ≤ skimmer2R)
            print("HIT\n")
            return true;
        end
        if (y ≥ skimmer3Y-0.01 && y ≤ skimmer3Y &&
            abs(x) ≤ skimmer3R)
            print("HIT\n")
            return true;
        end
        if (y ≥ knifeY-0.01 && y ≤ knifeY &&
            x ≤ knifeX)
            print("HIT\n")
            return true;
        end
        false
    end
    experimentPyrroleWater = CMInject.Experiment(;source=sourcePyrroleWater,
                                            n_particles=particles,
                                            fields=(field,),
                                            detectors=Tuple(CMInject.SectionDetector{Float64,:y}.(
                                                             -0.001:0.001:0.001,
                                                             true)),
                                            solver=CMInject.EulerHeun(),
                                            # TODO: Think about reasonable time spans
                                            time_span=(0.0, 1.0),
                                            time_step=1e-1,
                                            ensemble_alg=CMInject.EnsembleThreads(),
                                            solver_opts=(adaptive=false, dense=false, unstable_check=detectHits))

    simResPyrroleWater = CMInject.simulate(experimentPyrroleWater)
    dataPyrroleWater = simResPyrroleWater[1]

    @test simResPyrroleWater != nothing
    @test dataPyrroleWater != nothing

    for i in 1:particles
        # TODO: Validate results
        print(last(dataPyrroleWater[i].u).x, ", ", last(dataPyrroleWater[i].u).y, "\n")
    end
end

