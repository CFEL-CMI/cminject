using Plots
histogramData = [[], []]
# TODO: Reduce duplicates resulting from Pyrrole AND Pyrrole Water
@testset "Fly Stark Simulation" begin
    # TODO: There are contradictoray information for these values
    skimmer1Y = 0.065
    skimmer2Y = 0.302
    # This value comes from CMIfly and is not available in the paper
    deflectorY= 0.437
    knifeY    = 0.512
    knifeX    = 0
    skimmer3Y = 0.537
    skimmer1R = 0.003
    skimmer2R = 0.0015
    skimmer3R = 0.0015
    destination = 0.778

    fieldLines=readlines("test/example_field")
    initial = [parse(Float64, token) for token ∈ split(fieldLines[1], " ") if length(token) > 0]
    initial[2] += 0.367
    stepSizes = [parse(Float64, token) for token ∈ split(fieldLines[2], " ") if length(token) > 0]
    counts = [parse(Int, token) for token ∈ split(fieldLines[3], " ") if length(token) > 0]
    grid = [parse(Float64, fieldLines[4 + y + x * counts[2]]) for y = 0:(counts[2]-1), x = 0:(counts[1]-1)]
    itp = CMInject.interpolate(grid, CMInject.BSpline(CMInject.Linear()))
    ext = CMInject.extrapolate(itp, 0)
    itpScaled = CMInject.itpscale(ext, initial[1]:stepSizes[1]:(initial[1]+(counts[1]-1)*stepSizes[1]),
                                  (initial[2]+deflectorY):stepSizes[2]:(initial[2]+deflectorY+(counts[2]-1)*stepSizes[2]))
    # TODO: Read & use gradient

    distsPyrrole = Dict(:x => CMInject.Normal(0, 0.0001*0.0001),
                        :y => CMInject.Dirac(0),
                        :vx => CMInject.Normal(0, 1.5*1.5),
                        :vy => CMInject.Normal(670, 7*7),
                        # C4H5N
                        :m => CMInject.Dirac((4*12+5*1.0078+14.003)/9.223e18))
    distsPyrroleWater = Dict(:x => CMInject.Normal(0, 0.0001*0.0001),
                             :y => CMInject.Dirac(0),
                             :vx => CMInject.Normal(0, 1.5*1.5),
                             :vy => CMInject.Normal(670, 7*7),
                             # C4H7NO
                             :m => CMInject.Dirac((4*12+7*1.0078+14.003+15.995)/9.223e18))
    # TODO: Allow different states
    stateDists = Dict(:J => CMInject.DiscreteUniform(0,0), :M => CMInject.DiscreteUniform(0,0))
    sourcePyrroleWater = CMInject.StarkSamplingSource{CMInject.StarkParticle2D{Float64}, Float64}(
                                     distsPyrroleWater, stateDists, "test/pyrrole-water.molecule")
    sourcePyrrole = CMInject.StarkSamplingSource{CMInject.StarkParticle2D{Float64}, Float64}(
                                     distsPyrrole, stateDists, "test/pyrrole.molecule")

    @test CMInject.generate(sourcePyrroleWater, 1) != nothing

    field = CMInject.ElectricField2D(itpScaled)
    particles = 10000
    @inline function detectHits(args...)
        x = args[2].x
        y = args[2].y
        # TODO: Don't ignore deflector hits
        if (y ≥ skimmer1Y-0.01 && y ≤ skimmer1Y &&
            abs(x) ≥ skimmer1R)
            print("HIT skimmer 1\n")
            return true;
        end
        if (y ≥ skimmer2Y-0.01 && y ≤ skimmer2Y &&
            abs(x) ≥ skimmer2R)
            print("HIT skimmer 2\n")
            return true;
        end
        if (y ≥ skimmer3Y-0.01 && y ≤ skimmer3Y &&
            abs(x) ≥ skimmer3R)
            print("HIT skimmer 3\n")
            return true;
        end
        if (y ≥ knifeY-0.01 && y ≤ knifeY &&
            x ≥ knifeX)
            print("HIT knife\n")
            return true;
        end
        false
    end
    experimentPyrroleWater = CMInject.Experiment(;source=sourcePyrroleWater,
                                            n_particles=particles,
                                            fields=(field,),
                                            # TODO: Figure out how to use detectors
                                            detectors=Tuple(CMInject.SectionDetector{Float64,:y}.(
                                                             -0.001:0.001:0.001,
                                                             true)),
                                            solver=CMInject.EulerHeun(),
                                            time_span=(0.0, 2e-3),
                                            time_step=1e-6,
                                            ensemble_alg=CMInject.EnsembleThreads(),
                                            solver_opts=(adaptive=false, dense=false, unstable_check=detectHits))
    experimentPyrrole = CMInject.Experiment(;source=sourcePyrrole,
                                            n_particles=particles,
                                            fields=(field,),
                                            detectors=Tuple(CMInject.SectionDetector{Float64,:y}.(
                                                             -0.001:0.001:0.001,
                                                             true)),
                                            solver=CMInject.EulerHeun(),
                                            time_span=(0.0, 2e-3),
                                            time_step=1e-6,
                                            ensemble_alg=CMInject.EnsembleThreads(),
                                            solver_opts=(adaptive=false, dense=false, unstable_check=detectHits))

    simResPyrroleWater = CMInject.simulate(experimentPyrroleWater)
    simResPyrrole = CMInject.simulate(experimentPyrrole)
    dataPyrroleWater = simResPyrroleWater[1]
    dataPyrrole = simResPyrrole[1]

    @test simResPyrroleWater != nothing
    @test simResPyrrole != nothing
    @test dataPyrroleWater != nothing
    @test dataPyrrole != nothing

    print("PYRROLE WATER deflections:\n██████████████████████████\n")
    averageDeflectionPyrroleWater = 0
    finalCountPyrroleWater = 0
    for i in 1:particles
        # TODO: Validate results
        for u in dataPyrroleWater[i].u
            if (u.y ≥ destination-0.005 && u.y ≤ destination+0.005)
                finalCountPyrroleWater += 1
                averageDeflectionPyrroleWater += u.x
                print(u.x, ", ", u.y, "\n")
                push!(histogramData[1], u.x)
                break
            end
        end
    end
    averageDeflectionPyrroleWater /= finalCountPyrroleWater;
    print("-> average deflection: ", averageDeflectionPyrroleWater, ", total: ", finalCountPyrroleWater, "\n");
    print("PYRROLE deflections:\n██████████████████████████\n")
    averageDeflectionPyrrole = 0
    finalCountPyrrole = 0
    for i in 1:particles
        # TODO: Validate results
        for u in dataPyrrole[i].u
            if (u.y ≥ destination-0.005 && u.y ≤ destination+0.005)
                finalCountPyrrole += 1
                averageDeflectionPyrrole += u.x
                print(u.x, ", ", u.y, "\n")
                if (size(histogramData[2]) < size(histogramData[1]))
                    push!(histogramData[2], u.x)
                end
                break
            end
        end
    end
    averageDeflectionPyrrole /= finalCountPyrrole;
    print("-> average deflection: ", averageDeflectionPyrrole, ", total: ", finalCountPyrrole, "\n");
end
plot(histogramData[2], seriestype=:histogram, nbins=10, fillalpha=0.5, labels=["Pyrrole Water", "Pyrrole"][2])
plot!(histogramData[1], seriestype=:histogram, nbins=10, fillalpha=0.5, labels=["Pyrrole Water", "Pyrrole"][1])

