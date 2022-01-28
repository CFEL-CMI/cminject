using Plots
histogramData = [[], []]
trajectories = [[[]], [[]]]
# TODO: Reduce duplicates resulting from Pyrrole AND Pyrrole Water

original0 = 0.34572
skimmer1Z     = -0.28072 + original0
skimmer2Z     = -0.044 + original0
deflectorZ    = 0 + original0
# 154mm as Sebastian stated
deflectorEndZ = 0.154 + original0
knifeZ        = 0.16528 + original0
skimmer3Z     = 0.19041 + original0
destination   = 0.36675 + original0
skimmer1R     = 0.0015
skimmer2R     = 0.00075
skimmer3R     = 0.00075

sigmaPos      = 0.00052
limitPos      = 2e-3
sigmaVtrans   = 10
limitVtrans   = 10
sigmaVlong    = 20
limitVlong    = 50
meanVz        = 1860

# FIELD START

fieldLines=readlines("test/example_field")
initial = [parse(Float64, token) for token ∈ split(fieldLines[1], " ") if length(token) > 0]
stepSizes = [parse(Float64, token) for token ∈ split(fieldLines[2], " ") if length(token) > 0]
counts = [parse(Int, token) for token ∈ split(fieldLines[3], " ") if length(token) > 0]

initial[3] = deflectorZ
stepSizes[3] = deflectorEndZ - deflectorZ
counts[3] = 2

myGrid = [parse(Float64, fieldLines[4 + y + x * counts[2]]) *
          # The field was stored with 60kV instead of 14kV
          14.0 / 60.0
          # y (later to be x) is stored in reverse order
          for y = (counts[2]-1):-1:0, x = 0:(counts[1]-1), z = 0:(counts[3]-1)]
print("Grid is: ", size(myGrid), "\n")
itp = CMInject.interpolate(myGrid, CMInject.BSpline(CMInject.Linear()))
ext = CMInject.extrapolate(itp, 0)
itpScaled = CMInject.itpscale(ext, 
                              initial[2]:stepSizes[2]:(initial[2]+(counts[2]-1)*stepSizes[2]),
                              initial[1]:stepSizes[1]:(initial[1]+(counts[1]-1)*stepSizes[1]),
                              initial[3]:stepSizes[3]:deflectorEndZ)

# FIELD END

# GRADIENT START

gradFieldLines=readlines("test/example_gradient")
gradInitial = [parse(Float64, token) for token ∈ split(gradFieldLines[1], " ") if length(token) > 0]
gradStepSizes = [parse(Float64, token) for token ∈ split(gradFieldLines[2], " ") if length(token) > 0]
gradCounts = [parse(Int, token) for token ∈ split(gradFieldLines[3], " ") if length(token) > 0]

gradInitial[3] = deflectorZ
gradStepSizes[3] = deflectorEndZ - deflectorZ
gradCounts[3] = 2

gradGrid = [[map(x->parse(Float64, x) * 
                 # The field was stored with 60kV instead of 14kV
                 14.0 / 60.0,
                 split(gradFieldLines[4 + y + x * gradCounts[2]], " ", keepempty=false))[i]
             # y (later to be x) is stored in reverse order
             for y = (gradCounts[2]-1):-1:0, x = 0:(gradCounts[1]-1), z = 0:(gradCounts[3]-1)] for i = 1:3]
# TODO: Validate that I can just mirror (instead of rotate) it
# x and y are wrong :/
tmp = gradGrid[1]
gradGrid[1] = gradGrid[2]
gradGrid[2] = tmp
gradItps = @. CMInject.interpolate(gradGrid, CMInject.BSpline(CMInject.Linear()))
gradExts = @. CMInject.extrapolate(gradItps, 0)
gradItpScaleds = tuple([CMInject.itpscale(gradExts[i], 
                              initial[2]:stepSizes[2]:(initial[2]+(counts[2]-1)*stepSizes[2]),
                              initial[1]:stepSizes[1]:(initial[1]+(counts[1]-1)*stepSizes[1]),
                              initial[3]:stepSizes[3]:deflectorEndZ) for i = 1:3]...)


# GRADIENT END

print("Initialization done\n")

@testset "Fly Stark Simulation" begin

    u = 1.66053906660e-27
    truncatedNormal(mean, sigma, limit) = CMInject.Truncated(CMInject.Normal(mean, sigma), mean-limit, mean+limit)
    distsPyrrole = Dict(:x => truncatedNormal(0, sigmaPos, limitPos),
                        :y => truncatedNormal(0, sigmaPos, limitPos),
                        :z => CMInject.Dirac(0),
                        :vx => truncatedNormal(0, sigmaVtrans, limitVtrans),
                        :vy => truncatedNormal(0, sigmaVtrans, limitVtrans),
                        :vz => truncatedNormal(meanVz, sigmaVlong, limitVlong),
                        # C4H5N
                        :m => CMInject.Dirac((4*12+5*1.0078+14.003)*u))
    distsPyrroleWater = Dict(:x => truncatedNormal(0, sigmaPos, limitPos),
                             :y => truncatedNormal(0, sigmaPos, limitPos),
                             :z => CMInject.Dirac(0),
                             :vx => truncatedNormal(0, sigmaVtrans, limitVtrans),
                             :vy => truncatedNormal(0, sigmaVtrans, limitVtrans),
                             :vz => truncatedNormal(meanVz, sigmaVlong, limitVlong),
                             # C4H7NO
                             :m => CMInject.Dirac((4*12+7*1.0078+14.003+15.995)*u))
    # TODO: Allow different states
    stateDists = Dict(:J => CMInject.DiscreteUniform(0,0), :M => CMInject.DiscreteUniform(0,0))
    sourcePyrroleWater = CMInject.StarkSamplingSource{CMInject.StarkParticle{Float64}, Float64}(
                                     distsPyrroleWater, stateDists, "test/pyrrole-water.molecule")
    sourcePyrrole = CMInject.StarkSamplingSource{CMInject.StarkParticle{Float64}, Float64}(
                                     distsPyrrole, stateDists, "test/pyrrole.molecule")

    @test CMInject.generate(sourcePyrroleWater, 1) != nothing

    field = CMInject.ElectricField(itpScaled, gradItpScaleds)
    particles = 10000
    @inline function detectHits(args...)
        x = args[1].x
        y = args[1].y
        z = args[1].z
        if (z ≥ skimmer1Z-0.01 && z ≤ skimmer1Z &&
            sqrt(x*x + y*y) ≥ skimmer1R)
            #print("HIT skimmer 1\n")
            return true
        end
        if (z ≥ skimmer2Z-0.01 && z ≤ skimmer2Z &&
            sqrt(x*x + y*y) ≥ skimmer2R)
            #print("HIT skimmer 2\n")
            return true
        end
        if (z ≥ skimmer3Z-0.01 && z ≤ skimmer3Z &&
            sqrt(x*x + y*y) ≥ skimmer3R)
            #print("HIT skimmer 3\n")
            return true
        end
        if (z ≥ knifeZ-0.01 && z ≤ knifeZ &&
            y ≤ knifeY)
            #print("HIT knife\n")
            return true
        end
        if (z ≥ deflectorZ && z ≤ deflectorEndZ &&
            (x ≤ initial[1] || x ≥ initial[1]+(counts[1]-1)*stepSizes[1] ||
             y ≤ initial[2] || y ≥ initial[2]+(counts[2]-1)*stepSizes[2]))
            #print("HIT deflector\n")
            return true
        end
        false
    end
    solverOpts=(adaptive=false,
                dense=false,
                callback=CMInject.DiscreteCallback(detectHits, integrator -> CMInject.terminate!(integrator)))
    experimentPyrroleWater = CMInject.Experiment(;source=sourcePyrroleWater,
                                            n_particles=particles,
                                            fields=(field,),
                                            # TODO: Figure out how to use detectors
                                            detectors=Tuple(CMInject.SectionDetector{Float64,:y}.(
                                                             -0.001:0.001:0.001,
                                                             true)),
                                            # TODO: Try with Runge-Kutta 4/5
                                            solver=CMInject.EulerHeun(),
                                            # TODO: If necessary, adjust the time span
                                            time_span=(0.0, 2e-3),
                                            time_step=1e-6,
                                            ensemble_alg=CMInject.EnsembleThreads(),
                                            solver_opts=solverOpts)
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
                                            solver_opts=solverOpts)

    knifeY = 0#3.188356743548741e-6
    simResPyrroleWater = CMInject.simulate(experimentPyrroleWater)
    knifeY = 0#1.5357245757532834e-6
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
    averageKnifePyrroleWater = 0
    countKnifePyrroleWater = 0
    for i in 1:particles
        push!(trajectories[1], [])
        # TODO: Validate results
        for u in dataPyrroleWater[i].u
            if (u.z ≥ knifeZ-0.001 && u.z ≤ knifeZ+0.001)
                averageKnifePyrroleWater += u.y
                countKnifePyrroleWater += 1
            end
            if (u.z ≥ destination-0.001 && u.z ≤ destination+0.005)
                finalCountPyrroleWater += 1
                averageDeflectionPyrroleWater += u.y
                push!(histogramData[1], u.y)
                break
            end
            push!(trajectories[1][i], (u.z, u.y))
        end
    end
    averageDeflectionPyrroleWater /= finalCountPyrroleWater;
    averageKnifePyrroleWater /= countKnifePyrroleWater;
    print("-> average knife pyrrole water: ", averageKnifePyrroleWater, "\n")
    print("-> average deflection: ", averageDeflectionPyrroleWater, ", total: ", finalCountPyrroleWater, "\n");
    print("PYRROLE deflections:\n██████████████████████████\n")
    averageDeflectionPyrrole = 0
    finalCountPyrrole = 0
    averageKnifePyrrole = 0
    countKnifePyrrole = 0
    for i in 1:particles
        push!(trajectories[2], [])
        # TODO: Validate results
        for u in dataPyrrole[i].u
            if (u.z ≥ knifeZ-0.001 && u.z ≤ knifeZ+0.001)
                averageKnifePyrrole += u.y
                countKnifePyrrole += 1
            end
            if (u.z ≥ destination-0.001 && u.z ≤ destination+0.005)
                finalCountPyrrole += 1
                averageDeflectionPyrrole += u.y
                push!(histogramData[2], u.y)
                break
            end
            push!(trajectories[2][i], (u.z, u.y))
        end
    end
    averageDeflectionPyrrole /= finalCountPyrrole;
    averageKnifePyrrole /= countKnifePyrrole;
    print("-> average knife pyrrole: ", averageKnifePyrrole, "\n")
    print("-> average deflection: ", averageDeflectionPyrrole, ", total: ", finalCountPyrrole, "\n");
end
range=-0.002:0.0001:0.002
plot(histogramData[1], seriestype=:histogram, bins=range, fillalpha=0.5, labels=["Pyrrole Water", "Pyrrole"][1])
plot!(histogramData[2], seriestype=:histogram, bins=range, fillalpha=0.5, labels=["Pyrrole Water", "Pyrrole"][2])
#plot([[pair[1] for pair ∈ trajectory] for trajectory ∈ trajectories[1]], [[pair[2] for pair ∈ trajectory] for trajectory ∈ trajectories[1]], legend=false, title="Pyrrole Water")
#plot([[pair[1] for pair ∈ trajectory] for trajectory ∈ trajectories[2]], [[pair[2] for pair ∈ trajectory] for trajectory ∈ trajectories[2]], legend=false, title="Pyrrole")
#plot([trajectory[1] for trajectory ∈ trajectories[1][1]], [trajectory[2] for trajectory ∈ trajectories[1][1]], label = "Pyrrole Water")
#plot!([trajectory[1] for trajectory ∈ trajectories[2][1]], [trajectory[2] for trajectory ∈ trajectories[2][1]], label = "Pyrrole")

