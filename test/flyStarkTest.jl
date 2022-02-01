using Plots
############## Plotting stuff (for debugging)
# Note that in order to use this stuff, other stuff further down has to be commented in
histogramData = [[], []]
trajectories = [[[]], [[]]]

############## Constants
original0 = 0.34572
skimmer1Z     = -0.28072 + original0
skimmer2Z     = -0.044 + original0
deflectorZ    = 0 + original0
# 154mm as Sebastian stated
deflectorEndZ = 0.154 + original0
deflectorX = -0.00198571
deflectorEndX = deflectorX + (140-1) * 2.85714e-5
deflectorY = -0.00398571
deflectorEndY = deflectorY + (210-1) * 2.85714e-5
knifeZ        = 0.16528 + original0
knifeY = 0
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

u = 1.66053906660e-27
# Tolerance
ε = 0.01

# Use only 2 in z direction as there's no change along z
zCount = 2
# In kV
recordedVoltage = 60.0
paperVoltage = 14.0

particles = 10000
############## Distributions
truncatedNormal(mean, sigma, limit) = CMInject.Truncated(CMInject.Normal(mean, sigma), mean-limit, mean+limit)
distsPyrrole = Dict(:x => truncatedNormal(0, sigmaPos, limitPos),
                    :y => truncatedNormal(0, sigmaPos, limitPos),
                    :z => CMInject.Dirac(0),
                    :vx => truncatedNormal(0, sigmaVtrans, limitVtrans),
                    :vy => truncatedNormal(0, sigmaVtrans, limitVtrans),
                    :vz => truncatedNormal(meanVz, sigmaVlong, limitVlong),
                    # C4H5N
                    :m => CMInject.Dirac((4*12+5*1.0078+14.003)*u))
# Copy all values...
distsPyrroleWater = Dict(distsPyrrole)
# ... except for the mass
distsPyrroleWater[:m] = CMInject.Dirac((4*12+7*1.0078+14.003+15.995)*u)

# TODO: Eventually I might want to allow other states too
stateDists = Dict(:J => CMInject.DiscreteUniform(0,0), :M => CMInject.DiscreteUniform(0,0))

############## Simulation
@testset "Fly Stark Simulation" begin

    sourcePyrroleWater = CMInject.StarkSamplingSource{CMInject.StarkParticle{Float64}, Float64}(
                                     distsPyrroleWater, stateDists, "test/pyrrole-water.molecule")
    sourcePyrrole = CMInject.StarkSamplingSource{CMInject.StarkParticle{Float64}, Float64}(
                                     distsPyrrole, stateDists, "test/pyrrole.molecule")

    @test CMInject.generate(sourcePyrroleWater, 1) != nothing

    field = CMInject.ElectricField("test/example_field.h5")
    @inline function detectHits(args...)
        x = args[1].x
        y = args[1].y
        z = args[1].z
        hitSkimmer(x, y, z, sz, sr) = z ≥ sz-ε && z ≤ sz && sqrt(x*x + y*y) ≥ sr
        if (hitSkimmer(x, y, z, skimmer1Z, skimmer1R))
            return true
        end
        if (hitSkimmer(x, y, z, skimmer2Z, skimmer2R))
            return true
        end
        if (hitSkimmer(x, y, z, skimmer3Z, skimmer3R))
            return true
        end
        if (z ≥ knifeZ-ε && z ≤ knifeZ &&
            y ≤ knifeY)
            return true
        end
        if (z ≥ deflectorZ && z ≤ deflectorEndZ &&
            (x ≤ deflectorX || x ≥ deflectorEndX ||
             y ≤ deflectorY || y ≥ deflectorEndY))
            return true
        end
        # No need to continue simulation when particle has arrived
        if (z ≥ destination)
            return true
        end
        false
    end
    solverOpts=(adaptive=false,
                dense=false,
                callback=CMInject.DiscreteCallback(detectHits, integrator -> CMInject.terminate!(integrator)))
    getExperiment(source) = CMInject.Experiment(;source=source,
                                            n_particles=particles,
                                            fields=(field,),
                                            # TODO: Figure out how to use detectors
                                            detectors=Tuple(CMInject.SectionDetector{Float64,:y}.(
                                                             -0.001:0.001:0.001,
                                                             true)),
                                            # TODO: Try with Runge-Kutta 4/5
                                            solver=CMInject.EulerHeun(),
                                            time_span=(0.0, 2e-3),
                                            time_step=1e-6,
                                            ensemble_alg=CMInject.EnsembleThreads(),
                                            solver_opts=solverOpts)
    experimentPyrroleWater = getExperiment(sourcePyrroleWater)
    experimentPyrrole = getExperiment(sourcePyrrole)

    simResPyrroleWater = CMInject.simulate(experimentPyrroleWater)
    simResPyrrole = CMInject.simulate(experimentPyrrole)
    dataPyrroleWater = simResPyrroleWater[1]
    dataPyrrole = simResPyrrole[1]

    @test simResPyrroleWater != nothing
    @test simResPyrrole != nothing
    @test dataPyrroleWater != nothing
    @test dataPyrrole != nothing

    function getAverageDeflection(index, data)
        averageDeflection = 0
        finalCount = 0
        for i in 1:particles
            # Comment in if wanting to plot trajectories
            # push!(trajectories[index], [])
            for u in data[i].u
                if (u.z ≥ destination-0.001 && u.z ≤ destination+0.005)
                    finalCount += 1
                    averageDeflection += u.y
                    # Comment in if wanting to plot histograms
                    #push!(histogramData[index], u.y)
                    break
                end
                # Comment in if wanting to plot trajectories
                # push!(trajectories[index][i], (u.z, u.y))
            end
        end
        averageDeflection /= finalCount
    end
    averageDeflectionPyrroleWater = getAverageDeflection(1, dataPyrroleWater)
    averageDeflectionPyrrole = getAverageDeflection(2, dataPyrrole)

    expectedPyrroleWater = -6.27e-4
    expectedPyrrole = -5.92e-5
    @test abs(averageDeflectionPyrroleWater - expectedPyrroleWater) ≤ -expectedPyrroleWater
    @test abs(averageDeflectionPyrrole - expectedPyrrole) ≤ -expectedPyrrole
end
# Comment in to plot histogram
#range=-0.002:0.0001:0.002
#plot(histogramData[1], seriestype=:histogram, bins=range, fillalpha=0.5, labels=["Pyrrole Water", "Pyrrole"][1])
#plot!(histogramData[2], seriestype=:histogram, bins=range, fillalpha=0.5, labels=["Pyrrole Water", "Pyrrole"][2])
# Comment in to plot trajectories
#plot([[pair[1] for pair ∈ trajectory] for trajectory ∈ trajectories[1]], [[pair[2] for pair ∈ trajectory] for trajectory ∈ trajectories[1]], legend=false, title="Pyrrole Water")
#plot([[pair[1] for pair ∈ trajectory] for trajectory ∈ trajectories[2]], [[pair[2] for pair ∈ trajectory] for trajectory ∈ trajectories[2]], legend=false, title="Pyrrole")

