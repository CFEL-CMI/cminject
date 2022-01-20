using Plots
histogramData = [[], []]
# TODO: Reduce duplicates resulting from Pyrrole AND Pyrrole Water

# TODO: There are contradictoray information for these values
skimmer1Z = 0.065
skimmer2Z = 0.302
# This value comes from CMIfly and is not available in the paper
deflectorZ= 0.437
deflectorEndZ = skimmer2Z + 0.197
knifeZ    = deflectorEndZ + 0.015
skimmer3Z = knifeZ + 0.025
skimmer1R = 0.003
skimmer2R = 0.0015
skimmer3R = 0.0015
destination = skimmer3Z + 0.176

# FIELD START

fieldLines=readlines("test/example_field")
initial = [parse(Float64, token) for token ∈ split(fieldLines[1], " ") if length(token) > 0]
stepSizes = [parse(Float64, token) for token ∈ split(fieldLines[2], " ") if length(token) > 0]
counts = [parse(Int, token) for token ∈ split(fieldLines[3], " ") if length(token) > 0]

initial[3] = deflectorZ
stepSizes[3] = deflectorEndZ - deflectorZ
counts[3] = 2

myGrid = [parse(Float64, fieldLines[4 + y + x * counts[2]]) for y = (counts[2]-1):-1:0, x = (counts[1]-1):-1:0, z = 0:(counts[3]-1)]
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

gradGrid = [[map(x->parse(Float64, x), split(gradFieldLines[4 + y + x * gradCounts[2]], " ", keepempty=false))[i]
             for y = (gradCounts[2]-1):-1:0, x = (gradCounts[1]-1):-1:0, z = 0:(gradCounts[3]-1)] for i = 1:3]
# x and y are wrong :/
tmp = gradGrid[1]
gradGrid[1] = gradGrid[2]
gradGrid[2] = (-1) .* tmp
gradItps = @. CMInject.interpolate(gradGrid, CMInject.BSpline(CMInject.Linear()))
gradExts = @. CMInject.extrapolate(gradItps, 0)
gradItpScaleds = tuple([CMInject.itpscale(gradExts[i], 
                              initial[2]:stepSizes[2]:(initial[2]+(counts[2]-1)*stepSizes[2]),
                              initial[1]:stepSizes[1]:(initial[1]+(counts[1]-1)*stepSizes[1]),
                              initial[3]:stepSizes[3]:deflectorEndZ) for i = 1:3]...)


# GRADIENT END

@testset "Fly Stark Simulation" begin
    # TODO: Check if calculated gradient is much different to given one
    # TODO: Read & use given gradient

    distsPyrrole = Dict(:x => CMInject.Normal(0, 0.0001*0.0001),
                        :y => CMInject.Normal(0, 0.0001*0.0001),
                        :z => CMInject.Dirac(0),
                        :vx => CMInject.Normal(0, 1.5*1.5),
                        :vy => CMInject.Normal(0, 1.5*1.5),
                        :vz => CMInject.Normal(670, 7*7),
                        # C4H5N
                        :m => CMInject.Dirac((4*12+5*1.0078+14.003)/9.223e18))
    distsPyrroleWater = Dict(:x => CMInject.Normal(0, 0.0001*0.0001),
                             :y => CMInject.Normal(0, 0.0001*0.0001),
                             :z => CMInject.Dirac(0),
                             :vx => CMInject.Normal(0, 1.5*1.5),
                             :vy => CMInject.Normal(0, 1.5*1.5),
                             :vz => CMInject.Normal(670, 7*7),
                             # C4H7NO
                             :m => CMInject.Dirac((4*12+7*1.0078+14.003+15.995)/9.223e18))
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
        x = args[2].x
        y = args[2].y
        z = args[2].z
        if (z ≥ skimmer1Z-0.01 && z ≤ skimmer1Z &&
            sqrt(x*x + y*y) ≥ skimmer1R)
            print("HIT skimmer 1\n")
            return true
        end
        if (z ≥ skimmer2Z-0.01 && z ≤ skimmer2Z &&
            sqrt(x*x + y*y) ≥ skimmer2R)
            print("HIT skimmer 2\n")
            return true
        end
        if (z ≥ skimmer3Z-0.01 && z ≤ skimmer3Z &&
            sqrt(x*x + y*y) ≥ skimmer3R)
            print("HIT skimmer 3\n")
            return true
        end
        if (z ≥ knifeZ-0.01 && z ≤ knifeZ &&
            y ≤ knifeY)
            print("HIT knife\n")
            return true
        end
        if (z ≥ deflectorZ && z ≤ deflectorEndZ &&
            (x ≤ initial[1] || x ≥ initial[1]+(counts[1]-1)*stepSizes[1] ||
             y ≤ initial[2] || y ≥ initial[2]+(counts[2]-1)*stepSizes[2]))
            print("HIT deflector\n")
            return true
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
                                            # TODO: Try with Runge-Kutta 4/5
                                            solver=CMInject.EulerHeun(),
                                            # TODO: If necessary, adjust the time span
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

    knifeY = 3.188356743548741e-6
    simResPyrroleWater = CMInject.simulate(experimentPyrroleWater)
    knifeY = 1.5357245757532834e-6
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
        # TODO: Validate results
        for u in dataPyrroleWater[i].u
            if (u.z ≥ knifeZ-0.001 && u.z ≤ knifeZ+0.001)
                averageKnifePyrroleWater += u.y
                countKnifePyrroleWater += 1
            end
            if (u.z ≥ destination-0.001 && u.z ≤ destination+0.005)
                finalCountPyrroleWater += 1
                averageDeflectionPyrroleWater += u.y
                #print(u.x, ", ", u.y, ", ", u.z, "\n")
                push!(histogramData[1], u.y)
                break
            end
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
        # TODO: Validate results
        for u in dataPyrrole[i].u
            if (u.z ≥ knifeZ-0.001 && u.z ≤ knifeZ+0.001)
                averageKnifePyrrole += u.y
                countKnifePyrrole += 1
            end
            if (u.z ≥ destination-0.001 && u.z ≤ destination+0.005)
                finalCountPyrrole += 1
                averageDeflectionPyrrole += u.y
                #print(u.x, ", ", u.y, ", ", u.z, "\n")
                if (size(histogramData[2]) < size(histogramData[1]))
                    push!(histogramData[2], u.y)
                end
                break
            end
        end
    end
    averageDeflectionPyrrole /= finalCountPyrrole;
    averageKnifePyrrole /= countKnifePyrrole;
    print("-> average knife pyrrole: ", averageKnifePyrrole, "\n")
    print("-> average deflection: ", averageDeflectionPyrrole, ", total: ", finalCountPyrrole, "\n");
end
plot(histogramData[2], seriestype=:histogram, nbins=10, fillalpha=0.5, labels=["Pyrrole Water", "Pyrrole"][2])
plot!(histogramData[1], seriestype=:histogram, nbins=10, fillalpha=0.5, labels=["Pyrrole Water", "Pyrrole"][1])

