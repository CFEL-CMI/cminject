using HDF5

"""
    main(args)

`args` should contain the following elements (in the same order):
1. norm file name
2. gradient file name
3. destination HDF5 file name
4. start z in m
5. end z in m
6. scaling factor
Both input files should be in the format that was used by CMIfly.

Example call: julia --project=. ./scripts/FlyToHDF5.jl ./test/example_field ./test/example_gradient ./test/example_field.h5 0.34572 0.49972 0.23333333333333
"""
function main(args)
    expectedArgumentCount = 6
    if (size(args)[1] != expectedArgumentCount)
        error("Expected $expectedArgumentCount arguments, got ", size(args)[1], "\n")
    end

    fieldFile = args[1]
    gradientFile = args[2]
    destinationFile = args[3]
    startZ = parse(Float64, args[4])
    endZ = parse(Float64, args[5])
    scaling = parse(Float64, args[6])

    (normGrid, gradGrid, initial, stepSizes) = getGrids(fieldFile, gradientFile, startZ, endZ, scaling)
    saveToHDF5(normGrid, gradGrid, initial, stepSizes, destinationFile)
end

function getGrids(fieldFile, gradientFile, startZ, endZ, scaling)
    ############## Common stuff for the norm and the gradient
    function applyZ(initial, stepSizes, counts)
        initial[3] = startZ
        stepSizes[3] = endZ - startZ
        # Use only 2 in z direction as there's no change along z
        counts[3] = 2
    end
    getInitial(fieldLines) = [parse(Float64, token) for token ∈ split(fieldLines[1], " ") if length(token) > 0]
    getStepSizes(fieldLines) = [parse(Float64, token) for token ∈ split(fieldLines[2], " ") if length(token) > 0]
    getCounts(fieldLines) = [parse(Int, token) for token ∈ split(fieldLines[3], " ") if length(token) > 0]
    parseGrid(counts, parseFunction) = [parseFunction(x, y, z)
                                        # y (later to be x) is stored in reverse order
                                        for y ∈ (counts[2]-1):-1:0, x ∈ 0:(counts[1]-1), z ∈ 0:(counts[3]-1)]

    ############## Norm field

    fieldLines = readlines(fieldFile)
    initial = getInitial(fieldLines)
    stepSizes = getStepSizes(fieldLines)
    counts = getCounts(fieldLines)

    applyZ(initial, stepSizes, counts)

    normGrid = parseGrid(counts, (x,y,z) -> parse(Float64, fieldLines[4 + y + x * counts[2]]) * scaling)

    ############## Gradient field

    gradFieldLines = readlines(gradientFile)
    gradInitial = getInitial(gradFieldLines)
    gradStepSizes = getStepSizes(gradFieldLines)
    gradCounts = getCounts(gradFieldLines)

    applyZ(gradInitial, gradStepSizes, gradCounts)

    gradGrid = [parseGrid(gradCounts, (x,y,z) ->
                          parse(Float64, split(gradFieldLines[4 + y + x * gradCounts[2]], " ", keepempty=false)[i])
                          * scaling) for i ∈ 1:3]

    function swap(a, i, j)
        tmp = a[i]
        a[i] = a[j]
        a[j] = tmp
    end
    # Swap x and y
    swap(gradGrid, 1, 2)

    if (initial != gradInitial || stepSizes != gradStepSizes)
        error("Norm and gradient bases aren't equal!")
    end

    # Swap x and y for initial and step sizes too
    swap(initial, 1, 2)
    swap(stepSizes, 1, 2)
    
    (normGrid, gradGrid, initial, stepSizes)
end

function saveToHDF5(normGrid, gradGrid, initial, stepSizes, destination)
    h5open(destination, "w") do file
        fieldGroup = create_group(file, "field")
        fieldGroup["norm"] = normGrid
        fieldGroup["gradient"] = [gradGrid[i][x,y,z] for x ∈ 1:size(gradGrid[1])[1],
                                  y ∈ 1:size(gradGrid[1])[2],
                                  z ∈ 1:size(gradGrid[1])[3],
                                  i ∈ 1:size(gradGrid)[1]]
        metaGroup = create_group(file, "metadata")
        attributes(metaGroup)["description"] = "This field was automatically converted from a field type as is used in the CMIfly project"
        attributes(metaGroup)["dimension"] = 3
        attributes(metaGroup)["initial"] = initial
        attributes(metaGroup)["steps"] = stepSizes
    end
end

main(ARGS)
