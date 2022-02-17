using HDF5
using LabelledArrays
using NamedTupleTools  # Since HDF5.jl seems to need NamedTuples for writing to compound datatypes...


abstract type AbstractResultStorage end

struct HDF5ResultStorage{PT<:AbstractParticle} <: AbstractResultStorage
    path::String
end


function get_hdf5_datatype(::Type{T}) where T
    dtype = HDF5.h5t_create(HDF5.H5T_COMPOUND, sizeof(T))
    @inbounds for i in 1:length(T.types)
        if (applicable(datatype, fieldtype(T,i)))
            HDF5.h5t_insert(dtype, fieldname(T,i), fieldoffset(T,i), datatype(fieldtype(T,i)))
        end
    end
    HDF5.Datatype(dtype)
end


function store_properties!(file::HDF5.File, particles::AbstractVector{PT}) where {PT<:AbstractParticle}
    # Collect properties
    dtype = get_hdf5_datatype(associated_props_type(PT))
    props = [ntfromstruct(particle._p) for particle in particles]

    props_dset = create_dataset(file, "particles/properties", dtype, dataspace(props))
    write_dataset(props_dset, dtype, props)
end

"""
    pad(array, fill, length)

Pads the array with fill to the specified length
"""
function pad(array::V, fill::T, length::N) where {N<:Int, T, V<:AbstractVector{T}}
    vcat(array, repeat([fill], length-size(array)[1]))
end

function store_trajectories!(file::HDF5.File, solution)
    # Collect trajectories from `solution[i].u` and `solution[i].t` for every `i`
    # FIXME: this is terribly slow because of type instability from the NamedTuple call - we need to somehow force
    #        HDF5 to accept our structs instead of converting to a NamedTuple...
    #        I've seen this said to be working on some GitHub issue comment for HDF5.jl, but I just got random floats
    #        when I tried passing an array of structs rather than one of NamedTuples :/
    trajs = map(sol -> NamedTuple.(sol.u), solution)
    times = [[(t=t_,) for t_ in solution[i].t] for i ∈ 1:length(solution)]

    # For now, pad the trajectories with NaN to take different trajectory sizes into account
    allTrajectories = [merge.(traj, time) for (traj, time) in zip(trajs, times)]
    maxTrajectoryLength = maximum(size.(allTrajectories))[1]
    nanTuple = NamedTuple([(key, NaN) for key ∈ keys(allTrajectories[1][1])])
    paddedTrajectories = [pad(trajectory, nanTuple, maxTrajectoryLength) for trajectory ∈ allTrajectories]

    # Concatenate them together as one matrix that will be stored
    trajs_with_times = hcat(paddedTrajectories...)

    # Only store trajectories if there is at least one
    if !isempty(trajs)
        dtype = get_hdf5_datatype(eltype(trajs_with_times))
        trajs_dset = create_dataset(
            file, "particles/tracked/trajectories",
            dtype, dataspace(trajs)
        )
        write_dataset(trajs_dset, dtype, trajs)
    end

    nothing
end

function store_hits!(file::HDF5.File, detector_hits)
    for (i, detector) in enumerate(detector_hits)
        if !isempty(detector)
            hits = NamedTuple.(detector)
            dtype = get_hdf5_datatype(eltype(hits))
            hits_dset = create_dataset(
                file, "detectors/$(i)", dtype, dataspace(hits)
            )
            write_dataset(hits_dset, dtype, hits)
        end
    end

    nothing
end

function store_results!(
    stor::HDF5ResultStorage{PT}, solution, hits, particles::AbstractVector{PT};
    store_trajectories=true
) where {PT<:AbstractParticle}
    h5open(stor.path, "cw") do f
        store_properties!(f, particles)
        store_hits!(f, hits)
        if store_trajectories
            store_trajectories!(f, solution)
        end
    end

    nothing
end
