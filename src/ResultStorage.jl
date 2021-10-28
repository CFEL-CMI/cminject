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
        HDF5.h5t_insert(dtype, fieldname(T,i), fieldoffset(T,i), datatype(fieldtype(T,i)))
    end
    HDF5.Datatype(dtype)
end


function store_properties(file::HDF5.File, particles::AbstractVector{PT}) where {PT<:AbstractParticle}
    # Collect properties
    props_dtype = get_hdf5_datatype(associated_props_type(PT))
    props = [ntfromstruct(particle._p) for particle in particles]

    props_dset = create_dataset(file, "particles/properties", props_dtype, dataspace(props))
    write_dataset(props_dset, props_dtype, props)
end

function store_trajectories(file::HDF5.File, solution)
    # Collect trajectories from `solution[i].u` and `solution[i].t` for every `i`
    trajs = [NamedTuple.(solution[i].u) for i ∈ 1:length(solution)]
    times = [[(t=t_,) for t_ in solution[i].t] for i ∈ 1:length(solution)]
    # Concatenate them together as one matrix that will be stored
    # TODO: in the future we'll need to play around with variable lengths, or put each trajectory
    # in its own dataset: not all trajectories may be equally long!
    trajs_with_times = hcat([merge.(traj, time) for (traj, time) in zip(trajs, times)]...)

    # Only store trajectories if there is at least one
    if !isempty(trajs_with_times)
        trajs_with_times_dtype = get_hdf5_datatype(eltype(trajs_with_times))
        trajs_dset = create_dataset(
            file, "particles/trajectories",
            trajs_with_times_dtype, dataspace(trajs_with_times)
        )
        write_dataset(trajs_dset, trajs_with_times_dtype, trajs_with_times)
    end
end

function store_results(
    stor::HDF5ResultStorage{PT}, solution, hits, particles::AbstractVector{PT},
) where {PT<:AbstractParticle}
    h5open(stor.path, "cw") do f
        store_properties(f, particles)
        store_trajectories(f, solution)
        # TODO: store the hits as well
    end

    if !isempty(hits)
        @warn "Detector hits are currently not stored in HDF5 files - sorry!"
    end
end
