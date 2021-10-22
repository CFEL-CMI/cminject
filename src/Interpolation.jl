using Logging
using HDF5
using LabelledArrays
using Statistics
using Interpolations
import Interpolations: scale as itpscale


function _reduce_to_labelled_array(A::AbstractArray{T, N}, InnerT::Type{U}) where {T, N, U <: SLArray}
    sz = Base.front(size(A))
    out_A = Array{U, N - 1}(undef, sz)
    for i in CartesianIndices(sz)
        out_A[i] = InnerT(A[i, :])
    end
    return out_A
end


function _apply_default_symbols(symbols, filename)
	# If symbols were created from pure indices, attempt
	if all(symbols .== Symbol.(0:length(symbols)-1))
		defaults = nothing
		if length(symbols) == 3
			defaults = (:vx, :vz, :p)
		elseif length(symbols) == 4
			defaults = (:vx, :vy, :vz, :p)
		end

		if isnothing(defaults)
			@warn "No names were defined in your file $(filename). Could not assign default names, "
			       "so using $(symbols) - it would be better if the file contained names for the "
				   "columns!"
			return symbols
		else
			@warn "No names were defined in your file $(filename). Implicitly using names $(symbols) "
			       "- it would be better if the file contained names for the columns!"
			return defaults
		end
	else
		symbols
	end
end

function _grid_points_to_range(grid_points)
	a = grid_points[1]
	z = grid_points[end]
	# try to reduce local numerical errors by taking the mean of all neighbor differences
	spacing = mean(diff(grid_points))
	a:spacing:z
end

function hdf5_to_interpolator(filename::AbstractString, names=:default)
    h5open(filename) do h5f
		# Read the index (grid points) from the file
		idx = h5f["index"]
		index = [read(idx, k) for k in keys(idx)]
		index = Tuple(index)
		sizes = (size(idx)[1] for idx in reverse(index))

		# Read the field data from the file
		h5data = h5f["data"]
		N = length(keys(h5data))
		if names == :default
			names = _apply_default_symbols(Tuple(Symbol.(keys(h5data))), filename)
		end
		# Construct a LabelledArray type based on the names in the file / default names
		InnerT = @SLVector Float64 names
		data_ = cat([
			permutedims(reshape(read(h5data, k), sizes...), (2,1))
			for k in keys(h5data)
		]..., dims=3)
		# Convert the array 
        data = _reduce_to_labelled_array(data_, InnerT)

		grid_ranges = _grid_points_to_range.(index)
        itp = itpscale(interpolate(data, BSpline(Linear())), grid_ranges...)
        extrapolate(itp, InnerT(repeat([NaN], N)))
    end
end

function interpolateStarkCurve(ΔE, E_min,
        energies::Vector{T})::AbstractInterpolation where T<:Real
    # TODO: Validate that quadratic is reasonable (e.g. vs. cubic / linear)
    itp = interpolate(energies, BSpline(Quadratic(Natural(OnGrid()))))
    scaledItp = itpscale(itp, range(E_min; length=length(energies), step=ΔE))
    extrapolate(scaledItp, NaN)
end
