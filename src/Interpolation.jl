using HDF5
using LabelledArrays
import Interpolations: scale as itpscale, interpolate, extrapolate, BSpline,
                       Linear, Quadratic, AbstractInterpolation, OnGrid, Natural, gradient


GridT = @SLVector Float64 (:vx, :vz, :p)  # TODO generalize: read from file!
function hdf5_to_interpolator(filename::AbstractString)
    h5open(filename) do h5f
		idx = h5f["index"]
		index = [read(idx, k) for k in keys(idx)]
		index = tuple(index...)
		sizes = (size(idx)[1] for idx in reverse(index))

		dat = h5f["data"]
		data = cat([
			permutedims(reshape(read(dat, k), sizes...), (2,1))'
			for k in keys(dat)
		]..., dims=3)
        data_ = [
            GridT(data[k, l, :]) for l in 1:size(data)[2], k in 1:size(data)[1]
        ]

        # TODO generalize scale: read from file!
        itp = itpscale(interpolate(data_, BSpline(Linear())), -0.01:1e-4:0.01, -0.1285:1e-4:0.005)
        extrapolate(itp, GridT(NaN, NaN, NaN))
    end
end

function interpolateStarkCurve(ΔE, E_min,
        energies::Vector{T})::AbstractInterpolation where T<:Real
    # TODO: Validate that quadratic is reasonable (e.g. vs. cubic / linear)
    itp = interpolate(energies, BSpline(Quadratic(Natural(OnGrid()))))
    scaledItp = itpscale(itp, range(E_min; length=length(energies), step=ΔE))
    extrapolate(scaledItp, NaN)
end
