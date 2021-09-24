using HDF5
using LabelledArrays
using Interpolations

struct RegularGrid{D,N<:Any,T,A<:AbstractArray{T,D}}
    ranges::NTuple{D,StepRangeLen{T,T,T}}
    data::NTuple{N,A}
end

function interpol1!(out::V, grid::RegularGrid{2,N,T}, x::T, y::T) where {N, T, V<:AbstractVector{T}}
    xrange = grid.ranges[1]
    yrange = grid.ranges[2]

    if isnan(x) || isnan(y)
        out .= NaN
        return nothing
    end

    scaledx = (x - xrange[1]) / xrange.step
    scaledy = (y - yrange[1]) / yrange.step
    xw, xif = modf(scaledx)
    yw, yif = modf(scaledy)
    xi = Int(xif)+1  # TODO: +1 correct? seems good at first glance, also in 2D, but hmm...
    yi = Int(yif)+1

    outofbounds = (xi < 1) || (xi >= xrange.len) || (yi < 1) || (yi >= yrange.len)

    @inbounds begin
        if outofbounds
            out .= NaN
        else
            for i in eachindex(grid.data)
                dat = grid.data[i]'
                a0 = dat[xi,yi  ]*(1-xw) + dat[xi+1,yi  ]*xw
                a1 = dat[xi,yi+1]*(1-xw) + dat[xi+1,yi+1]*xw
                out[i] = a0*(1-yw) + a1*yw
            end
        end
    end

    nothing
end

function interpol1(grid::RegularGrid{2,N,T,A}, x::T, y::T) where {T, N, A}
	out = zeros(T, N)
	interpol1!(out, grid, x, y)
	out
end

function nodes_to_range(nodes::AbstractArray{T,N}; delta=nothing) where {T, N}
    xα = nodes[1]
    xβ = nodes[2]
	Δx = isnothing(delta) ? xβ - xα : delta
    StepRangeLen{T, T, T}(xα, Δx, length(nodes))
end

function hdf5_to_regulargrid(filename::AbstractString)
    h5open(filename) do h5f
		idx = h5f["index"]
		index = [read(idx, k) for k in keys(idx)]
		index = tuple(index...)
		sizes = (size(idx)[1] for idx in reverse(index))
        ranges = nodes_to_range.(index)

		dat = h5f["data"]
		data = [
			permutedims(reshape(read(dat, k), sizes...), (2,1))'
			for k in keys(dat)
		]
        RegularGrid(ranges, Tuple(data))
    end
end


GridT = @SLVector Float64 (:vx, :vz, :p)
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

        itp = scale(interpolate(data_, BSpline(Linear())), -0.01:1e-4:0.01, -0.1285:1e-4:0.005)
        extrapolate(itp, GridT(NaN, NaN, NaN))
    end
end