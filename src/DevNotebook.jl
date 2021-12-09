### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 8712d4d9-6edf-4531-990b-e5f698937690
using Pkg; cd("/home/chip/Documents"); Pkg.activate("CMInject"); Pkg.status()

# ╔═╡ de326458-4d1e-44ed-b6af-cf69221e3e1a
using HDF5, Interpolations

# ╔═╡ ec4eccbf-cc53-4dcd-8320-67c420253b08
using PlutoUI, ImageView, Plots

# ╔═╡ a25196c5-97de-452a-9e82-10f2dd1c85b0
using CUDA

# ╔═╡ 70b66fdc-592d-42aa-801d-a799b7f47a86
using DifferentialEquations, DifferentialEquations.EnsembleAnalysis

# ╔═╡ 54f9eece-cc60-4bb7-ae2c-3e74669e46b4
using DiffEqGPU

# ╔═╡ d08a7856-40f1-490f-ae83-a08b990a5293
using BenchmarkTools

# ╔═╡ 0e428388-2360-46fb-9606-39f1885ea954
import Distributions: Uniform

# ╔═╡ 0c68f428-f741-497f-849e-dd8d3f841da4
plotly()

# ╔═╡ db03ba1b-4d09-485e-a88c-36f7d9851299
struct FluidPoint
	vx::Float64
	vz::Float64
	p::Float64
end

# ╔═╡ 4ffe9ddf-1ef9-47a2-b923-c129a4b5f6b8
isbitstype(FluidPoint)

# ╔═╡ a65d2c89-0fde-41a0-be24-02eba4e10380
begin
	Base.:+(a::FluidPoint, b::FluidPoint) = FluidPoint(a.vx+b.vx, a.vz+b.vz, a.p+b.p)
	Base.:*(a::FluidPoint, b::FluidPoint) = FluidPoint(a.vx*b.vx, a.vz*b.vz, a.p*b.p)

	Base.:*(a::FluidPoint, λ::Number) = FluidPoint(a.vx*λ, a.vz*λ, a.p*λ)
	Base.:*(λ::Number, a::FluidPoint) = a * λ
end

# ╔═╡ 45a9b688-1f64-437e-ae14-22bab58efee3
replaceNaNs = true

# ╔═╡ 0c9a84ce-5048-4ca2-8c64-b48363c972fd
nodes, data = h5open("/home/chip/Documents/cminject-materials/goldspheres-27nm/nitrogen_1.8mbar_extended_nan.h5") do h5f
	idx = h5f["index"]

	index = [read(idx, k) for k in keys(idx)]
	index = tuple(index...)
	sizes = (size(idx)[1] for idx in reverse(index))

	dat = h5f["data"]
	data = [permutedims(reshape(read(dat, k), sizes...), (2,1)) for k in keys(dat)]

	if replaceNaNs
		for dataset in data
			dataset[isnan.(dataset)] .= 0.0
		end
	end

	index, data
end;

# ╔═╡ a9ed88c7-2e72-4ccb-ad5e-7980dd9f083e
begin
	points = FluidPoint.(data...)
	points[100, 750]
end

# ╔═╡ 3108b39e-7acb-4f46-8d55-0a2baa476a28
vr, vz, p = data;

# ╔═╡ 4810691d-7878-465e-869b-58a5e849ba44
heatmap(nodes[1], nodes[2], p')

# ╔═╡ 2d7d6df4-aa19-498d-b09c-4071485cf5e1
grid = [[x, y] for x in nodes[1] for y in nodes[2]]

# ╔═╡ 5342159f-4b64-4df0-8f97-fe7f49b16509
begin
	r0 = nodes[1][1]
	r1 = nodes[1][end]
	dr = 1e-4
	r0, r1, dr
end

# ╔═╡ 4983bf76-1f46-4008-a8ad-7a8eda0389e1
begin
	z0 = nodes[2][1]
	z1 = nodes[2][end]
	dz = 1e-4
	z0, z1, dz
end

# ╔═╡ 7b5a55ff-cbf0-4d76-ac77-37366c112030
length(z0:dz:z1)

# ╔═╡ e48d86ed-bed9-4ed6-9a4f-6f10fb301632
Gray.(nodes[1] .== collect(r0:dr:r1))

# ╔═╡ 608679d3-e09f-41a7-84d1-54847e27e7db
begin
	rs = -0.01:0.000002:0.01
	zs = -0.1285:0.000012:0.005
end

# ╔═╡ 28a8a290-5e95-47c0-bb3a-4f040f0e0047
v = hcat(map(t -> [t[1], t[2]], collect(zip(rs, zs))))

# ╔═╡ 678095be-c515-4464-a977-17710ed5e3ab
outside_point = FluidPoint(NaN, NaN, NaN)

# ╔═╡ df5ae5d1-affd-4897-be6b-1c755d9cff4d
itp = extrapolate(interpolate(nodes, points, Gridded(Linear())), outside_point)

# ╔═╡ 3c6d05e9-3431-43fb-a676-6ee9a8e18868
begin
	RS = rand(rs[1]:1e-10:rs[end], 10000)
	ZS = rand(zs[1]:1e-10:zs[end], 10000)
	@btime IP = itp.(RS, ZS)
end

# ╔═╡ 90f47568-4058-4a33-bc56-60bc5609d80a
kB = 1.380649e-23

# ╔═╡ e65b92a3-7ec8-42c2-8e48-003b0b117491
T = 293.15

# ╔═╡ a22ed30a-b1fa-4615-a0da-32d42e608cc5
μ=1.76e-5

# ╔═╡ 7fa6d680-88c3-4205-a06c-d25a7c3da785
mf=4.7e-26

# ╔═╡ 1c124bba-24b8-4d9a-a307-a8bb77881aa0
Knf = μ * √(π * kB * T / (2mf))

# ╔═╡ ec23bcc3-7ed3-4566-a5f7-1f503fc4dc51
function stokes!(ddu, du, u, params, t)
	x, z = u
	vx, vz = du

	ρ, rp = params
	mp = 4/3 * rp^3. * π * ρ

	fluid = itp(x, z)
	Δvx = fluid.vx - vx
	Δvz = fluid.vz - vz

	Kn = Knf / (fluid.p * rp)
	Cc = 1 + Kn * (1.231 + 0.4695 * exp(-1.1783 / Kn))

	a₀ = 6π * μ * rp / Cc / mp
	ax = a₀ * Δvx
	az = a₀ * Δvz

	ddu[1] = ax
	ddu[2] = az
	nothing
end

# ╔═╡ 8cba414e-4eef-49ad-929d-79a4171770a1
params = (1050, 13.5e-9)

# ╔═╡ 775a0e63-92dc-421a-9f68-eb428332aaa4
let
	tspan = (0.0, 1e-3)
	prob = SecondOrderODEProblem{true}(stokes!, [0.0, 10.0], [0.005, -0.128], tspan, params)
	sol = solve(prob, Tsit5(), dt=1e-5, reltol=1e-8, abstol=1e-8)
	plot(sol)
end

# ╔═╡ aa4879b1-c77f-49ce-acc9-b3358dedf4f3
function stokesdirect!(du, u, params, t)
	x, vx, z, vz = u
	ρ, rp = params
	mp = 4/3 * rp^3. * π * ρ

	fluid = itp(x, z)
	Δvx = fluid.vx - vx
	Δvz = fluid.vz - vz

	Kn = Knf / (fluid.p * rp)
	Cc = 1 + Kn * (1.231 + 0.4695 * exp(-1.1783 / Kn))

	a₀ = 6π * μ * rp / Cc / mp
	ax = a₀ * Δvx
	az = a₀ * Δvz

	du[1] = u[2]
	du[2] = ax
	du[3] = u[4]
	du[4] = az
	nothing
end

# ╔═╡ 4e01f492-0037-4d92-b7f2-4bd0ffa2a2e1
let
	tspan = (0.0, 0.1)
	prob = ODEProblem{true}(stokesdirect!, [0.000, 10.0, -0.128, 0.0], tspan, params)
	sol = solve(prob, Tsit5(), dt=1e-5, reltol=1e-8, abstol=1e-8)
	plot(sol[1,:])
end

# ╔═╡ 94a7761d-3507-457b-8b4e-2b259fbe5b0b
function ahmadidirect!(du, u, params, t)
	x, vx, z, vz = u

	ρ, rp = params
	s₀fπ = π * 216μ * kB * T / (π^2 * (2rp)^5 * ρ^2)
	mp = 4/3 * rp^3. * π * ρ

	fluid = itp(x, z)
	Δvx = fluid.vx - vx
	Δvz = fluid.vz - vz

	Kn = Knf / (fluid.p * rp)
	Cc = 1 + Kn * (1.231 + 0.4695 * exp(-1.1783 / Kn))
	s₀ = √(s₀fπ / Cc)

	du[1] = 0
	du[2] = s₀
	du[3] = 0
	du[4] = s₀
	nothing
end

# ╔═╡ 2a8e116f-6537-4878-9422-d9f2111097b7
@bind rp2 Slider(10e-9:10e-9:100e-9, show_value=true)

# ╔═╡ 8899754e-f126-4b95-9dd1-fbfb522a119c
@bind rho2 Slider(1000:1000:20000, show_value=true, default=19320)

# ╔═╡ 78a816d5-8856-4f0e-a41c-e0d3f24ec9ec
params2 = [rho2, rp2]

# ╔═╡ 157d1492-4502-422f-95a9-abf831a3bb33
begin
	tspan = (0.0, 0.03)
	prob = SDEProblem{true}(stokesdirect!, ahmadidirect!,
		[0.0, 0.0, -0.128, 10.0], tspan, params2)
	ensembleprob = EnsembleProblem(prob, safetycopy=false)

	@time sol = solve(ensembleprob, EulerHeun(),
		trajectories=1000, dt=1e-5, unstable_check=((_...) -> false))
	println("doneski")
	nothing
end

# ╔═╡ 46bbd1eb-d2af-447f-80b3-4870cd8b63e3
begin
	rxs = [solu[1,:] for solu in sol]
	rzs = [solu[3,:] for solu in sol]
	plot(rxs, rzs)
end

# ╔═╡ 8cd737fa-c269-47d1-be86-0cac1fe9761a


# ╔═╡ 2c63f8b2-1e2b-4854-af89-01131f2cb4f7


# ╔═╡ 6d3a16be-754c-4dd6-883d-ed9de5f2ef0b


# ╔═╡ 375dfc43-e92f-4c49-8309-89f0ff7a1914


# ╔═╡ e8e9ba5f-1dd7-48c8-a569-3a19c730f369


# ╔═╡ 2019a71a-7b1c-4342-8dd7-8a71a101f8ee


# ╔═╡ 254686e2-3968-469f-87b8-c84734d19c21


# ╔═╡ a5c202ec-b975-4bc2-9cdd-8dd2b4354daf


# ╔═╡ 24936751-c85b-470a-b6df-f15cdee06120
function nodes_to_range(nodes::AbstractArray{T,N}; delta=nothing) where {T, N}
	xα = nodes[1]
	xβ = nodes[2]
	xω = nodes[end]
	if isnothing(delta)
		Δx = xβ - xα
	else
		Δx = delta
	end
	StepRangeLen{T, T, T}(xα, Δx, length(nodes))
end

# ╔═╡ d44d0a6c-ae33-4bb0-82d6-b8d260803931
ranges = nodes_to_range.(nodes; delta=1e-3)

# ╔═╡ f6c281ba-ec80-437a-97e9-78807fc4b8cc
rng = ranges[1]

# ╔═╡ 8b8c23b6-6bf9-4273-b89f-dd631e80f87d
rng.step

# ╔═╡ 390f47dd-4d58-4995-9458-4b32304b0eff
data[1]

# ╔═╡ cce8afcc-4d02-4ef7-979f-32f0521e0ccb
findall(≠(0), data[1])

# ╔═╡ dda0047d-9751-478f-b39c-9ff0ddec4470
StepRange

# ╔═╡ dba6b812-115e-48d9-ac21-f8b548e22c3f
# D = spatial dimensions, N = interpolated quantities
struct RegularGrid{D,N,T,A <: AbstractArray{T, D}}
	ranges::NTuple{D, StepRangeLen{T, T, T}}
	data::NTuple{N, A}
end

# ╔═╡ 0868fbe6-e096-422b-a010-7fa4ea71bc8d
rgrid = RegularGrid(ranges, Tuple(data))

# ╔═╡ bdde88b0-3097-4e62-9b33-875761bf64ec
typeof(rgrid)

# ╔═╡ 927dfd7f-1005-4ba3-861d-a33f2951644a
function interpol!(out::M, grid::RegularGrid{2,N,T}, xs::V, ys::V) where {N,T,M <: AbstractMatrix{T}, V <: AbstractVector{T}}
	xrange = grid.ranges[1]
	yrange = grid.ranges[2]

	for (k, (x, y)) in enumerate(zip(xs, ys))
		scaledx = (x - xrange[1]) ./ xrange.step
		scaledy = (y - yrange[1]) ./ yrange.step
		xw, xi = modf(scaledx)
		yw, yi = modf(scaledy)
		xi = Int(xi)
		yi = Int(yi)

		outofbounds = (xi < 1) || (xi >= xrange.len) || (yi < 1) || (yi >= yrange.len)

		if outofbounds
			for i in 1:N
				out[k,i] = NaN
			end
		else
			for i in 1:N
				dat = grid.data[i]'
				a0 = dat[xi,yi  ]*xw + dat[xi+1,yi  ]*(1-xw)
				a1 = dat[xi,yi+1]*xw + dat[xi+1,yi+1]*(1-xw)
				out[k,i] = a0*yw + a1*(1-yw)
			end
		end
	end

	nothing
end

# ╔═╡ e2c3fac2-dfd0-4cc7-b8ad-f8ca714134ba
function interpol!(out::V, grid::RegularGrid{2,N,T}, x::T, y::T)	where {N,T,V <: AbstractVector{T}}
	xrange = grid.ranges[1]
	yrange = grid.ranges[2]

	scaledx = (x - xrange[1]) ./ xrange.step
	scaledy = (y - yrange[1]) ./ yrange.step
	xw, xi = modf(scaledx)
	yw, yi = modf(scaledy)
	xi = Int(xi)
	yi = Int(yi)

	outofbounds = (xi < 1) || (xi >= xrange.len) || (yi < 1) || (yi >= yrange.len)

	if outofbounds
		for i in 1:N
			out[i] = NaN
		end
	else
		for i in 1:N
			dat = grid.data[i]'
			a0 = dat[xi,yi  ]*xw + dat[xi+1,yi  ]*(1-xw)
			a1 = dat[xi,yi+1]*xw + dat[xi+1,yi+1]*(1-xw)
			out[i] = a0*yw + a1*(1-yw)
		end
	end

	nothing
end

# ╔═╡ 669f1e32-b291-41ac-a1c7-cce7239121ef
function interpol(grid::RegularGrid{2,N,T}, x::T, y::T) where {N,T}
	out = similar(1:N, T)
	interpol!(out, grid, x, y)
	out
end

# ╔═╡ 214c4aed-ddda-4f38-a3e2-a64f16d3f00a
function to_interpol(grid::RegularGrid{D,N,T}) where {D,N,T}
	(x::T, y::T) -> interpol(grid, x, y)
end

# ╔═╡ 2d8f79bf-5020-4d60-8c2e-f5573b5067db
ip = to_interpol(rgrid); typeof(ip)

# ╔═╡ 7e941aae-4286-4bc2-93b0-1823c0556b6f
begin
	N = 10000
	randx = rand(Uniform(-0.01, 0.01), N)
	randy = rand(Uniform(-0.1285, 0.0), N)
end

# ╔═╡ 86a8ae8f-ec96-4932-bfc2-0cc71addf659
@btime ip.(randx, randy)

# ╔═╡ 1f4fedb2-9109-4a6b-a6a5-4ce6dda120a0
function to_reused_interpol(grid::RegularGrid{D,N,T}, M::Integer) where {D,N,T}
	matrix = Matrix{T}(undef, M, N)
	(xs::Vector{T}, ys::Vector{T}) -> begin
		interpol!(matrix, grid, xs, ys)
		matrix
	end
end

# ╔═╡ b8950bdd-0a79-4fb1-89b2-8754c61d09aa
ip2 = to_reused_interpol(rgrid, 10000)

# ╔═╡ Cell order:
# ╠═8712d4d9-6edf-4531-990b-e5f698937690
# ╠═de326458-4d1e-44ed-b6af-cf69221e3e1a
# ╠═ec4eccbf-cc53-4dcd-8320-67c420253b08
# ╠═0e428388-2360-46fb-9606-39f1885ea954
# ╠═a25196c5-97de-452a-9e82-10f2dd1c85b0
# ╠═70b66fdc-592d-42aa-801d-a799b7f47a86
# ╠═54f9eece-cc60-4bb7-ae2c-3e74669e46b4
# ╠═0c68f428-f741-497f-849e-dd8d3f841da4
# ╠═d08a7856-40f1-490f-ae83-a08b990a5293
# ╠═db03ba1b-4d09-485e-a88c-36f7d9851299
# ╠═4ffe9ddf-1ef9-47a2-b923-c129a4b5f6b8
# ╠═a65d2c89-0fde-41a0-be24-02eba4e10380
# ╠═45a9b688-1f64-437e-ae14-22bab58efee3
# ╠═0c9a84ce-5048-4ca2-8c64-b48363c972fd
# ╠═a9ed88c7-2e72-4ccb-ad5e-7980dd9f083e
# ╠═3108b39e-7acb-4f46-8d55-0a2baa476a28
# ╠═4810691d-7878-465e-869b-58a5e849ba44
# ╠═2d7d6df4-aa19-498d-b09c-4071485cf5e1
# ╠═5342159f-4b64-4df0-8f97-fe7f49b16509
# ╠═4983bf76-1f46-4008-a8ad-7a8eda0389e1
# ╠═7b5a55ff-cbf0-4d76-ac77-37366c112030
# ╠═e48d86ed-bed9-4ed6-9a4f-6f10fb301632
# ╠═608679d3-e09f-41a7-84d1-54847e27e7db
# ╠═28a8a290-5e95-47c0-bb3a-4f040f0e0047
# ╠═678095be-c515-4464-a977-17710ed5e3ab
# ╠═df5ae5d1-affd-4897-be6b-1c755d9cff4d
# ╠═3c6d05e9-3431-43fb-a676-6ee9a8e18868
# ╠═90f47568-4058-4a33-bc56-60bc5609d80a
# ╠═e65b92a3-7ec8-42c2-8e48-003b0b117491
# ╠═a22ed30a-b1fa-4615-a0da-32d42e608cc5
# ╠═7fa6d680-88c3-4205-a06c-d25a7c3da785
# ╠═1c124bba-24b8-4d9a-a307-a8bb77881aa0
# ╠═ec23bcc3-7ed3-4566-a5f7-1f503fc4dc51
# ╠═8cba414e-4eef-49ad-929d-79a4171770a1
# ╠═775a0e63-92dc-421a-9f68-eb428332aaa4
# ╠═aa4879b1-c77f-49ce-acc9-b3358dedf4f3
# ╠═4e01f492-0037-4d92-b7f2-4bd0ffa2a2e1
# ╠═94a7761d-3507-457b-8b4e-2b259fbe5b0b
# ╠═2a8e116f-6537-4878-9422-d9f2111097b7
# ╠═8899754e-f126-4b95-9dd1-fbfb522a119c
# ╠═78a816d5-8856-4f0e-a41c-e0d3f24ec9ec
# ╠═157d1492-4502-422f-95a9-abf831a3bb33
# ╠═46bbd1eb-d2af-447f-80b3-4870cd8b63e3
# ╠═8cd737fa-c269-47d1-be86-0cac1fe9761a
# ╠═2c63f8b2-1e2b-4854-af89-01131f2cb4f7
# ╠═6d3a16be-754c-4dd6-883d-ed9de5f2ef0b
# ╠═375dfc43-e92f-4c49-8309-89f0ff7a1914
# ╠═e8e9ba5f-1dd7-48c8-a569-3a19c730f369
# ╠═2019a71a-7b1c-4342-8dd7-8a71a101f8ee
# ╠═254686e2-3968-469f-87b8-c84734d19c21
# ╠═a5c202ec-b975-4bc2-9cdd-8dd2b4354daf
# ╠═24936751-c85b-470a-b6df-f15cdee06120
# ╠═d44d0a6c-ae33-4bb0-82d6-b8d260803931
# ╠═f6c281ba-ec80-437a-97e9-78807fc4b8cc
# ╠═8b8c23b6-6bf9-4273-b89f-dd631e80f87d
# ╠═390f47dd-4d58-4995-9458-4b32304b0eff
# ╠═cce8afcc-4d02-4ef7-979f-32f0521e0ccb
# ╠═dda0047d-9751-478f-b39c-9ff0ddec4470
# ╠═dba6b812-115e-48d9-ac21-f8b548e22c3f
# ╠═0868fbe6-e096-422b-a010-7fa4ea71bc8d
# ╠═bdde88b0-3097-4e62-9b33-875761bf64ec
# ╠═927dfd7f-1005-4ba3-861d-a33f2951644a
# ╠═e2c3fac2-dfd0-4cc7-b8ad-f8ca714134ba
# ╠═669f1e32-b291-41ac-a1c7-cce7239121ef
# ╠═214c4aed-ddda-4f38-a3e2-a64f16d3f00a
# ╠═2d8f79bf-5020-4d60-8c2e-f5573b5067db
# ╠═7e941aae-4286-4bc2-93b0-1823c0556b6f
# ╠═86a8ae8f-ec96-4932-bfc2-0cc71addf659
# ╠═1f4fedb2-9109-4a6b-a6a5-4ce6dda120a0
# ╠═b8950bdd-0a79-4fb1-89b2-8754c61d09aa
